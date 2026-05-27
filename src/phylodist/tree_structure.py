from __future__ import annotations

import os
import re
from typing import FrozenSet, Iterable, List, Optional

from ete4 import Tree
from . import utils

Taxon = str
Part = FrozenSet[Taxon]


class TreeStructure:
    """
    Thin wrapper around an ete4.Tree that exposes rooted clades or unrooted splits.

    Robust behaviour:
        1. Try ETE4 first.
        2. If ETE4 fails, fall back to a simple topology-only Newick parser.

    Accepted input sources:
        - Newick string
        - path to a file containing Newick
        - path to a NEXUS file containing a tree
        - file-like object with .read()
        - ete4.Tree
    """

    @staticmethod
    def _extract_rooted_flag(text: str) -> tuple[str, Optional[bool]]:
        rooted = None
        stripped = text.lstrip()

        if stripped.startswith("[&R]"):
            rooted = True
            stripped = re.sub(r"^\[&R\]\s*", "", stripped)
        elif stripped.startswith("[&U]"):
            rooted = False
            stripped = re.sub(r"^\[&U\]\s*", "", stripped)

        return stripped, rooted

    def __init__(
        self,
        source: str | Tree,
        is_rooted: Optional[bool] = None,
        parser: int | None = None,
    ):
        self.source = source
        self.t = None
        self._fallback_root = None
        self._fallback_taxa: FrozenSet[Taxon] | None = None

        inferred_rooted = None

        if isinstance(source, Tree):
            self.t = source.copy()
            self.taxa: FrozenSet[Taxon] = frozenset(self.t.leaf_names())
            self.is_rooted: Optional[bool] = is_rooted
            return

        text = self._coerce_tree_text(source)
        text, inferred_rooted = self._extract_rooted_flag(text)
        text = self._clean_tree_text(text)

        ete_errors = []

        parsers = [parser] if parser is not None else [1, 0, 2, 3]

        for p in parsers:
            if p is None:
                continue

            try:
                self.t = Tree(text, parser=p)
                self.taxa = frozenset(self.t.leaf_names())
                self.is_rooted = is_rooted if is_rooted is not None else inferred_rooted
                return
            except Exception as e:
                ete_errors.append(f"ETE parser={p}: {e}")

        # ------------------------------------------------------------
        # Fallback: parse topology directly.
        # This avoids ETE4 choking on support labels / long branch lengths.
        # ------------------------------------------------------------
        try:
            topology_text = self._to_topology_only_newick(text)
            self._fallback_root = self._parse_topology_newick(topology_text)
            self._fallback_taxa = frozenset(self._collect_fallback_taxa(self._fallback_root))
            self.taxa = self._fallback_taxa
            self.is_rooted = is_rooted if is_rooted is not None else inferred_rooted
            return

        except Exception as e:
            msg = (
                "Could not parse tree with ETE4 or phylodist fallback parser.\n"
                f"Original source: {source}\n"
                f"Fallback error: {e}\n"
                "ETE errors:\n  "
                + "\n  ".join(ete_errors)
            )
            raise ValueError(msg)

    @classmethod
    def from_newick(
        cls,
        newick: str,
        is_rooted: Optional[bool] = None,
        parser: int | None = None,
    ) -> "TreeStructure":
        return cls(newick, is_rooted=is_rooted, parser=parser)

    @classmethod
    def from_tree(
        cls,
        tree: Tree,
        is_rooted: Optional[bool] = None,
        parser: int | None = None,
    ) -> "TreeStructure":
        return cls(tree, is_rooted=is_rooted, parser=parser)

    @staticmethod
    def _coerce_tree_text(source) -> str:
        if hasattr(source, "read"):
            return source.read()

        text = str(source)

        if os.path.exists(text):
            with open(text, "r", encoding="utf-8") as f:
                return f.read()

        return text

    @staticmethod
    def _clean_tree_text(text: str) -> str:
        """
        Normalize tree text into one Newick string.
        """
        text = text.strip().lstrip("\ufeff").strip()

        # NEXUS detection.
        if text.lower().startswith("#nexus") or "begin trees" in text.lower():
            match = re.search(r"=\s*(?:\[&[^\]]+\]\s*)?(.+?);", text, re.S | re.I)
            if not match:
                raise ValueError("Could not find tree in NEXUS content.")
            text = match.group(1).strip()

        # If non-NEXUS but contains tree assignment:
        # tree xxx = [&R] (...);
        elif "=" in text and "(" in text:
            match = re.search(r"=\s*(?:\[&[^\]]+\]\s*)?(.+?);", text, re.S | re.I)
            if match:
                text = match.group(1).strip()

        # Remove rootedness annotation such as [&R] / [&U].
        text = re.sub(r"^\[&[^\]]+\]\s*", "", text)

        # Remove comments anywhere in the tree.
        text = re.sub(r"\[.*?\]", "", text)

        # Remove whitespace/newlines.
        text = re.sub(r"\s+", "", text)

        if not text.endswith(";"):
            text += ";"

        return text

    @staticmethod
    def _to_topology_only_newick(text: str) -> str:
        """
        Convert Newick to topology-only Newick.

        Examples:
            ((A:0.1,B:0.2)100:0.3,C:0.4);
        becomes:
            ((A,B),C);

        This is safe for topology-only distances.
        """
        text = text.strip()

        # Remove branch lengths, including scientific notation.
        text = re.sub(r":[^,();]+", "", text)

        # Remove internal node labels/supports after closing parentheses.
        # After branch length removal, )100:0.1 has become )100.
        text = re.sub(r"\)([^,();]+)(?=[,);])", ")", text)

        if not text.endswith(";"):
            text += ";"

        return text

    @staticmethod
    def _parse_topology_newick(text: str):
        """
        Parse a topology-only Newick string into nested tuples.

        Leaf:
            "A"

        Internal node:
            ("A", "B", ("C", "D"))
        """
        text = text.strip()
        if text.endswith(";"):
            text = text[:-1]

        i = 0
        n = len(text)

        def parse_name():
            nonlocal i

            if i < n and text[i] == "'":
                i += 1
                start = i
                while i < n and text[i] != "'":
                    i += 1
                name = text[start:i]
                if i < n and text[i] == "'":
                    i += 1
                return name

            start = i
            while i < n and text[i] not in ",()":
                i += 1

            return text[start:i].strip()

        def parse_subtree():
            nonlocal i

            if i >= n:
                raise ValueError("Unexpected end of Newick string")

            if text[i] == "(":
                i += 1
                children = []

                while True:
                    child = parse_subtree()
                    children.append(child)

                    if i >= n:
                        raise ValueError("Unexpected end inside internal node")

                    if text[i] == ",":
                        i += 1
                        continue

                    if text[i] == ")":
                        i += 1
                        break

                    raise ValueError(f"Unexpected character {text[i]!r} at position {i}")

                # Ignore any internal node label that survived.
                while i < n and text[i] not in ",()":
                    i += 1

                return tuple(children)

            name = parse_name()
            if not name:
                raise ValueError(f"Empty leaf name at position {i}")

            return name

        root = parse_subtree()

        if i != n:
            raise ValueError(f"Unexpected trailing text at position {i}: {text[i:i+50]!r}")

        return root

    @staticmethod
    def _collect_fallback_taxa(node) -> set[Taxon]:
        if isinstance(node, str):
            return {node}

        taxa = set()
        for child in node:
            taxa.update(TreeStructure._collect_fallback_taxa(child))

        return taxa

    @staticmethod
    def _fallback_internal_parts(node, *, include_root=False):
        """
        Return descendant leaf sets for every internal node.
        """
        out = []

        def rec(x, is_root):
            if isinstance(x, str):
                return frozenset([x])

            leaves = set()
            for child in x:
                leaves.update(rec(child, False))

            fp = frozenset(leaves)

            if include_root or not is_root:
                out.append(fp)

            return fp

        rec(node, True)
        return out

    def __repr__(self) -> str:
        kind = (
            "rooted" if self.is_rooted is True
            else "unrooted" if self.is_rooted is False
            else "unknown-rootedness"
        )
        source = "ete4" if self.t is not None else "fallback"
        return f"TreeStructure(n_taxa={len(self.taxa)}, {kind}, source={source})"

    def __len__(self) -> int:
        return self.n_taxa

    def __iter__(self):
        if self.t is None:
            return iter(())
        return self.t.traverse()

    def write(self, *args, **kwargs):
        if self.t is None:
            raise ValueError("Cannot write fallback-parsed tree as ete4 tree.")
        return self.t.write(*args, **kwargs)

    @property
    def n_taxa(self) -> int:
        return len(self.taxa)

    def get_tree_taxa(self) -> FrozenSet[Taxon]:
        return self.taxa

    def copy(self) -> "TreeStructure":
        if self.t is not None:
            return TreeStructure(self.t.copy(), is_rooted=self.is_rooted)

        # Rebuild from fallback topology string is unnecessary for comparison use.
        new = object.__new__(TreeStructure)
        new.source = self.source
        new.t = None
        new._fallback_root = self._fallback_root
        new._fallback_taxa = self._fallback_taxa
        new.taxa = self.taxa
        new.is_rooted = self.is_rooted
        return new

    def _iter_internal_nodes(self):
        if self.t is None:
            return iter(())

        for node in self.t.traverse():
            if not node.is_leaf:
                yield node

    def get_rooted_clades(self) -> List[Part]:
        if self.is_rooted is False:
            raise ValueError("Cannot extract rooted clades from a tree marked unrooted.")

        if self.t is not None:
            clades: List[Part] = []
            for node in self._iter_internal_nodes():
                leaves = frozenset(node.leaf_names())
                if 1 < len(leaves) < self.n_taxa - 1:
                    clades.append(leaves)
            return clades

        clades = []
        for leaves in self._fallback_internal_parts(self._fallback_root, include_root=False):
            if 1 < len(leaves) < self.n_taxa - 1:
                clades.append(leaves)

        return clades

    def get_unrooted_splits(self) -> List[Part]:
        if self.is_rooted is True:
            raise ValueError("Cannot extract unrooted splits from a tree marked rooted.")

        if self.t is not None:
            splits: set[Part] = set()
            for node in self._iter_internal_nodes():
                if node.is_root:
                    continue

                leaves = frozenset(node.leaf_names())
                canon = utils._canon_unrooted_split(leaves, self.taxa)
                if canon is not None:
                    splits.add(canon)

            return list(splits)

        splits: set[Part] = set()

        for leaves in self._fallback_internal_parts(self._fallback_root, include_root=False):
            canon = utils._canon_unrooted_split(leaves, self.taxa)
            if canon is not None:
                splits.add(canon)

        return list(splits)

    def normalize_parts(
        self,
        parts: Iterable[Iterable[Taxon]],
        *,
        rooted: bool,
    ) -> List[Part]:
        out: List[Part] = []
        seen: set[Part] = set()

        for p in parts:
            fp = utils._to_fset(p)
            canon = (
                utils._canon_rooted_clade(fp, self.taxa)
                if rooted else
                utils._canon_unrooted_split(fp, self.taxa)
            )
            if canon is not None and canon not in seen:
                seen.add(canon)
                out.append(canon)

        return out

    @staticmethod
    def project_parts_to_shared_taxa(
        parts: Iterable[Part],
        old_taxa: Iterable[Taxon],
        shared_taxa: Iterable[Taxon],
        rooted: bool,
    ) -> List[Part]:
        old_taxa = utils._to_fset(old_taxa)
        shared_taxa = utils._to_fset(shared_taxa)
        del old_taxa

        projected: List[Part] = []
        seen: set[Part] = set()

        for p in parts:
            p = utils._to_fset(p)
            q = p & shared_taxa
            canon = (
                utils._canon_rooted_clade(q, shared_taxa)
                if rooted else
                utils._canon_unrooted_split(q, shared_taxa)
            )
            if canon is not None and canon not in seen:
                seen.add(canon)
                projected.append(canon)

        return projected