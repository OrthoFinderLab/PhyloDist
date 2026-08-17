from __future__ import annotations

import numpy as np

from .tree_structure import TreeStructure


VALID_FROBENIUS_TYPES = frozenset({"raw", "relative"})


def _leaf_distance_matrix(tree) -> tuple[np.ndarray, list[str]]:
    """Return ETE4's triangular distance matrix as a symmetric array."""
    names = list(tree.leaf_names())
    triangular = tree.distance_matrix()
    size = len(names)

    if len(triangular) != size:
        raise ValueError("ETE4 returned an invalid leaf-distance matrix.")

    matrix = np.zeros((size, size), dtype=float)
    for row_index, row in enumerate(triangular):
        values = np.asarray(row, dtype=float)
        expected_size = size - row_index - 1
        if values.size != expected_size:
            raise ValueError("ETE4 returned an invalid leaf-distance matrix.")

        matrix[row_index, row_index + 1:] = values
        matrix[row_index + 1:, row_index] = values

    return matrix, names


def frobenius_distance(
    tree1: TreeStructure,
    tree2: TreeStructure,
    distance_type: str = "raw",
) -> float:
    """Compare tree path lengths using their leaf-distance matrices.

    Only the upper triangle is used, so every unordered leaf pair contributes
    once.  For ``distance_type="relative"``, tree2 is the reference and the
    raw distance is divided by the Euclidean norm of tree2's pair distances.
    """
    distance_type = distance_type.lower().strip()
    if distance_type not in VALID_FROBENIUS_TYPES:
        raise ValueError(
            "frobenius_type must be one of {'raw', 'relative'}"
        )

    if tree1.t is None or tree2.t is None:
        raise ValueError(
            "Frobenius distance requires trees parsed by ETE4 because the "
            "fallback parser does not retain branch lengths."
        )

    matrix1, names1 = _leaf_distance_matrix(tree1.t)
    matrix2, names2 = _leaf_distance_matrix(tree2.t)

    names1 = list(names1)
    names2 = list(names2)
    if len(names1) != len(set(names1)) or len(names2) != len(set(names2)):
        raise ValueError("Frobenius distance requires unique leaf names.")

    taxa1 = set(names1)
    taxa2 = set(names2)
    if taxa1 != taxa2:
        raise ValueError(
            "Trees do not contain the same leaves. "
            f"Missing from tree1: {sorted(taxa2 - taxa1)}; "
            f"missing from tree2: {sorted(taxa1 - taxa2)}"
        )

    index2 = {name: i for i, name in enumerate(names2)}
    order2 = [index2[name] for name in names1]
    matrix2 = matrix2[np.ix_(order2, order2)]

    upper = np.triu_indices_from(matrix1, k=1)
    distances1 = matrix1[upper]
    distances2 = matrix2[upper]
    raw = float(np.linalg.norm(distances1 - distances2))

    if distance_type == "raw":
        return raw

    reference_norm = float(np.linalg.norm(distances2))
    if reference_norm == 0.0:
        return float("nan")
    return raw / reference_norm
