import numpy as np
import pytest

from phylodist import PhyloDist, Tree


TREE1 = "((A:1,B:2)AB:3,(C:1,D:1)CD:2);"
TREE2 = "((D:1,C:2)CD:3,(B:1,A:1)AB:2);"


def test_raw_frobenius_uses_unique_leaf_pairs_and_aligns_names():
    result = PhyloDist.from_newick(TREE1, TREE2, parser=1).compare(
        method="frobenius"
    )
    assert result.score == pytest.approx(np.sqrt(4.0))
    assert result.shared_taxa == frozenset({"A", "B", "C", "D"})


def test_relative_frobenius_uses_tree2_as_reference():
    comparison = PhyloDist.from_newick(TREE1, TREE2, parser=1)
    raw = comparison.compare(method="frobenius").score
    relative = comparison.compare(
        method="frobenius", frobenius_type="relative"
    ).score
    expected_reference_norm = np.linalg.norm([2, 8, 7, 8, 7, 3])
    assert relative == pytest.approx(raw / expected_reference_norm)


def test_frobenius_rejects_different_taxa():
    comparison = PhyloDist(
        Tree("((A:1,B:1):1,(C:1,D:1):1);"),
        Tree("((A:1,B:1):1,(C:1,E:1):1);"),
    )
    with pytest.raises(ValueError, match="same leaves"):
        comparison.compare(method="frobenius")


def test_frobenius_rejects_normalise_option():
    comparison = PhyloDist.from_newick(TREE1, TREE2, parser=1)
    with pytest.raises(ValueError, match="not defined"):
        comparison.compare(method="frobenius", normalise=True)
