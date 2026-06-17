"""Unit tests for _spatial.py geometric primitives."""

from __future__ import annotations

import numpy as np
import pytest

from gbsa_pipeline._spatial import _find_clashing_residues, contact_pairs


def test_contact_pairs_returns_empty_when_no_atoms() -> None:
    """Empty coordinate arrays produce no pairs."""
    assert contact_pairs(np.empty((0, 3)), np.array([[1.0, 0.0, 0.0]]), 2.0) == []
    assert contact_pairs(np.array([[1.0, 0.0, 0.0]]), np.empty((0, 3)), 2.0) == []


def test_contact_pairs_detects_close_pair() -> None:
    """Atoms within cutoff are returned with their distance."""
    a = np.array([[0.0, 0.0, 0.0]])
    b = np.array([[1.0, 0.0, 0.0]])
    result = contact_pairs(a, b, 2.0)
    assert len(result) == 1
    i, j, dist = result[0]
    assert i == 0
    assert j == 0
    assert abs(dist - 1.0) < 1e-6


def test_contact_pairs_boundary_inclusive() -> None:
    """A pair exactly at the cutoff distance is included."""
    a = np.array([[0.0, 0.0, 0.0]])
    b = np.array([[2.0, 0.0, 0.0]])
    result = contact_pairs(a, b, 2.0)
    assert len(result) == 1


def test_contact_pairs_excludes_distant_pair() -> None:
    """Atoms beyond cutoff are not returned."""
    a = np.array([[0.0, 0.0, 0.0]])
    b = np.array([[10.0, 0.0, 0.0]])
    assert contact_pairs(a, b, 2.0) == []


def test_contact_pairs_multiple_pairs() -> None:
    """All pairs within cutoff are returned; those outside are not."""
    a = np.array([[0.0, 0.0, 0.0], [5.0, 0.0, 0.0]])
    b = np.array([[1.0, 0.0, 0.0], [10.0, 0.0, 0.0]])
    result = contact_pairs(a, b, 2.0)
    assert len(result) == 1
    assert result[0][:2] == (0, 0)


def test_find_clashing_residues_detects_clash() -> None:
    """A candidate within cutoff of a reference atom is returned."""
    candidates = [("WAT1", (1.0, 0.0, 0.0))]
    reference = [(0.0, 0.0, 0.0)]
    result = _find_clashing_residues(candidates, reference, cutoff=1.5)
    assert result == {"WAT1"}


def test_find_clashing_residues_no_clash() -> None:
    """A candidate beyond cutoff is not returned."""
    candidates = [("WAT1", (10.0, 0.0, 0.0))]
    reference = [(0.0, 0.0, 0.0)]
    result = _find_clashing_residues(candidates, reference, cutoff=1.5)
    assert result == set()


def test_find_clashing_residues_boundary_inclusive() -> None:
    """A candidate exactly at the cutoff distance is included."""
    candidates = [("WAT1", (1.5, 0.0, 0.0))]
    reference = [(0.0, 0.0, 0.0)]
    result = _find_clashing_residues(candidates, reference, cutoff=1.5)
    assert result == {"WAT1"}


def test_find_clashing_residues_multiple_candidates_mixed() -> None:
    """Only candidates within cutoff of any reference atom are returned."""
    candidates = [
        ("WAT1", (0.5, 0.0, 0.0)),  # close
        ("WAT2", (10.0, 0.0, 0.0)),  # far
        ("WAT3", (0.0, 0.5, 0.0)),  # close to a different reference atom
    ]
    reference = [(0.0, 0.0, 0.0), (0.0, 1.0, 0.0)]
    result: set[str] = _find_clashing_residues(candidates, reference, cutoff=1.0)
    assert result == {"WAT1", "WAT3"}


def test_find_clashing_residues_empty_candidates() -> None:
    """Empty candidate list returns an empty set."""
    result: set[str] = _find_clashing_residues([], [(0.0, 0.0, 0.0)], cutoff=1.0)
    assert result == set()


def test_find_clashing_residues_empty_reference() -> None:
    """Empty reference list means nothing clashes."""
    candidates = [("WAT1", (0.0, 0.0, 0.0))]
    result = _find_clashing_residues(candidates, [], cutoff=1.0)
    assert result == set()


def test_find_clashing_residues_invalid_cutoff() -> None:
    """Non-positive cutoff raises ValueError."""
    with pytest.raises(ValueError, match="cutoff must be positive"):
        _find_clashing_residues([("WAT1", (0.0, 0.0, 0.0))], [(0.0, 0.0, 0.0)], cutoff=0.0)
