"""Unit tests for beckmann_core.geometry -- pure (sym,x,y,z)-tuple math, no
Gaussian-log dependency."""
import math

from beckmann_core.geometry import displace_leaving_group, find_leaving_group, no_distance

# A minimal synthetic "protonated oxime tail": N(1) - O(2) with two H's(3,4)
# bonded to O, plus one unrelated atom(5) far away. 1-based indices throughout,
# matching beckmann_core.geometry's convention.
ATOMS = [
    ("N", 0.0, 0.0, 0.0),
    ("O", 1.40, 0.0, 0.0),
    ("H", 1.70, 0.90, 0.0),
    ("H", 1.70, -0.90, 0.0),
    ("C", -1.50, 0.0, 0.0),
]
NI, OI = 1, 2


def test_no_distance_matches_euclidean():
    assert math.isclose(no_distance(ATOMS, NI, OI), 1.40, rel_tol=1e-9)


def test_find_leaving_group_includes_o_and_its_hydrogens_only():
    group = find_leaving_group(ATOMS, OI)
    assert group == {2, 3, 4}  # O + its two H's, not N or the far carbon


def test_displace_leaving_group_moves_only_leaving_group():
    delta = 0.30
    displaced = displace_leaving_group(ATOMS, NI, OI, delta)
    # N and the unrelated carbon are untouched
    assert displaced[0] == ATOMS[0]
    assert displaced[4] == ATOMS[4]
    # O and its H's moved
    assert displaced[1] != ATOMS[1]
    assert displaced[2] != ATOMS[2]
    assert displaced[3] != ATOMS[3]


def test_displace_leaving_group_increases_no_distance_by_exactly_delta():
    delta = 0.30
    displaced = displace_leaving_group(ATOMS, NI, OI, delta)
    new_r = no_distance(displaced, NI, OI)
    assert math.isclose(new_r, 1.40 + delta, rel_tol=1e-9)


def test_displace_leaving_group_preserves_leaving_group_internal_geometry():
    """The two H's should keep their distance to O unchanged -- a rigid
    translation, not a distortion."""
    delta = 0.30
    displaced = displace_leaving_group(ATOMS, NI, OI, delta)
    o_before, h1_before = ATOMS[1], ATOMS[2]
    o_after, h1_after = displaced[1], displaced[2]
    dist_before = math.sqrt(sum((a - b) ** 2 for a, b in zip(o_before[1:], h1_before[1:])))
    dist_after = math.sqrt(sum((a - b) ** 2 for a, b in zip(o_after[1:], h1_after[1:])))
    assert math.isclose(dist_before, dist_after, rel_tol=1e-9)
