"""Layer 1: slots of w, group enumeration, and the signed-permutation action, checked against the notes."""
import pytest

from symembed import GroupSpec
from symembed.slots import build_slots, iter_group_elements, signed_permutation_table, slot_orbits


def images(spec, k):
    """For each group element, g.w as a list of signed slot names."""
    slots = build_slots(spec, k)
    elements = tuple(iter_group_elements(spec))
    index, sign = signed_permutation_table(slots, elements)
    names = [s.name(spec.labels) for s in slots]
    return names, {g.cycle_notation(spec.labels): [("-" if sign[r][j] < 0 else "") + names[index[r][j]]
                                                   for j in range(len(slots))]
                   for r, g in enumerate(elements)}


def test_worked_example_slot_order():
    names, _ = images(GroupSpec((("a", "b"), ("p",))), 2)
    assert names == ["d_aa", "d_bb", "d_pp", "d_ab", "d_ap", "d_bp", "e_ab", "e_ap", "e_bp"]


def test_worked_example_orbits():
    spec = GroupSpec((("a", "b"), ("p",)))
    slots = build_slots(spec, 2)
    assert slot_orbits(spec, slots) == ((0, 1), (2,), (3,), (4, 5), (6,), (7, 8))


def test_worked_example_without_parity():
    _, gw = images(GroupSpec((("a", "b"), ("p",))), 2)
    assert gw == {
        "e": ["d_aa", "d_bb", "d_pp", "d_ab", "d_ap", "d_bp", "e_ab", "e_ap", "e_bp"],
        "(a b)": ["d_bb", "d_aa", "d_pp", "d_ab", "d_bp", "d_ap", "-e_ab", "e_bp", "e_ap"],
    }


def test_worked_example_with_parity():
    _, gw = images(GroupSpec((("a", "b"), ("p",)), parity=True), 2)
    assert gw == {
        "e": ["d_aa", "d_bb", "d_pp", "d_ab", "d_ap", "d_bp", "e_ab", "e_ap", "e_bp"],
        "(a b)": ["d_bb", "d_aa", "d_pp", "d_ab", "d_bp", "d_ap", "-e_ab", "e_bp", "e_ap"],
        "p": ["d_aa", "d_bb", "d_pp", "d_ab", "d_ap", "d_bp", "-e_ab", "-e_ap", "-e_bp"],
        "(a b)p": ["d_bb", "d_aa", "d_pp", "d_ab", "d_bp", "d_ap", "e_ab", "-e_bp", "-e_ap"],
    }


def test_three_cycle_sign_of_epsilon():
    # Under the cycle a->b->c->a, e_abc -> e_bca = +e_abc (an even permutation).
    spec = GroupSpec((("a", "b", "c"),))
    _, gw = images(spec, 3)
    e_slot = gw["e"].index("e_abc")
    assert gw["(a b c)"][e_slot] == "e_abc"
    assert gw["(a b)"][e_slot] == "-e_abc"


@pytest.mark.parametrize("species,parity,order", [
    ((("a",),), False, 1),
    ((("a", "b", "c"),), False, 6),
    ((("a", "b", "c"), ("p", "q")), False, 12),
    ((("a", "b", "c"), ("p", "q")), True, 24),
])
def test_group_order_and_identity_first(species, parity, order):
    spec = GroupSpec(species, parity)
    elements = tuple(iter_group_elements(spec))
    assert len(elements) == spec.order == order
    assert elements[0].is_identity()
    assert len(set(elements)) == order


def test_slot_counts():
    from math import comb
    for n, k in [(3, 2), (4, 3), (5, 4), (5, 6)]:
        spec = GroupSpec((tuple("abcdefgh"[:n]),))
        assert len(build_slots(spec, k)) == n + comb(n, 2) + comb(n, k)


def test_group_action_is_a_signed_permutation():
    spec = GroupSpec((("a", "b", "c"), ("p", "q")), parity=True)
    slots = build_slots(spec, 3)
    index, sign = signed_permutation_table(slots, tuple(iter_group_elements(spec)))
    for row_index, row_sign in zip(index, sign):
        assert sorted(row_index) == list(range(len(slots)))
        assert set(row_sign) <= {-1, 1}


def test_bad_specs():
    with pytest.raises(ValueError):
        GroupSpec((("a", "b"), ("b",)))
    with pytest.raises(ValueError):
        GroupSpec(((),))
    with pytest.raises(TypeError):
        GroupSpec((["a", "b"],))
