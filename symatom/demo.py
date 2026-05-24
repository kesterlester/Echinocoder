"""
symatom demo — run with either:
    python -m symatom.demo          (from repo root)
    python symatom/demo.py          (from repo root)
    python demo.py                  (from inside symatom/)
"""
from __future__ import annotations
import os, sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from symatom.atoms   import ArgumentSymmetry, Operation, VectorType
from symatom.context import Context, Plan
from symatom.rep     import repS, canonical_pair_flavours, FlavouredOperator, Flavour


# ---------------------------------------------------------------------------
# A typical two-family context
# ---------------------------------------------------------------------------

mass = Operation("m",   rank=1, odd_parity=False, argument_symmetry=ArgumentSymmetry.SYMMETRIC,    eval_fn=lambda v: 0.0)
dot  = Operation("dot", rank=2, odd_parity=False, argument_symmetry=ArgumentSymmetry.SYMMETRIC,    eval_fn=lambda v: 0.0)
eps3 = Operation("eps", rank=3, odd_parity=True,  argument_symmetry=ArgumentSymmetry.ANTISYMMETRIC, eval_fn=lambda v: 0.0)

electrons = VectorType("electrons", labels=("A", "B", "C", ))
muons     = VectorType("muons",     labels=("P", "Q", ))

ctx  = Context(types=(electrons, muons))
plan = Plan(context=ctx, operations=(mass, dot, eps3))


def _section(title: str) -> None:
    print(f"\n{'='*60}")
    print(f"  {title}")
    print(f"{'='*60}")


def demo():
    _section("Context and Operations")
    print(f"  context : {ctx!r}")
    print(f"  plan    : {plan!r}")
    for op in plan.operations:
        print(f"  op      : {op!r}")

    _section("VectorTypes")
    for g in ctx.types:
        print(f"  {g!r}")

    _section("repS — FlavouredOperators")
    fo_list = repS(ctx, plan.operations)
    for fo in fo_list:
        rep = fo.canonical_representative()
        print(f"  {fo!r}  →  {fo.count_of_atoms_one_per_sign()} atoms-one-per-sign, in the orbit of {rep!r}")

    _section(f"repS — all atoms")
    all_atoms_and_fos = [(atom, fo) for fo in fo_list for atom in fo.atoms_one_per_sign()]
    last_fo = None
    for i, (atom,fo) in enumerate(all_atoms_and_fos):
        if last_fo != fo and i !=0:
            print("---------------------")
        last_fo = fo
        print(f"  [{i:3d}]  {atom!r}")
    print(f"\n  total atoms: {len(all_atoms_and_fos)}")

    _section("Canonical PairFlavours")
    pfs = canonical_pair_flavours(fo_list, ctx)
    type_names = [g.name for g in ctx.types]
    for pf in pfs:
        n = pf.count(tuple(g.size for g in ctx.types))
        fl_u = pf.flavour_u.describe(type_names)
        fl_v = pf.flavour_v.describe(type_names)
        print(
            f"  {pf!r}\n"
            f"       u={pf.op_u.name}({fl_u}), v={pf.op_v.name}({fl_v})"
            f"  →  {n} ordered pairs"
        )
    print(f"\n  distinct pair-flavours: {len(pfs)}")


def demo_scrambled_labels():
    """
    Show that atoms_one_per_sign() can yield atoms with MIXED signs when the
    VectorType's label tuple is not in sorted order.

    Root cause: combinations() preserves the declaration order of g.labels.
    If that order is not sorted, the Atom constructor's label-sorting step
    requires an odd permutation for some choices, flipping their sign to -1.

    Example: labels=("zebra","apple","toast") with a rank-2 ANTISYMMETRIC op.
    combinations(("zebra","apple","toast"), 2) yields:
        ("zebra","apple")  →  Atom sorts to ("apple","zebra"), 1 swap → sign=-1
        ("zebra","toast")  →  Atom sorts to ("toast","zebra"), 1 swap → sign=-1
        ("apple","toast")  →  already sorted                         → sign=+1
    Output: a mix of + and - atoms despite all being passed sign=+1.
    """
    eps2 = Operation("eps2", rank=2, odd_parity=True, argument_symmetry=ArgumentSymmetry.ANTISYMMETRIC, eval_fn=lambda v: 0.0)
    scrambled = VectorType("scrambled", labels=("zebra", "apple", "toast"))
    ctx_s = Context(types=(scrambled,))
    fo = FlavouredOperator(operation=eps2, flavour=Flavour((2,)), context=ctx_s)

    _section("atoms_one_per_sign with scrambled labels — mixed sign output")
    print("  VectorType labels (declaration order): ('zebra', 'apple', 'toast')")
    print("  combinations yield pairs in that order, not in sorted order.")
    print()
    atoms = list(fo.atoms_one_per_sign())
    for atom in atoms:
        print(f"  {atom!r}")
    signs = [a.sign for a in atoms]
    print(f"\n  signs: {signs}  ← mix of +1 and -1 despite all constructed with sign=+1")


if __name__ == "__main__":
    demo()
    demo_scrambled_labels()
