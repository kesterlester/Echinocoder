"""
symatom demo — run with:  python -m symatom.demo
Prints atoms, flavoured operators, and pair flavours for a small physics context.
"""
from .atoms   import ArgumentSymmetry, Operation, VectorType
from .context import Context, Plan
from .canon   import SimpleCanonicaliser
from .rep     import repS, canonical_pair_flavours


# ---------------------------------------------------------------------------
# A typical two-family context
# ---------------------------------------------------------------------------

mass = Operation("m",  rank=1, parity=+1, argument_symmetry=ArgumentSymmetry.SYMMETRIC)
dot  = Operation("dot",  rank=2, parity=+1, argument_symmetry=ArgumentSymmetry.SYMMETRIC)
eps3 = Operation("eps", rank=3, parity=-1, argument_symmetry=ArgumentSymmetry.ANTISYMMETRIC)

electrons = VectorType("electrons", labels=("A", "B", "C", ))
muons     = VectorType("muons",     labels=("P", "Q", ))

ctx  = Context(types=(electrons, muons))
plan = Plan(context=ctx, canonicaliser=SimpleCanonicaliser(), operations=(mass, dot, eps3))


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
        rep = fo.canonical_representative(plan.canonicaliser)
        print(f"  {fo!r}  →  {fo.count()} atoms, in the orbit of {rep!r}")

    _section(f"repS — all atoms")
    all_atoms_and_fos = [(atom, fo) for fo in fo_list for atom in fo.atoms()]
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


if __name__ == "__main__":
    demo()
