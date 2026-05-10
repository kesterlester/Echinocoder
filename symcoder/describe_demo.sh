#!/bin/sh

cd /Users/lester/github/Echinocoder && venv/bin/python -c "
import numpy as np
from symatom import ArgumentSymmetry, VectorGroup, Context, Plan
from symcoder import EvaluableOperation, encode, describe_encoding

dot  = EvaluableOperation('dot',  rank=2, parity=+1, argument_symmetry=ArgumentSymmetry.SYMMETRIC,
                           eval_fn=lambda v: float(np.dot(v[0], v[1])))
eps3 = EvaluableOperation('eps3', rank=3, parity=-1, argument_symmetry=ArgumentSymmetry.ANTISYMMETRIC,
                           eval_fn=lambda v: float(np.dot(v[0], np.cross(v[1], v[2]))))

electrons = VectorGroup('electrons', ('a','b','c'))
muons     = VectorGroup('muons',     ('p','q'))
ctx  = Context((electrons, muons, ))
plan = Plan(context=ctx, operations=(dot, eps3))

segs = describe_encoding(plan)
for s in segs:
    print(s)

print()
event = {l: np.random.randn(3) for l in ('a','b','c','p','q',)}
out = encode(plan, event)
print(f'Total output length: {len(out)}  (sum of segment lengths: {sum(s.length for s in segs)})')
print()
import json
print(json.dumps([s.to_dict() for s in segs[:2]], indent=2))
"
