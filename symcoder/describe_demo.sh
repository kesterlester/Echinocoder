#!/bin/sh

cd /Users/lester/github/Echinocoder && venv/bin/python -c "
import numpy as np
from symatom import ArgumentSymmetry, VectorGroup, Context, Plan
from symcoder import EvaluableOperation, encode, describe_encoding

dot  = EvaluableOperation('dot',  rank=2, parity=+1, argument_symmetry=ArgumentSymmetry.SYMMETRIC,
                           eval_fn=lambda v: float(np.dot(v[0], v[1])))
eps3 = EvaluableOperation('eps3', rank=3, parity=-1, argument_symmetry=ArgumentSymmetry.ANTISYMMETRIC,
                           eval_fn=lambda v: float(np.dot(v[0], np.cross(v[1], v[2]))))

ctx = Context((
   VectorGroup('electrons', ('a','b','c')),
   VectorGroup('muons',     ('p','q')),
   VectorGroup('jets',      ('u','v','w','x')),
))
plan = Plan(context=ctx, operations=(dot, eps3))

segs = describe_encoding(plan)
for s in segs:
    print(s)

print()
event = {l: np.random.randn(3) for l in ctx.all_labels}
out = encode(plan, event)
print(f'Total output length: {len(out)}  (sum of segment lengths: {sum(s.length for s in segs)})')
print()
for s in segs[:100]:
    vals = out[s.start:s.stop]
    vals_str = ', '.join(f'{v:.4f}' for v in vals) if len(vals) else '(empty)'
    print(f'{s}')
    print(f'    [{vals_str}]')
print()
import json
print(json.dumps([s.to_dict() for s in segs[:2]], indent=2))
"
