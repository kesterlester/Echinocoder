#!/bin/sh

cd /Users/lester/github/Echinocoder && venv/bin/python -c "
import numpy as np
from symatom import ArgumentSymmetry, VectorType, Context, Plan
from symcoder import EvaluableOperation, encode, describe_encoding

mag  = EvaluableOperation('mag',  rank=1, parity=+1, argument_symmetry=ArgumentSymmetry.SYMMETRIC,
                           eval_fn=lambda v: float(np.sqrt(np.dot(v[0], v[0]))))
dot  = EvaluableOperation('dot',  rank=2, parity=+1, argument_symmetry=ArgumentSymmetry.SYMMETRIC,
                           eval_fn=lambda v: float(np.dot(v[0], v[1])))
eps3 = EvaluableOperation('eps3', rank=3, parity=-1, argument_symmetry=ArgumentSymmetry.ANTISYMMETRIC,
                           eval_fn=lambda v: float(np.dot(v[0], np.cross(v[1], v[2]))))

ctx = Context((
   VectorType('electrons', ('a','b','c')),
   VectorType('muons',     ('p','q')),
   #VectorType('jets',      ('u','v','w','x')),
))
plan = Plan(context=ctx, operations=(mag, dot, eps3))

segs = describe_encoding(plan)
print('=== START ===')
for s in segs:
    print(s)
print('=== STOP ===')

print()
event = {l: np.random.randn(3) for l in ctx.all_labels}

# --- stub: pass a registry into encode() ---
# To exercise the encoder registry path, create a registry, register at least
# one encoder, and pass it as the third argument to encode().  Currently
# commented out because the concrete encoders (SortEncoder, PolyEncoder) are
# not yet implemented (they raise NotImplementedError).
#
from symcoder.encoders import AtomOrbitEncoderRegistry, SortEncoder
registry = AtomOrbitEncoderRegistry()
registry.register(SortEncoder())
out = encode(plan, event, registry)
# --- end stub ---

#out = encode(plan, event)

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
