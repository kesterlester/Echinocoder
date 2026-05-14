# Atom Orbit Encoder Registry

## What this is

The `symcoder.encoders` subpackage provides a **pluggable registry** of objects
that can embed a symbolic *atom orbit* into an ordered vector of real numbers.

An **atom orbit** is the set of all `Atom` objects obtained by applying the
symmetry group (defined by the `Context`) to a single representative atom.
For example, if the context has two indistinguishable particles {p1, p2} and
the operation is a dot product, the orbit of `dot(p1,p2)` contains just one
element (since dot is symmetric and p1,p2 are interchangeable); whereas the
orbit of `eps(p1,p2)` under a 2-particle antisymmetric group contains
`{+eps(p1,p2), -eps(p2,p1)}`.

The registry lets you register many different encoders, each specialised for
certain kinds of orbit, without the high-level manager needing to know what
encoders exist in advance.

---

## Core concepts

### OrbitSpec — three ways to describe an orbit

The orbit can be communicated to an encoder in any of three forms, wrapped
in an `OrbitSpec` object:

```python
from symcoder.encoders import OrbitSpec

# Form 1: a single representative Atom
spec = OrbitSpec.from_atom(my_atom)

# Form 2: a FlavouredOperator (structural description, no specific atoms needed)
spec = OrbitSpec.from_flavoured_operator(fo)

# Form 3: the full enumerated orbit
spec = OrbitSpec.from_explicit_orbit([atom1, atom2, atom3])
```

Each encoder declares which forms it can handle by returning `can_encode=True`
or `False` from its `assess()` method.

### The two-step protocol

Every encoder supports two operations:

1. **`assess(spec, plan) → EncodingCapability`** — cheap probe; asks "can you
   handle this orbit, and if so, what would the output look like?"  Never
   performs any evaluation of event data.

2. **`encode(spec, event, plan) → EncodingResult`** — the actual embedding;
   returns a float64 numpy array.

### EncodingCapability

```python
@dataclass
class EncodingCapability:
    can_encode:  bool       # True iff this encoder can handle the spec
    output_dim:  int | None # how many reals the embedding would produce
    method_name: str | None # e.g. "sort", "poly_antisym"
    priority:    float      # encoder's self-assessment of fit quality
    metadata:    dict       # anything extra the encoder wants to say
```

### EncodingResult

```python
@dataclass
class EncodingResult:
    values:   np.ndarray   # float64 array of shape (output_dim,)
    metadata: dict         # diagnostic info (method used, orbit size, etc.)
```

---

## The registry

```python
from symcoder.encoders import AtomOrbitEncoderRegistry, SortEncoder, PolyEncoder

registry = AtomOrbitEncoderRegistry()
registry.register(SortEncoder())   # lower priority fallback
registry.register(PolyEncoder())   # higher priority, handles ANTISYMMETRIC orbits
```

### Querying all capable encoders (primary path)

The high-level manager calls `query_all()` to get every capable encoder and
its self-reported capability, then applies its own selection logic:

```python
spec = OrbitSpec.from_flavoured_operator(fo)
candidates = registry.query_all(spec, plan)
# candidates: [(encoder, capability), ...]  — only capable encoders, in
#             registration order

for enc, cap in candidates:
    print(f"  {cap.method_name}: output_dim={cap.output_dim}, "
          f"priority={cap.priority}, metadata={cap.metadata}")

# Manager picks one — here by largest output_dim as a toy criterion:
enc, cap = max(candidates, key=lambda pair: pair[1].output_dim)
result = enc.encode(spec, event, plan)
# result.values is a float64 numpy array of length cap.output_dim
```

### Convenience path (testing / simple pipelines)

```python
# Automatically picks the highest-priority capable encoder:
result = registry.encode(spec, event, plan)

# Or inspect the winner without encoding:
enc, cap = registry.best_for(spec, plan)
```

---

## Writing a new encoder

Subclass `AtomOrbitEncoder` and implement `assess()` and `encode()`:

```python
from symcoder.encoders import (
    AtomOrbitEncoder, OrbitSpec, OrbitSpecForm,
    EncodingCapability, EncodingResult,
)
from symatom.context import Plan
import numpy as np

class MyEncoder(AtomOrbitEncoder):

    def assess(self, spec: OrbitSpec, plan: Plan) -> EncodingCapability:
        # Return can_encode=False for any form or orbit type you cannot handle.
        if spec.form != OrbitSpecForm.FLAVOURED_OPERATOR:
            return EncodingCapability(can_encode=False, output_dim=None,
                                      method_name=None, priority=0.0)
        fo = spec.payload
        # ... inspect fo.operation, fo.flavour, fo.context ...
        return EncodingCapability(
            can_encode=True,
            output_dim=2 * fo.count_of_atoms_one_per_sign(),
            method_name="my_method",
            priority=1.5,
        )

    def encode(self, spec: OrbitSpec, event: dict, plan: Plan) -> EncodingResult:
        fo = spec.payload
        # ... compute embedding ...
        values = np.zeros(2 * fo.count_of_atoms_one_per_sign(), dtype=np.float64)  # placeholder
        return EncodingResult(values=values, metadata={"method": "my_method"})
```

Then register it:

```python
registry.register(MyEncoder())
```

The high-level manager does not need to be modified.  It will find `MyEncoder`
the next time it calls `query_all()`.

---

## Internal delegation between encoders

An encoder may hold another encoder instance and delegate part of its work.
This is useful when one encoder is a refinement of another:

```python
class RefinedEncoder(AtomOrbitEncoder):
    def __init__(self):
        self._base = SortEncoder()   # holds another encoder internally

    def assess(self, spec, plan):
        base_cap = self._base.assess(spec, plan)
        if not base_cap.can_encode:
            return base_cap   # can't even do the base step
        # further checks ...
        return EncodingCapability(can_encode=True, output_dim=base_cap.output_dim // 2,
                                  method_name="refined", priority=2.0)

    def encode(self, spec, event, plan):
        base_result = self._base.encode(spec, event, plan)
        # ... compress base_result.values further ...
        return EncodingResult(values=compressed, metadata={"method": "refined"})
```

`RefinedEncoder` is not registered via the registry's knowledge of `SortEncoder`;
it simply holds a private instance.

---

## Extending to higher-level objects

The same pattern can be repeated for objects that are not single Atoms — for
example, atom pairs, or composite objects that hold Atoms internally.  The
steps are:

1. Define a new `XyzOrbitSpec` class (or reuse `OrbitSpec` with new form values).
2. Define `XyzOrbitEncoder(ABC)` with the same `assess / encode` contract.
3. Define `XyzOrbitEncoderRegistry` (near-identical copy of `AtomOrbitEncoderRegistry`).
4. Write concrete encoder subclasses and register them.

The `AtomOrbitEncoderRegistry` does not need to be modified or subclassed; the
new registry is an independent object.

---

## Files in this subpackage

| File | Contents |
|---|---|
| `_base.py` | `OrbitSpecForm`, `OrbitSpec`, `EncodingCapability`, `EncodingResult`, `AtomOrbitEncoder` ABC |
| `_registry.py` | `AtomOrbitEncoderRegistry` |
| `sort_encoder.py` | `SortEncoder` — sorted evaluation, all orbit forms, low priority |
| `poly_encoder.py` | `PolyEncoder` — polynomial compression for ANTISYMMETRIC orbits |
| `__init__.py` | Public re-exports |
| `README.md` | This file |
