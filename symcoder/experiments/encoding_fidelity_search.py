"""
symcoder.experiments.encoding_fidelity_search
=============================================

Empirical search for "shadow events" — pairs (E1, E2) of physical events that:

  (a) produce *identical* Echinocoder Phase 1 + Phase 2 encodings, AND
  (b) are NOT equivalent under the action of (S_n × SO(3)).

If any such pair exists, it proves that the current pair-level multiset
encoding loses information beyond the desired SO(3) × label-permutation
invariances — i.e. the encoding itself is information-incomplete, and
faithful encoding would require Phase 3 (triple-level) data.

Equivalence oracle (independent of Echinocoder, derived from classical
invariant theory):

  E1 ~ E2  iff  there exists σ ∈ S_n (respecting particle types) such that
                  σ(Gram(E1))     == Gram(E2)        (SO(3)-invariant)
                  σ(SignedDet(E1)) == SignedDet(E2)  (SO(3)-invariant, det=+1)

Run with:
    source venv/bin/activate
    python -m symcoder.experiments.encoding_fidelity_search
"""
from __future__ import annotations

import itertools as it
import time
from collections import defaultdict
from typing import Dict, Tuple

import numpy as np

from symatom import ArgumentSymmetry, Operation, VectorType, Context, Plan
import symcoder.operations.euclidean2
import symcoder.operations.euclidean3
from symcoder.encoders import (
    OrbitEncoderFactory, SortEncoderFactory, HalfSortEncoderFactory,
    standard_row_pair_factories, OverlapBlockEncoderFactory, Phase2EncoderFactory,
)


# ---------------------------------------------------------------------------
# Plan setup (same as the adversarial test)
# ---------------------------------------------------------------------------

ELECTRON_LABELS = ("a", "b", "c")
MUON_LABELS     = ("p", "q")
ALL_LABELS      = ELECTRON_LABELS + MUON_LABELS
N_PARTICLES     = len(ALL_LABELS)


def make_plan():
    """Pure 3D Euclidean plan: mag + 3D dot + 3D eps.

    Deliberately omits ``euclidean2.dot`` so the plan's intended invariance
    is exactly SO(3) × (S_electrons × S_muons), which matches the
    equivalence oracle (``canonical_invariant``).  Adding euclidean2.dot
    would break SO(3) symmetry (the operation only sees the first 2 coords)
    and require a different, narrower oracle.
    """
    mymag = Operation(
        "Mymag", rank=1, odd_parity=False,
        argument_symmetry=ArgumentSymmetry.SYMMETRIC,
        eval_fn=lambda v: float(np.sqrt(np.dot(v[0], v[0]))),
        mass_dimension=1,
    )
    electrons = VectorType("electrons", ELECTRON_LABELS)
    muons     = VectorType("muons",     MUON_LABELS)
    ctx       = Context((electrons, muons))
    plan      = Plan(context=ctx, operations=(
        mymag,
        symcoder.operations.euclidean3.dot,
        symcoder.operations.euclidean3.eps,
    ))
    return plan, ctx


# ---------------------------------------------------------------------------
# Encoders (built once, reused)
# ---------------------------------------------------------------------------

def build_encoders(plan):
    orbit_fac  = OrbitEncoderFactory([HalfSortEncoderFactory(), SortEncoderFactory()])
    phase2_fac = Phase2EncoderFactory([
        OverlapBlockEncoderFactory(standard_row_pair_factories(),
                                   use_complementarity_drop=False)
    ])
    return orbit_fac.build(plan), phase2_fac.build(plan)


def encode(orbit_enc, phase2_enc, event):
    """Return the full encoding as a single np.ndarray (Phase 1 then Phase 2)."""
    ph1 = orbit_enc.encode(event).values
    ph2 = phase2_enc.encode(event).values
    return np.concatenate([np.asarray(ph1, dtype=float),
                           np.asarray(ph2, dtype=float)])


# ---------------------------------------------------------------------------
# Independent equivalence oracle (Gram + signed determinants up to S_n)
# ---------------------------------------------------------------------------

def gram_and_dets(event) -> Tuple[np.ndarray, np.ndarray]:
    """Return (G, D) where G[i,j] = v_i·v_j and D[i,j,k] = det([v_i,v_j,v_k]).

    Particles are ordered as ELECTRON_LABELS + MUON_LABELS.  Vectors of length
    other than 3 are zero-padded so the determinant is well defined; for our
    setup every vector is already in R^3.
    """
    V = np.array([np.asarray(event[lbl], dtype=float) for lbl in ALL_LABELS])
    # Pad to 3 dimensions if needed (not expected for this plan).
    if V.shape[1] < 3:
        pad = np.zeros((V.shape[0], 3 - V.shape[1]))
        V = np.concatenate([V, pad], axis=1)
    G = V @ V.T
    n = V.shape[0]
    D = np.zeros((n, n, n))
    for i in range(n):
        for j in range(n):
            for k in range(n):
                D[i, j, k] = float(np.linalg.det(np.stack([V[i], V[j], V[k]])))
    return G, D


def typed_permutations() -> list:
    """All σ ∈ S_electrons × S_muons, expressed as index permutations of ALL_LABELS."""
    n_e = len(ELECTRON_LABELS)
    n_m = len(MUON_LABELS)
    perms = []
    for pe in it.permutations(range(n_e)):
        for pm in it.permutations(range(n_e, n_e + n_m)):
            perm = list(pe) + list(pm)
            perms.append(tuple(perm))
    return perms


TYPED_PERMS = typed_permutations()  # 3! × 2! = 12 of them


def _normalise_zeros(arr: np.ndarray) -> np.ndarray:
    """Replace -0.0 with +0.0 (otherwise byte-comparison spuriously differs)."""
    out = arr.copy()
    out[out == 0.0] = 0.0
    return out


def canonical_invariant(event, decimals: int = 6, group: str = "SO3") -> bytes:
    """Lex-min over typed permutations of (rounded G, rounded D), as bytes.

    Two events have equal canonical_invariant iff they are S_n × G equivalent
    (modulo the chosen rounding tolerance), where G ∈ {"SO3", "O3"}.

    For G = "O3", the orientation can flip globally: each candidate considers
    both D and -D and picks the lex-smaller, so any σ that matches G and gives
    ±D-matched events collapses to the same canonical key.
    """
    G_mat, D = gram_and_dets(event)
    Gr = _normalise_zeros(np.round(G_mat, decimals))
    Dr = _normalise_zeros(np.round(D, decimals))
    Dr_neg = _normalise_zeros(-Dr)
    best: bytes | None = None
    for perm in TYPED_PERMS:
        p = np.array(perm)
        Gp = Gr[np.ix_(p, p)]
        Dp = Dr[np.ix_(p, p, p)]
        candidate = Gp.tobytes() + b"|" + Dp.tobytes()
        if best is None or candidate < best:
            best = candidate
        if group == "O3":
            Dpn = Dr_neg[np.ix_(p, p, p)]
            cand2 = Gp.tobytes() + b"|" + Dpn.tobytes()
            if cand2 < best:
                best = cand2
    assert best is not None
    return best


def encoding_key(vec: np.ndarray, decimals: int = 6) -> bytes:
    """Hashable, tolerance-rounded encoding key.

    Phase 1 and Phase 2 outputs are already canonicalised by the encoders
    (e.g. sorting), so identical events produce identical arrays modulo
    floating-point noise.
    """
    return _normalise_zeros(np.round(vec, decimals)).tobytes()


# ---------------------------------------------------------------------------
# Event generators
# ---------------------------------------------------------------------------

def _to_event(coords: np.ndarray) -> Dict[str, np.ndarray]:
    """Build an event dict from a (N_PARTICLES, 3) coord array."""
    return {lbl: coords[i].copy() for i, lbl in enumerate(ALL_LABELS)}


def random_events(n: int, rng: np.random.Generator, low: float = -1.0, high: float = 1.0):
    for _ in range(n):
        coords = rng.uniform(low, high, size=(N_PARTICLES, 3))
        yield _to_event(coords)


def random_antipodal_muons(n: int, rng: np.random.Generator):
    """Adversarial muon pattern: q = -p, with p on unit sphere; electrons random."""
    for _ in range(n):
        # Random unit vector for p
        v = rng.standard_normal(3)
        v /= np.linalg.norm(v)
        coords = np.empty((N_PARTICLES, 3))
        coords[:3] = rng.uniform(-1.0, 1.0, size=(3, 3))  # electrons
        coords[3]  = v   # p
        coords[4]  = -v  # q
        yield _to_event(coords)


def discrete_events_small(values=(-1, 0, 1)):
    """All events with coords drawn from a small discrete set.

    For 5 particles × 3 dims × |values|=3, this is 3^15 ≈ 14M events.
    Too many to enumerate; we sample structured subsets instead.
    """
    raise NotImplementedError("Use sampled discrete generators below.")


def random_discrete_events(n: int, rng: np.random.Generator,
                           values=(-1, 0, 1)):
    """Sample events with coords from a small discrete set."""
    vals = np.array(values, dtype=float)
    for _ in range(n):
        idx = rng.integers(0, len(vals), size=(N_PARTICLES, 3))
        coords = vals[idx]
        yield _to_event(coords)


def random_discrete_antipodal(n: int, rng: np.random.Generator,
                              values=(-1, 0, 1)):
    """Discrete coords for electrons; muons antipodal unit vectors along
    one of the 3 coordinate axes."""
    vals = np.array(values, dtype=float)
    axes = np.eye(3)
    for _ in range(n):
        idx = rng.integers(0, len(vals), size=(3, 3))   # electrons
        coords = np.empty((N_PARTICLES, 3))
        coords[:3] = vals[idx]
        ax = axes[rng.integers(0, 3)]
        coords[3] = ax
        coords[4] = -ax
        yield _to_event(coords)


def random_collinear_electrons(n: int, rng: np.random.Generator):
    """Electrons collinear (b, c are multiples of a)."""
    for _ in range(n):
        a = rng.uniform(-1, 1, size=3)
        scales = rng.uniform(-2, 2, size=2)
        coords = np.empty((N_PARTICLES, 3))
        coords[0] = a
        coords[1] = scales[0] * a
        coords[2] = scales[1] * a
        # Antipodal muons on z
        coords[3] = np.array([0.0, 0.0,  1.0])
        coords[4] = np.array([0.0, 0.0, -1.0])
        yield _to_event(coords)


def coplanar_electrons(n: int, rng: np.random.Generator):
    """Electrons all have z=0 (lie in xy-plane)."""
    for _ in range(n):
        coords = np.empty((N_PARTICLES, 3))
        for i in range(3):
            coords[i] = np.array([rng.uniform(-1, 1), rng.uniform(-1, 1), 0.0])
        coords[3] = np.array([0.0, 0.0,  1.0])
        coords[4] = np.array([0.0, 0.0, -1.0])
        yield _to_event(coords)


# ---------------------------------------------------------------------------
# Search driver
# ---------------------------------------------------------------------------

def run_search(generators, decimals: int = 6, verbose: bool = True,
               oracle_group: str = "SO3"):
    """Drive a search: encode every event from each generator, bucket by
    encoding key, then within each bucket check whether all events share
    the same canonical_invariant.

    A "shadow pair" is a bucket containing >=2 events with distinct
    canonical invariants.
    """
    plan, _ = make_plan()
    orbit_enc, phase2_enc = build_encoders(plan)

    buckets: Dict[bytes, list] = defaultdict(list)
    total = 0
    t0 = time.time()

    for label, gen in generators:
        n_local = 0
        for event in gen:
            try:
                vec = encode(orbit_enc, phase2_enc, event)
            except Exception as e:
                if verbose:
                    print(f"[{label}] encode failed: {e}")
                continue
            key = encoding_key(vec, decimals=decimals)
            inv = canonical_invariant(event, decimals=decimals, group=oracle_group)
            buckets[key].append((label, event, inv))
            n_local += 1
            total  += 1
        if verbose:
            print(f"[{label}] encoded {n_local} events"
                  f" (cumulative {total}, {time.time()-t0:.1f}s)")

    if verbose:
        print(f"\nTotal events: {total}")
        print(f"Distinct encoding buckets: {len(buckets)}")
        bucket_sizes = sorted((len(v) for v in buckets.values()), reverse=True)
        print(f"Top bucket sizes: {bucket_sizes[:10]}")

    shadow_pairs = []
    for key, entries in buckets.items():
        if len(entries) < 2:
            continue
        # Collect distinct canonical invariants in this bucket
        by_inv: Dict[bytes, list] = defaultdict(list)
        for label, event, inv in entries:
            by_inv[inv].append((label, event))
        if len(by_inv) >= 2:
            shadow_pairs.append((key, by_inv))

    return shadow_pairs, buckets


def report(shadow_pairs, buckets, verbose: bool = True):
    print("\n" + "=" * 70)
    if not shadow_pairs:
        print("NO SHADOW PAIRS FOUND.")
        print(f"  ({len(buckets)} distinct encoding buckets searched.)")
        print(f"  This is *evidence* that pair-level encoding is faithful")
        print(f"  for events in the searched families, but NOT a proof.")
        return

    print(f"FOUND {len(shadow_pairs)} SHADOW BUCKET(S)")
    print("Each bucket has events that encode identically but are NOT")
    print("(S_n × SO(3))-equivalent.  This would be a real counterexample")
    print("to pair-level encoding faithfulness.")
    print("=" * 70)
    for i, (key, by_inv) in enumerate(shadow_pairs[:5]):
        print(f"\n--- Shadow bucket #{i+1} ---")
        print(f"  encoding key hash: {hash(key) & 0xffffffff:08x}")
        print(f"  distinct (S_n × SO(3)) classes in this bucket: {len(by_inv)}")
        for j, (inv, evs) in enumerate(by_inv.items()):
            label, ev = evs[0]
            print(f"  class #{j+1} (from family '{label}'):")
            for lbl in ALL_LABELS:
                print(f"    {lbl} = {np.round(ev[lbl], 4).tolist()}")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    rng = np.random.default_rng(20260525)

    # Generators for the main search.
    rng2 = np.random.default_rng(99)
    bigger_gens = lambda: [
        ("random_uniform",           random_events(3000, rng2)),
        ("random_antipodal_muons",   random_antipodal_muons(3000, rng2)),
        ("random_discrete_011",      random_discrete_events(3000, rng2,
                                                            values=(-1, 0, 1))),
        ("random_discrete_antipodal",random_discrete_antipodal(3000, rng2,
                                                               values=(-1, 0, 1))),
        ("collinear_electrons",      random_collinear_electrons(500, rng2)),
        ("coplanar_electrons",       coplanar_electrons(500, rng2)),
    ]

    print("=== STAGE 1: SO(3) × S_n oracle ===")
    rng2 = np.random.default_rng(99)
    sp_so3, buckets_so3 = run_search(bigger_gens(), decimals=6,
                                     verbose=True, oracle_group="SO3")
    report(sp_so3, buckets_so3)

    print("\n=== STAGE 2: O(3) × S_n oracle (reflections allowed) ===")
    rng2 = np.random.default_rng(99)
    sp_o3, buckets_o3 = run_search(bigger_gens(), decimals=6,
                                    verbose=True, oracle_group="O3")
    report(sp_o3, buckets_o3)

    print("\n" + "=" * 70)
    print("INTERPRETATION")
    print("=" * 70)
    print(f"Shadow buckets under SO(3) × S_n oracle : {len(sp_so3)}")
    print(f"Shadow buckets under O(3)  × S_n oracle : {len(sp_o3)}")
    print()
    if len(sp_so3) > 0 and len(sp_o3) == 0:
        print("  All shadows under SO(3) collapse under O(3): the encoder")
        print("  is silently O(3)-invariant (loses chirality / handedness).")
        print("  This is a real fidelity claim that contradicts the inclusion")
        print("  of eps3 — worth flagging as an encoder property, not a bug")
        print("  in pair-level multiset representation.")
    elif len(sp_o3) > 0:
        print("  Some shadows survive even O(3)-equivalence: those are")
        print("  GENUINE pair-fidelity counterexamples — events that are")
        print("  not related by any orthogonal transformation + label perm")
        print("  yet encode identically.  These would prove pair-level")
        print("  encoding is informationally incomplete.")
    else:
        print("  No shadows under either oracle: encoding is faithful at")
        print("  pair level for all event families searched.")
