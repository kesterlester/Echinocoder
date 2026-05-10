from __future__ import annotations
import importlib.util
import numpy as np
from pathlib import Path
from symatom.atoms import ArgumentSymmetry
from symatom.rep import canonical_pair_flavours
from symatom import repL
from .pairs import eval_pair_orbit

# Load the Echinocoder zip embedder from the repo root.  It is not a proper
# package, so we locate it relative to this file (repo_root/symcoder/encode.py
# → repo_root/Cinf_numpy_polynomial_embedder_for_list_of_reals_as_multiset.py).
_EMBEDDER_PATH = (
    Path(__file__).resolve().parent.parent
    / "Cinf_numpy_polynomial_embedder_for_list_of_reals_as_multiset.py"
)
_spec = importlib.util.spec_from_file_location("_echinocoder", _EMBEDDER_PATH)
_echinocoder = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_echinocoder)
_zip_embed = _echinocoder.embed


def _force_positive_side_key(op, flavour):
    """
    Canonical key for one side of a pair, collapsing ANTISYMMETRIC/SYMMETRIC
    distinction into the same key (forcePositive: negating all atoms in an
    ANTISYMMETRIC row gives an equivalent z-vector, so both the signed and
    unsigned variants of that row are treated as the same 'side').
    """
    return (op.name, tuple(flavour.counts))


def _encoding_canonical_key(pf):
    """
    Canonical key for a PairFlavour under encoding-level equivalence.

    Two pairs of rows are encoding-equivalent if they produce z-vectors whose
    polynomial embeddings carry the same information, up to:
      - u↔v swap        (z_k → i·conj(z_k)):   handled by sorting the sides
      - negating row u  (z_k → −Re + i·Im):     handled by forcePositive
      - negating row v  (z_k → conj(z_k)):       handled by forcePositive

    The key is (sorted pair of side-keys, overlap), which is invariant under
    all three operations.
    """
    u_key = _force_positive_side_key(pf.op_u, pf.flavour_u)
    v_key = _force_positive_side_key(pf.op_v, pf.flavour_v)
    return (tuple(sorted([u_key, v_key])), tuple(pf.overlap))


def encode(plan, event: dict) -> np.ndarray:
    """
    Encode a physics event as a permutation-invariant vector.

    Algorithm:
      1. Enumerate all canonical PairFlavours for this plan.
      2. Deduplicate under encoding equivalence: pairs related by u↔v swap,
         negating row u, or negating row v carry the same information — keep
         only the first representative from each equivalence class (seen_keys).
      3. For each kept PairFlavour, evaluate the G-orbit of its representative
         atom-pair: z_k = eval(u'_k, E) + i * eval(v'_k, E).
      4. Polynomial-embed {z_k} and concatenate coefficient vectors.

    Returns a 1-D numpy array of complex dtype.
    If the plan has no registered operations, returns an empty array.
    """
    fo_list = repL(plan.context, plan.operations)
    pf_list = canonical_pair_flavours(fo_list, plan.context)
    if not pf_list:
        return np.array([], dtype=complex)

    # Encoding-level deduplication: skip PairFlavours whose z-vectors carry
    # the same information as one already embedded (related by u↔v swap or
    # negation of a row for ANTISYMMETRIC operations).
    #
    # NOTE: It is believed — but not yet formally proved — that this check
    # never fires given the current symatom machinery.  The argument is:
    # _encoding_canonical_key is injective on PairFlavours (same-named
    # operations are equal, same-counted Flavours are equal), and
    # canonical_pair_flavours already returns a duplicate-free list — so
    # seen_keys would never contain a collision.  The check is retained
    # because (a) the proof has not been written down rigorously, and
    # (b) if the upstream generation ever changes, this safety net matters.
    #
    # Longer exposition in docstring at end of method

    seen_keys = set()
    parts = []
    for pf in pf_list:
        key = _encoding_canonical_key(pf)
        if key in seen_keys:
            print("BEWARE: a deduplication event fired in encode.py despite it being thought that this could/should never happen with the implementation at time of coding.")
            continue
        seen_keys.add(key)
        z = np.array(eval_pair_orbit(pf, plan, event), dtype=complex)
        coeffs, _, _ = _zip_embed(z)
        parts.append(coeffs)

    return np.concatenate(parts)

    """
    Why 8-fold deduplication 

        (+u,+v),(-u,+v),(-u,-v),(+u,-v),(+v,+u),(-v,+u),(-v,-u),(+v,-u)

    might not be needed.

    _encoding_canonical_key(pf) produces (tuple(sorted([(op_u.name, fl_u.counts), (op_v.name, fl_v.counts)])), overlap).
    
    Two Operation objects with the same .name are equal by the Operation dataclass equality (name, rank, parity, argument_symmetry all match). So different operations → different names → different side-keys.
    Two Flavour objects with the same .counts are equal. So different flavours → different count-tuples → different side-keys.
    Therefore: if two PairFlavours differ in any of (op_u, fl_u, op_v, fl_v, overlap), their canonical keys also differ.
    Since canonical_pair_flavours already returns a duplicate-free list, seen_keys is injective over that list — it never fires.
    The argument rests on the assumption that no two registered operations share a name with different behaviour (which would be a user error). Under normal use, the deduplication is provably redundant given the current symatom machinery. Adding the comment you suggested is exactly the right thing:
    
    
    Saying all the above again but in different words:

    _encoding_canonical_key(pf) = (sorted([(op_u.name, fl_u.counts), (op_v.name, fl_v.counts)]), overlap)
    canonical_pair_flavours returns PairFlavours distinguished by (op_u, fl_u, op_v, fl_v, overlap) after canonical ordering — no two entries are equal.
    Two equal Operations always have the same .name (and vice versa, under normal use). Two equal Flavours always have the same .counts.
    Therefore different PairFlavours → different canonical keys → seen_keys is collision-free over the list.
    The one hole: if someone registers two different EvaluableOperation objects with the same name and different eval_fn, those would be considered equal by symatom (which is intentional — eval_fn is excluded from equality) but would also produce the same side-key. In that case seen_keys would suppress a PairFlavour incorrectly. But that's a pre-existing design constraint on how operations are registered, not a new risk introduced here.
    """
