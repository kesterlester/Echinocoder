from __future__ import annotations
import importlib.util
import numpy as np
from pathlib import Path
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


def encode(plan, event: dict) -> np.ndarray:
    """
    Encode a physics event as a permutation-invariant vector.

    Algorithm:
      1. Enumerate all canonical PairFlavours for this plan.
      2. For each PairFlavour, evaluate the full G-orbit of its representative
         atom-pair: z_k = eval(u'_k, E) + i * eval(v'_k, E).
      3. Polynomial-embed {z_k} via the Echinocoder zip embedder to get a
         fixed-length coefficient vector for this PairFlavour.
      4. Concatenate all coefficient vectors.

    Returns a 1-D numpy array of complex dtype.  Each PairFlavour contributes
    exactly orbit_size(group_sizes) coefficients (the polynomial has degree
    orbit_size, and we drop the leading 1).

    If the plan has no registered operations, returns an empty array.
    """
    fo_list = repL(plan.context, plan.operations)
    pf_list = canonical_pair_flavours(fo_list, plan.context)
    if not pf_list:
        return np.array([], dtype=complex)

    parts = []
    for pf in pf_list:
        z = np.array(eval_pair_orbit(pf, plan, event), dtype=complex)
        coeffs, _, _ = _zip_embed(z)
        parts.append(coeffs)

    return np.concatenate(parts)
