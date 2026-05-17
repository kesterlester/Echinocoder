"""
symcoder.encoders._embedder
============================
Load the external Echinocoder zip embedder once at import time and expose it
as `zip_embed`.

This module exists so that the embedder path resolution and loading happens in
exactly one place, and both encode.py (legacy path) and row_pair_encoders.py
(registry path) can import from here without duplication.
"""
import importlib.util
from pathlib import Path

_EMBEDDER_PATH = (
    Path(__file__).resolve().parent.parent.parent
    / "Cinf_numpy_polynomial_embedder_for_list_of_reals_as_multiset.py"
)
_spec = importlib.util.spec_from_file_location("_echinocoder", _EMBEDDER_PATH)
_echinocoder = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_echinocoder)

zip_embed = _echinocoder.embed
