"""Parity tests: verify that _2a_ and _2b_ produce identical embeddings.

The two implementations of Algorithm 2 must produce the same embedding vector
for the same input when _2b_ is configured with its defaults
(preserve_scale_in_step_1=False, preserve_scale_in_step_2=True), which match
_2a_'s fixed configuration.

Also checks that _2b_ with non-default configurations still runs without error
and produces embeddings of the correct shape (the numerical values will differ
from _2a_, which is expected and correct).
"""

import numpy as np
import pytest

from C0HomDeg1_simplicialComplex_embedder_2a_for_array_of_reals_as_multiset import (
    Embedder as Embedder2a,
)
from C0HomDeg1_simplicialComplex_embedder_2b_for_array_of_reals_as_multiset import (
    Embedder as Embedder2b,
)


# ── Shared test inputs ────────────────────────────────────────────────────────

# The canonical worked example from the documentation (n=4, k=3).
DATA_4x3 = np.asarray([[4, 2, 3],
                        [-3, 5, 1],
                        [8, 9, 2],
                        [2, 7, 2]])

# A slightly different matrix to check a second case.
DATA_4x3_B = np.asarray([[4, 3, 3],
                          [-3, 5, 1],
                          [8, 9, 2],
                          [2, 6, 2]])

# Smaller matrix (n=3, k=2) to stress a different shape.
DATA_3x2 = np.asarray([[1, 5],
                        [3, 2],
                        [7, 4]])

# Matrix with negative and zero entries.
DATA_4x2_NEGZERO = np.asarray([[0, -1],
                                [-3, 2],
                                [5, 0],
                                [1, -2]])

# A permutation of DATA_4x3; embeddings must agree (multiset-invariance is
# already tested elsewhere, but confirming parity holds under permutation is cheap).
DATA_4x3_PERMUTED = DATA_4x3[[2, 0, 3, 1], :]

ALL_INPUTS = [DATA_4x3, DATA_4x3_B, DATA_3x2, DATA_4x2_NEGZERO, DATA_4x3_PERMUTED]
INPUT_IDS  = ["4x3", "4x3_B", "3x2", "4x2_negzero", "4x3_permuted"]


# ── Core parity tests ─────────────────────────────────────────────────────────

@pytest.mark.parametrize("data,label", zip(ALL_INPUTS, INPUT_IDS))
def test_2a_2b_default_parity(data, label):
    """_2b_ with default config must produce the same embedding as _2a_."""
    emb2a = Embedder2a()
    emb2b = Embedder2b()  # default: preserve_scale_step1=False, step2=True

    out_2a, _, _meta2a = emb2a.embed(data)
    out_2b, _, _meta2b = emb2b.embed(data)

    assert out_2a.shape == out_2b.shape, (
        f"[{label}] shape mismatch: 2a={out_2a.shape}, 2b={out_2b.shape}"
    )
    np.testing.assert_allclose(
        out_2b, out_2a, rtol=1e-12, atol=1e-12,
        err_msg=f"[{label}] _2a_ and _2b_ (default config) produced different embeddings"
    )


@pytest.mark.parametrize("data,label", zip(ALL_INPUTS, INPUT_IDS))
def test_2a_2b_default_parity_exact(data, label):
    """_2b_ default output should be bit-for-bit identical to _2a_ output."""
    emb2a = Embedder2a()
    emb2b = Embedder2b()

    out_2a, _, _meta2a = emb2a.embed(data)
    out_2b, _, _meta2b = emb2b.embed(data)

    assert np.array_equal(out_2a, out_2b), (
        f"[{label}] _2a_ and _2b_ outputs differ:\n"
        f"  max abs diff = {np.max(np.abs(out_2a - out_2b))}"
    )


# ── Output-size tests ─────────────────────────────────────────────────────────

@pytest.mark.parametrize("data,label", zip(ALL_INPUTS, INPUT_IDS))
def test_2b_output_size(data, label):
    """_2b_ output length must equal size_from_n_k(n, k) = 2nk+1-k."""
    n, k = data.shape
    emb = Embedder2b()
    out, _, _meta = emb.embed(data)
    expected_size = emb.size_from_n_k(n, k)
    assert len(out) == expected_size, (
        f"[{label}] expected size {expected_size}, got {len(out)}"
    )


# ── Non-default configurations run without error ──────────────────────────────

@pytest.mark.parametrize("ps1,ps2", [
    (False, False),
    (True,  False),
    (True,  True),
    # (False, True) is the default, already covered above
])
def test_2b_non_default_configs_run(ps1, ps2):
    """Non-default _2b_ configurations must complete without error."""
    data = DATA_4x3
    n, k = data.shape
    emb = Embedder2b(preserve_scale_in_step_1=ps1, preserve_scale_in_step_2=ps2)
    out, _, _meta = emb.embed(data)
    assert out.shape == (emb.size_from_n_k(n, k),), (
        f"Wrong shape for ps1={ps1}, ps2={ps2}: {out.shape}"
    )


@pytest.mark.parametrize("ps1,ps2", [
    (False, False),
    (True,  False),
    (True,  True),
])
def test_2b_non_default_differs_from_2a(ps1, ps2):
    """Non-default _2b_ configurations should produce different output from _2a_."""
    data = DATA_4x3
    emb2a = Embedder2a()
    emb2b = Embedder2b(preserve_scale_in_step_1=ps1, preserve_scale_in_step_2=ps2)

    out_2a, _, _meta2a = emb2a.embed(data)
    out_2b, _, _meta2b = emb2b.embed(data)

    # These are expected to differ; we just confirm they're not accidentally equal.
    assert not np.array_equal(out_2a, out_2b), (
        f"ps1={ps1}, ps2={ps2}: unexpectedly equal to _2a_ output"
    )


# ── Wrapper selects 2a ────────────────────────────────────────────────────────

def test_wrapper_matches_2a():
    """The _2_ wrapper (defaulting to 2a) must match _2a_ directly."""
    from C0HomDeg1_simplicialComplex_embedder_2_for_array_of_reals_as_multiset import (
        Embedder as EmbedderWrapper,
    )
    data = DATA_4x3
    out_wrapper, _, _ = EmbedderWrapper().embed(data)
    out_2a,      _, _ = Embedder2a().embed(data)
    np.testing.assert_array_equal(out_wrapper, out_2a)


if __name__ == "__main__":
    # Quick smoke-run without pytest
    for data, label in zip(ALL_INPUTS, INPUT_IDS):
        out_2a, _, _ = Embedder2a().embed(data)
        out_2b, _, _ = Embedder2b().embed(data)
        match = np.array_equal(out_2a, out_2b)
        print(f"{label}: {'MATCH' if match else 'DIFFER'}  "
              f"max_abs_diff={np.max(np.abs(out_2a - out_2b)):.2e}")
