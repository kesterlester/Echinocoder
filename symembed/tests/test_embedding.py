"""Layer 2 and the evaluator: values, invariance, shapes, gradients, and metadata."""
import json

import numpy as np
import pytest
import sympy
import torch

from symembed import EmbeddingOptions, GroupSpec, build_invariant_embedding

F64 = EmbeddingOptions(dtype=torch.float64)

EXAMPLE = GroupSpec((("a", "b"), ("p",)))
EXAMPLE_P = GroupSpec((("a", "b"), ("p",)), parity=True)
TWO_PAIRS = GroupSpec((("a", "b"), ("p", "q")))

# (k, spec) cases small enough to expand symbolically.
SMALL_CASES = [(2, EXAMPLE), (2, EXAMPLE_P), (3, TWO_PAIRS), (1, GroupSpec((("a", "b", "c"),))),
               (3, GroupSpec((("a", "b"), ("p",)), parity=True)), (4, GroupSpec((("a", "b"), ("p",))))]


def random_rotation(k, generator, det=1.0):
    q, r = torch.linalg.qr(torch.randn(k, k, dtype=torch.float64, generator=generator))
    q = q * torch.sign(torch.diagonal(r))
    if torch.linalg.det(q) * det < 0:
        q[:, 0] = -q[:, 0]
    return q


@pytest.fixture
def gen():
    return torch.Generator().manual_seed(1234)


# ----- values ----------------------------------------------------------------------------------------------------

def test_worked_example_values(gen):
    embed, meta = build_invariant_embedding(2, EXAMPLE, F64)
    assert meta.D == 26
    X = torch.randn(5, 3, 2, dtype=torch.float64, generator=gen)
    Y = embed(X)
    w = dict(zip(meta.slot_names(), embed.w(X).unbind(-1)))
    expected_E1 = [w["d_aa"] + w["d_bb"], w["d_aa"] + w["d_bb"], 2 * w["d_pp"], 2 * w["d_ab"],
                   w["d_ap"] + w["d_bp"], w["d_ap"] + w["d_bp"], 0 * w["d_aa"],
                   w["e_ap"] + w["e_bp"], w["e_ap"] + w["e_bp"]]
    torch.testing.assert_close(Y[:, :9], torch.stack(expected_E1, -1))
    E_2_12 = 2 * (w["d_ap"] * w["e_ap"] + w["d_bp"] * w["e_bp"]) - w["e_ab"] ** 2
    torch.testing.assert_close(Y[:, 9 + 12], E_2_12)


def test_w_matches_definitions(gen):
    embed, meta = build_invariant_embedding(3, TWO_PAIRS, F64)
    X = torch.randn(3, 4, 3, dtype=torch.float64, generator=gen)
    w = embed.w(X)
    pos = {label: i for i, label in enumerate(TWO_PAIRS.labels)}
    for j, slot in enumerate(meta.slots):
        vecs = X[:, list(slot.positions), :]
        if slot.kind == "d":
            expected = (vecs[:, 0] * vecs[:, 1]).sum(-1)
        else:
            expected = torch.linalg.det(vecs)
        torch.testing.assert_close(w[:, j], expected)
    assert pos  # labels are positional


@pytest.mark.parametrize("k,spec", SMALL_CASES)
def test_numeric_matches_symbolic_expansion(k, spec, gen):
    embed, meta = build_invariant_embedding(k, spec, F64)
    X = torch.randn(spec.n, k, dtype=torch.float64, generator=gen)
    Y = embed(X)
    w_values = dict(zip(meta.slot_symbols(), embed.w(X).tolist()))
    expanded = meta.expanded_outputs()
    assert len(expanded) == meta.D == Y.shape[-1]
    expected = torch.tensor([float(sympy.sympify(e).subs(w_values)) for e in expanded], dtype=torch.float64)
    torch.testing.assert_close(Y, expected, rtol=1e-10, atol=1e-10)


def test_matches_legacy_bursar_embedder(gen):
    """
    Cross-check layer 2 against Echinocoder's original sympy Bursarial embedder, fed the multiset {g.w}. The legacy
    embedder builds sympy Poly(vector, x), which puts component 0 on the HIGHEST power of x, whereas q_g puts slot 0 on
    x^0; so it is fed each g.w reversed.
    """
    legacy = pytest.importorskip("Cinf_sympy_bursar_embedder_for_array_of_reals_as_multiset")
    embed, meta = build_invariant_embedding(3, TWO_PAIRS, F64)
    X = torch.randn(4, 3, dtype=torch.float64, generator=gen)
    w = embed.w(X)
    gw = (w[embed.g_index] * embed.g_sign).flip(-1).tolist()
    expected, _ = legacy.Embedder().embed_generic(np.array(gw))
    torch.testing.assert_close(embed(X), torch.tensor([float(v) for v in expected], dtype=torch.float64))


# ----- invariance and separation ---------------------------------------------------------------------------------

@pytest.mark.parametrize("parity", [False, True])
def test_invariance(parity, gen):
    spec = GroupSpec((("a", "b", "c"), ("p", "q")), parity=parity)
    k = 3
    embed, meta = build_invariant_embedding(k, spec, F64)
    X = torch.randn(spec.n, k, dtype=torch.float64, generator=gen)
    Y = embed(X)
    scale = Y.abs().max()

    # Permutations within species leave the output unchanged.
    for perm in ([1, 2, 0, 3, 4], [0, 1, 2, 4, 3], [2, 0, 1, 4, 3]):
        torch.testing.assert_close(embed(X[perm]), Y, rtol=1e-9, atol=1e-9 * scale)
    # Rotations leave it unchanged.
    R = random_rotation(k, gen, det=1.0)
    torch.testing.assert_close(embed(X @ R.T), Y, rtol=1e-9, atol=1e-9 * scale)
    # Reflections leave it unchanged only if parity is in G.
    S = random_rotation(k, gen, det=-1.0)
    Y_reflected = embed(X @ S.T)
    if parity:
        torch.testing.assert_close(Y_reflected, Y, rtol=1e-9, atol=1e-9 * scale)
    else:
        assert (Y_reflected - Y).abs().max() > 1e-6 * scale
    # Swapping vectors between species changes it.
    assert (embed(X[[3, 1, 2, 0, 4]]) - Y).abs().max() > 1e-6 * scale


def test_distinct_orbits_give_distinct_outputs(gen):
    # A smoke test only (random inputs), not evidence of injectivity on every input.
    embed, _ = build_invariant_embedding(2, EXAMPLE, F64)
    X = torch.randn(50, 3, 2, dtype=torch.float64, generator=gen)
    Y = embed(X)
    distances = torch.cdist(Y, Y) + torch.eye(50, dtype=torch.float64)
    assert distances.min() > 1e-6


# ----- shapes, dtypes, interface ---------------------------------------------------------------------------------

def test_shapes_and_input_types(gen):
    embed, meta = build_invariant_embedding(2, EXAMPLE, F64)
    assert embed.input_shape.tolist() == [3, 2] and embed.output_shape.tolist() == [meta.D]
    X = torch.randn(2, 4, 3, 2, dtype=torch.float64, generator=gen)
    Y = embed(X)
    assert Y.shape == (2, 4, meta.D)
    torch.testing.assert_close(embed(X[0, 0]), Y[0, 0])
    torch.testing.assert_close(embed(X.numpy()), Y)
    with pytest.raises(ValueError):
        embed(torch.randn(4, 2, 3))


def test_flattened_variant(gen):
    embed, meta = build_invariant_embedding(2, EXAMPLE, F64)
    flat = embed.flattened()
    assert flat.input_shape.tolist() == [6] and flat.output_shape.tolist() == [meta.D]
    X = torch.randn(7, 3, 2, dtype=torch.float64, generator=gen)
    torch.testing.assert_close(flat(X.reshape(7, 6)), embed(X))


def test_dtype_option_and_conversion(gen):
    embed32, _ = build_invariant_embedding(2, EXAMPLE)
    assert embed32.g_sign.dtype == torch.float32
    Y32 = embed32(np.random.default_rng(0).normal(size=(3, 3, 2)))
    assert Y32.dtype == torch.float32
    assert embed32.to(torch.float64).g_sign.dtype == torch.float64


def test_gradients(gen):
    embed, _ = build_invariant_embedding(2, EXAMPLE, F64)
    X = torch.randn(2, 3, 2, dtype=torch.float64, generator=gen, requires_grad=True)
    assert torch.autograd.gradcheck(embed, (X,))


def test_unsupported_options():
    for kwargs in ({"embedder": "other"}, {"reductions": "safe"}, {"homogenise": True},
                   {"root_normalise": True}, {"field": "complex"}):
        with pytest.raises(NotImplementedError):
            EmbeddingOptions(**kwargs)


def test_refuses_infeasible_sizes():
    with pytest.raises(ValueError, match="max_outputs"):
        build_invariant_embedding(5, GroupSpec((tuple("abcdefghij"),)))


# ----- metadata --------------------------------------------------------------------------------------------------

def test_meta_serialises(gen):
    _, meta = build_invariant_embedding(2, EXAMPLE_P)
    d = json.loads(meta.to_json())
    assert d["D"] == meta.D == 4 + 8 * 10
    assert d["group_order"] == 4
    assert d["certification"].startswith("Level 0")
    assert [s["name"] for s in d["slots"]][:3] == ["d_aa", "d_bb", "d_pp"]
    assert d["outputs"][0]["name"] == "E_{1,0}"


def test_possible_L_degrees():
    _, meta = build_invariant_embedding(3, GroupSpec((("a", "b", "c"),)))
    d = meta.to_dict()
    # c = 1: each output is linear in one slot orbit, so its degree is that slot's degree.
    for out in d["outputs"][:meta.m]:
        assert out["possible_L_degrees"] == [d["slots"][out["j"]]["L_degree"]]
    # Some c = 2 output mixes d*d (degree 4) with d*e (degree 5).
    assert any(out["possible_L_degrees"] == [4, 5] for out in d["outputs"] if out["c"] == 2)


def test_latex_mentions_everything():
    _, meta = build_invariant_embedding(2, EXAMPLE)
    text = meta.latex(expanded=True)
    assert "d_{aa}" in text and "E_{2,12}" in text and "(a b)" in text
    assert "Level 0" in text
