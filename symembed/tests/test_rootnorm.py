"""Root normalisation (option root_normalise)."""
import pytest
import torch

from symembed import EmbeddingOptions, GroupSpec, build_invariant_embedding
from symembed.rootnorm import RootNormaliseMeta

P_CASE = GroupSpec((("a", "b", "c"), ("p", "q")), parity=True)


@pytest.fixture
def gen():
    return torch.Generator().manual_seed(2468)


def test_all_outputs_normalised_and_linear_in_scale(gen):
    embed, meta = build_invariant_embedding(4, P_CASE, EmbeddingOptions(embedder="bursar_graded", reductions="safe",
                                                                         root_normalise=True, dtype=torch.float64))
    assert isinstance(meta, RootNormaliseMeta) and meta.D == 99
    zero = [i for i, z in enumerate(meta.skip) if z]
    assert len(meta.normalised) == 99 - len(zero)
    X = torch.randn(4, 5, 4, dtype=torch.float64, generator=gen)
    keep = meta.normalised
    torch.testing.assert_close(embed(2.5 * X)[:, keep], 2.5 * embed(X)[:, keep])
    torch.testing.assert_close(embed(-X)[:, keep], embed(X)[:, keep])    # even degree in L: invariant under L -> -L
    # ... and the unnormalised embedding is recovered output by output.
    raw, _ = build_invariant_embedding(4, P_CASE, EmbeddingOptions(embedder="bursar_graded", reductions="safe",
                                                                   dtype=torch.float64))
    r = torch.tensor([d[0] if i in meta.normalised else 1 for i, d in enumerate(meta.degrees)], dtype=torch.float64)
    Y = embed(X)
    torch.testing.assert_close(torch.sign(Y) * Y.abs() ** r, raw(X), rtol=1e-9, atol=1e-9 * raw(X).abs().max())


def test_partial_normalisation_without_parity(gen):
    spec = GroupSpec(P_CASE.species)
    embed, meta = build_invariant_embedding(4, spec, EmbeddingOptions(embedder="bursar_graded", root_normalise=True,
                                                                      dtype=torch.float64))
    assert 0 < len(meta.normalised) < meta.D        # outputs mixing d's (degree 2) and e's (degree 4) stay as they are
    X = torch.randn(3, 5, 4, dtype=torch.float64, generator=gen)
    Y, Y2 = embed(X), embed(2.0 * X)
    torch.testing.assert_close(Y2[:, meta.normalised], 2.0 * Y[:, meta.normalised])


def test_gradient_finite_at_zero_outputs():
    embed, _ = build_invariant_embedding(2, GroupSpec((("a", "b"), ("p",)), parity=True),
                                         EmbeddingOptions(embedder="bursar_graded", root_normalise=True,
                                                          reductions="safe", dtype=torch.float64))
    X = torch.zeros(1, 3, 2, dtype=torch.float64, requires_grad=True)
    embed(X).sum().backward()
    assert torch.isfinite(X.grad).all()


def test_with_parity_split(gen):
    spec = GroupSpec(P_CASE.species)
    embed, meta = build_invariant_embedding(3, spec, EmbeddingOptions(embedder="bursar_graded", root_normalise=True,
                                                                      parity_split=True, reductions="safe",
                                                                      dtype=torch.float64))
    X = torch.randn(3, 5, 3, dtype=torch.float64, generator=gen)
    reflection = torch.tensor([1.0, 1.0, -1.0], dtype=torch.float64)
    Y, Yr = embed(X), embed(X * reflection)
    names = meta.output_names()
    # Outputs that are identically zero as functions of L (left unnormalised; mathematically both even and odd) carry
    # only rounding noise, so they are excluded from the parity check.
    zero_positions, position = set(), 0
    for base, ps in enumerate(meta.parities):
        width = 1 if len(ps) <= 1 else 2
        if meta.inner.skip[base]:
            zero_positions.update(range(position, position + width))
        position += width
    odd = [i for i, nm in enumerate(names) if (nm.startswith("odd part") or nm.endswith("[odd]")) and i not in zero_positions]
    even = [i for i in range(len(names)) if i not in odd and i not in zero_positions]
    scale = Y[:, even + odd].abs().max()
    assert (Yr[:, even] - Y[:, even]).abs().max() <= 1e-9 * scale
    assert (Yr[:, odd] + Y[:, odd]).abs().max() <= 1e-9 * scale


def test_not_available_for_nested():
    with pytest.raises(NotImplementedError):
        build_invariant_embedding(4, GroupSpec((tuple("abcde"),)),
                                  EmbeddingOptions(embedder="bursar_nested", root_normalise=True))
