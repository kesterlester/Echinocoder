"""Outputs of definite parity (option parity_split)."""
import pytest
import torch

from symembed import EmbeddingOptions, GroupSpec, build_invariant_embedding
from symembed.parity import ParitySplitMeta

S3S2 = GroupSpec((("a", "b", "c"), ("p", "q")))


def options(**kw):
    return EmbeddingOptions(dtype=torch.float64, parity_split=True, **kw)


@pytest.fixture
def gen():
    return torch.Generator().manual_seed(8642)


@pytest.mark.parametrize("embedder,k", [("bursar_graded", 3), ("bursar_graded", 4), ("bursar_involution", 3),
                                        ("bursar", 2)])
def test_outputs_have_definite_parity(embedder, k, gen):
    spec = S3S2 if embedder != "bursar" else GroupSpec((("a", "b"), ("p",)))
    embed, meta = build_invariant_embedding(k, spec, options(embedder=embedder))
    assert isinstance(meta, ParitySplitMeta)
    X = torch.randn(5, spec.n, k, dtype=torch.float64, generator=gen)
    reflection = torch.ones(k, dtype=torch.float64)
    reflection[-1] = -1.0
    Y, Yr = embed(X), embed(X * reflection)
    names = meta.output_names()
    odd = [i for i, nm in enumerate(names) if nm.startswith("odd part") or nm.endswith("[odd]")]
    even = [i for i in range(len(names)) if i not in odd]
    scale = Y.abs().max()
    assert (Yr[:, even] - Y[:, even]).abs().max() <= 1e-9 * scale
    assert (Yr[:, odd] + Y[:, odd]).abs().max() <= 1e-9 * scale


def test_parts_add_up_to_the_original(gen):
    k = 3
    split, meta = build_invariant_embedding(k, S3S2, options(embedder="bursar_graded", reductions="safe"))
    plain, _ = build_invariant_embedding(k, S3S2, EmbeddingOptions(embedder="bursar_graded", reductions="safe",
                                                                   dtype=torch.float64))
    X = torch.randn(3, 5, k, dtype=torch.float64, generator=gen)
    Ys, Yo = split(X), plain(X)
    position = 0
    for i, ps in enumerate(meta.parities):
        if len(ps) <= 1:
            torch.testing.assert_close(Ys[:, position], Yo[:, i])
            position += 1
        else:
            torch.testing.assert_close(Ys[:, position] + Ys[:, position + 1], Yo[:, i])
            position += 2
    assert position == meta.D == 307


def test_counts():
    _, m3 = build_invariant_embedding(3, S3S2, options(embedder="bursar_graded", reductions="safe"))
    _, m4 = build_invariant_embedding(4, S3S2, options(embedder="bursar_graded", reductions="safe"))
    assert m3.to_dict()["parity_split"] == {"split": 143, "pure_even": 21, "pure_odd": 0}
    assert m4.D == 184                     # = 100 + 84: every mixed output is split, no pure one


def test_noop_with_parity_in_g():
    spec = GroupSpec(S3S2.species, parity=True)
    _, meta = build_invariant_embedding(4, spec, options(embedder="bursar_graded", reductions="safe"))
    assert not isinstance(meta, ParitySplitMeta) and meta.D == 99


def test_not_available_for_nested():
    spec = GroupSpec(S3S2.species, parity=False)
    with pytest.raises(NotImplementedError):
        build_invariant_embedding(4, GroupSpec((tuple("abcde"),)), options(embedder="bursar_nested"))
    assert spec is not None


@pytest.mark.parametrize("k", [2, 3, 4])
def test_declared_purity_holds_output_by_output(k, gen):
    # Includes the outputs found pure through the conjugation mode (e.g. E_{c,0} when position 0 holds
    # (d's) + i (e's), on which parity acts as complex conjugation).
    from symembed.parity import output_parities
    embed, meta = build_invariant_embedding(k, S3S2, EmbeddingOptions(embedder="bursar_graded", reductions="safe",
                                                                      dtype=torch.float64))
    parities = output_parities(meta)
    X = torch.randn(6, 5, k, dtype=torch.float64, generator=gen)
    reflection = torch.ones(k, dtype=torch.float64)
    reflection[-1] = -1.0
    Y, Yr = embed(X), embed(X * reflection)
    scale = Y.abs().max(dim=0).values
    for i, ps in enumerate(parities):
        if len(ps) == 1:
            sign = 1.0 if ps == {0} else -1.0
            assert (Yr[:, i] - sign * Y[:, i]).abs().max() <= 1e-9 * scale[i]
    if k == 4:
        assert sum(len(ps) == 1 for ps in parities) == 16
