"""Tests for beam-frame embeddings of collider events."""
import json
import math

import pytest
import torch

from symembed import EmbeddingOptions, GroupSpec
from symembed.beam import (BeamSpec, beam_inputs_from_beams, beam_inputs_from_lab, beam_inputs_from_pt_eta_phi_m,
                           beam_slots, build_beam_embedding)

THREE = BeamSpec((("a", "b", "c"),))
MIXED = BeamSpec((("a", "b"), ("c", "d", "e")))
F64 = dict(dtype=torch.float64)


@pytest.fixture
def gen():
    return torch.Generator().manual_seed(2468)


def events(gen, n, batch=4):
    x = torch.randn(batch, n, 4, dtype=torch.float64, generator=gen)
    x[..., 0] = x[..., 0].abs()                                   # masses >= 0
    return x


def rotate(x, theta):                                             # about the beam axis
    c, s = math.cos(theta), math.sin(theta)
    y = x.clone()
    y[..., 1], y[..., 2] = c * x[..., 1] - s * x[..., 2], s * x[..., 1] + c * x[..., 2]
    return y


def flip(x, *columns):
    y = x.clone()
    for col in columns:
        y[..., col] = -y[..., col]
    return y


def close(a, b):
    return torch.allclose(a, b, rtol=1e-9, atol=1e-9 * max(1.0, a.abs().max().item()))


@pytest.mark.parametrize("spec,options,swap", [
    (THREE, EmbeddingOptions(embedder="bursar_graded", **F64), "group"),
    (THREE, EmbeddingOptions(embedder="bursar_graded", reductions="safe", **F64), "group"),
    (THREE, EmbeddingOptions(embedder="bursar", **F64), "group"),
    (THREE, EmbeddingOptions(**F64), "projective"),
    (MIXED, EmbeddingOptions(embedder="bursar_nested", **F64), "group"),
    (BeamSpec((("a", "b"),), met=True), EmbeddingOptions(embedder="bursar_graded", **F64), "group"),
    (BeamSpec((("a", "b", "c"),), masses=(0.0,)), EmbeddingOptions(embedder="bursar_graded", **F64), "group"),
])
def test_symmetries_and_non_symmetries(spec, options, swap, gen):
    embed, meta = build_beam_embedding(spec, options, swap=swap)
    x = events(gen, spec.n_objects)
    y = embed(x)
    assert y.shape == (4, meta.D)
    assert close(embed(rotate(x, 0.7)), y)                         # rotation about the beam axis
    assert close(embed(flip(x, 2, 3)), y)                          # beam swap: rotation by pi about x
    blocks = spec.group_spec().species_index_blocks[:-1]
    for block in blocks:
        if len(block) > 1:
            perm = list(range(spec.n_objects))
            perm[block[0]], perm[block[1]] = perm[block[1]], perm[block[0]]
            assert close(embed(x[:, perm]), y)                     # identical particles
    assert not close(embed(flip(x, 2)), y)                         # a reflection: not a symmetry without parity
    assert not close(embed(flip(x, 3)), y)                         # z -> -z alone: not a symmetry either
    if len(blocks) > 1:
        perm = list(range(spec.n_objects))
        perm[blocks[0][0]], perm[blocks[1][0]] = perm[blocks[1][0]], perm[blocks[0][0]]
        assert not close(embed(x[:, perm]), y)                     # different species are not exchangeable


@pytest.mark.parametrize("reductions", ["none", "safe"])
def test_parity(reductions, gen):
    spec = BeamSpec((("a", "b", "c"),), parity=True)
    embed, meta = build_beam_embedding(spec, EmbeddingOptions(embedder="bursar_graded", reductions=reductions, **F64))
    assert meta.wedges_dropped == (reductions == "safe")
    assert all(s.kind != "wedge" for s in meta.slots) == (reductions == "safe")
    x = events(gen, 3)
    y = embed(x)
    assert close(embed(flip(x, 2)), y) and close(embed(flip(x, 3)), y) and close(embed(flip(x, 2, 3)), y)


def test_fixed_masses_are_parameters(gen):
    spec = BeamSpec((("a", "b"), ("c",)), masses=(None, 0.105))
    embed, meta = build_beam_embedding(spec, EmbeddingOptions(embedder="bursar_graded", **F64))
    assert [s.name() for s in meta.slots if s.kind == "m"] == ["m_a", "m_b"]
    x = events(gen, 3)
    x2 = x.clone()
    x2[..., 2, 0] = 7.0                                            # the fixed-mass particle's mass column is ignored
    assert close(embed(x2), embed(x))


@pytest.mark.parametrize("spec,expected", [
    (THREE, 11), (BeamSpec((("a", "b", "c"),), masses=(0.0,)), 8), (MIXED, 19),
    (BeamSpec((("a", "b"),), met=True), 9), (BeamSpec((("a", "b", "c"),), parity=True), 11),
])
def test_local_injectivity(spec, expected, gen):
    # Generic evidence only (the proof is in the certification): the Jacobian with respect to the inputs that matter
    # has rank 4 per particle (3 if its mass is fixed) + 2 for the missing transverse momentum - 1 (the azimuth).
    for swap in ("group", "projective"):
        options = EmbeddingOptions(embedder="bursar_graded" if swap == "group" else "bursar", **F64)
        embed, _ = build_beam_embedding(spec, options, swap=swap)
        x = events(gen, spec.n_objects, batch=1)[0]
        J = torch.autograd.functional.jacobian(embed, x).reshape(-1, spec.n_objects * 4)
        J = J / J.norm(dim=1, keepdim=True).clamp_min(1e-300)       # outputs differ in scale by many orders
        assert torch.linalg.matrix_rank(J, rtol=1e-8).item() == expected


def test_mass_squared_option(gen):
    a, _ = build_beam_embedding(THREE, EmbeddingOptions(embedder="bursar_graded", **F64))
    b, meta = build_beam_embedding(BeamSpec(THREE.species, mass_variable="m2"),
                                   EmbeddingOptions(embedder="bursar_graded", **F64))
    assert [s.name() for s in meta.slots][:3] == ["m2_a", "m2_b", "m2_c"]
    x = events(gen, 3)
    assert not close(a(x), b(x))
    x2 = x.clone()
    x2[..., 0] = x[..., 0] ** 2
    assert close(a(x2), b(x))       # masses squared stored = masses stored, applied to the squared masses


def random_lorentz(gen):
    A = torch.randn(4, 4, dtype=torch.float64, generator=gen) * 0.3
    A = A - A.T
    eta = torch.diag(torch.tensor([1.0, -1.0, -1.0, -1.0], dtype=torch.float64))
    return torch.linalg.matrix_exp(eta @ A)


def test_front_ends_agree(gen):
    embed, _ = build_beam_embedding(THREE, EmbeddingOptions(embedder="bursar_graded", **F64))
    x = events(gen, 3, batch=2)
    E = torch.sqrt(x[..., 0] ** 2 + (x[..., 1:] ** 2).sum(-1))
    lab = torch.cat([E[..., None], x[..., 1:]], dim=-1)
    assert close(beam_inputs_from_lab(lab), x)
    pt, phi = torch.hypot(x[..., 1], x[..., 2]), torch.atan2(x[..., 2], x[..., 1])
    eta = torch.asinh(x[..., 3] / pt)
    assert close(beam_inputs_from_pt_eta_phi_m(pt, eta, phi, x[..., 0]), x)
    # Any frame: transform beams and particles by a random Lorentz transformation (a crossing angle is a special case)
    mp, Eb = 0.938, 6500.0
    Pb = math.sqrt(Eb ** 2 - mp ** 2)
    p = torch.tensor([Eb, 0.0, 0.0, Pb], dtype=torch.float64)
    q = torch.tensor([Eb, 0.0, 0.0, -Pb], dtype=torch.float64)
    Lam = random_lorentz(gen)
    x_frame, root_s = beam_inputs_from_beams((Lam @ p).expand(2, 4), (Lam @ q).expand(2, 4), lab @ Lam.T)
    assert torch.allclose(root_s, torch.full((2,), 2 * Eb, dtype=torch.float64))
    assert close(embed(x_frame), embed(x))
    assert torch.allclose(x_frame[..., 0], x[..., 0], atol=1e-6) and torch.allclose(x_frame[..., 3], x[..., 3])


def test_meta_and_counts():
    _, meta = build_beam_embedding(THREE, EmbeddingOptions(embedder="bursar_graded"))
    d = json.loads(meta.to_json())
    assert d["D"] == meta.D == 84 and d["beam"]["swap"] == "group" and "Beam frame" in d["certification"]
    assert [s.name() for s in beam_slots(THREE)] == [
        "m_a", "m_b", "m_c", "pt_a", "pt_b", "pt_c", "dot_ab", "dot_ac", "dot_bc",
        "pz_a", "pz_b", "pz_c", "wedge_ab", "wedge_ac", "wedge_bc"]
    assert d["slot_orbits"] == [[0, 1, 2], [3, 4, 5], [6, 7, 8], [9, 10, 11], [12, 13, 14]]
    _, proj = build_beam_embedding(THREE, swap="projective")
    assert proj.D == 405 and "up to sign" in proj.certification


def test_orbit_runs_matches_slot_orbits():
    from symembed.slots import build_slots, iter_group_elements, orbit_runs, signed_permutation_table, slot_orbits
    for spec in (GroupSpec((("a", "b", "c"), ("p", "q"))), GroupSpec((("a", "b"),), parity=True)):
        slots = build_slots(spec, 3)
        index, _ = signed_permutation_table(slots, tuple(iter_group_elements(spec)))
        assert orbit_runs(index) == slot_orbits(spec, slots)
