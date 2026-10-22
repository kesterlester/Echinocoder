"""
Beam-frame embeddings of collider events (notes, Section "Beam-frame embeddings for collider events").

The two incoming protons p, q are timelike with fixed masses and fixed p.q (that is, fixed sqrt(s)). In their centre-
of-mass frame, with p along +z (the *beam frame*), the Lorentz transformations that keep the set {p, q} fixed are the
rotations about z and the *swap*: a rotation by pi about any transverse axis, (px, py, pz) -> (px, -py, -pz), which
exchanges p and q. Up to rotations about z, the visible momenta are described exactly by

    w = (m_i, pt_i, i.j | pz_i, i^j),

with masses m_i (omitted for species of fixed common mass), transverse momenta pt_i, transverse dot products
i.j = px_i px_j + py_i py_j and wedge products i^j = px_i py_j - py_i px_j (i < j), and longitudinal momenta pz_i.
The swap fixes the first block and negates the second (the *bracket*). Permuting identical particles permutes the
slots, with a sign on the wedges. So the symmetry group acts on w by signed permutations, and the multiset machinery
applies unchanged. sqrt(s) is not part of w: it is the same for every event of a run, and is stored separately.

Optional extras: ``parity`` (blind to handedness: a reflection containing the beam axis, which negates the wedges and
nothing else, joins the group); ``met`` (a transverse-only object, for missing transverse momentum).
"""
from __future__ import annotations

import json
from dataclasses import asdict, dataclass
from itertools import combinations

import torch
from torch import nn, Tensor

from .bursar import bursar_coefficients, bursar_output_labels, bursar_size
from .slots import GroupElement, iter_group_elements
from .spec import EmbeddingOptions, GroupSpec

_KINDS = ("m", "pt", "dot", "pz", "wedge", "m2")
_MOMENTUM_DEGREE = {"m": 1, "m2": 2, "pt": 1, "pz": 1, "dot": 2, "wedge": 2}
BEAM_LABELS = ("p", "q")


@dataclass(frozen=True)
class BeamSpec:
    """
    The visible particles of a collider event, and which symmetries the embedding should respect.

    Parameters
    ----------
    species
        A tuple of species, each a tuple of distinct labels, e.g. ``(("j1", "j2", "j3"),)`` for three jets. Particles
        of one species are exchangeable. The flattened label order defines the object axis of the input.
    masses
        One entry per species: ``None`` if the masses of its particles vary from event to event (each particle then
        has a mass slot), or a number, the fixed common mass of its particles (no mass slots: it is a parameter of the
        embedding, like sqrt(s)). Default: all masses vary.
    parity
        Whether the embedding should also be blind to handedness: then a reflection containing the beam axis is a
        symmetry. It negates every wedge and nothing else.
    met
        Whether there is also one transverse-only object (missing transverse momentum), the last row of the input.
    mass_variable
        ``"m"`` (store masses) or ``"m2"`` (store masses squared). Either determines the event, for physical
        particles; see the notes for the trade-off.
    """
    species: tuple
    masses: tuple | None = None
    parity: bool = False
    met: bool = False
    mass_variable: str = "m"

    def __post_init__(self):
        GroupSpec(self.species)                                   # validates the species
        if set(self.labels) & set(BEAM_LABELS) or "met" in self.labels:
            raise ValueError(f"the labels {BEAM_LABELS + ('met',)} are reserved")
        if self.masses is not None and len(self.masses) != len(self.species):
            raise ValueError(f"masses needs one entry per species ({len(self.species)}), got {self.masses!r}")
        if self.mass_variable not in ("m", "m2"):
            raise ValueError(f"mass_variable must be 'm' or 'm2', got {self.mass_variable!r}")

    @property
    def labels(self) -> tuple:
        return tuple(label for s in self.species for label in s)

    @property
    def n_visible(self) -> int:
        return len(self.labels)

    @property
    def n_objects(self) -> int:
        """Rows of the input: the visible particles, then the missing transverse momentum if any."""
        return self.n_visible + int(self.met)

    @property
    def object_labels(self) -> tuple:
        return self.labels + (("met",) if self.met else ())

    def mass_varies(self, species_index: int) -> bool:
        return self.masses is None or self.masses[species_index] is None

    def group_spec(self, parity: bool | None = None) -> GroupSpec:
        """The group as a GroupSpec: the visible species, the beams (p, q) whose exchange is the swap, and parity."""
        return GroupSpec(self.species + (BEAM_LABELS,), parity=self.parity if parity is None else parity)

    def to_dict(self) -> dict:
        d = asdict(self)
        d["species"] = [list(s) for s in self.species]
        return d


@dataclass(frozen=True)
class BeamSlot:
    """One entry of w. ``positions`` are object indices; ``names`` the labels used to display it."""
    kind: str
    positions: tuple
    names: tuple

    def name(self, labels: tuple = ()) -> str:
        sep = "" if all(len(x) == 1 for x in self.names) else ","
        return f"{self.kind}_{sep.join(self.names)}"

    def latex(self, labels: tuple = ()) -> str:
        sep = "" if all(len(x) == 1 for x in self.names) else ","
        body = sep.join(x if x != "met" else r"\mathrm{met}" for x in self.names)
        if self.kind == "dot":
            return rf"({self.names[0]}\cdot {self.names[1]})_T".replace("met", r"\mathrm{met}")
        if self.kind == "wedge":
            return rf"({self.names[0]}\wedge {self.names[1]})_T".replace("met", r"\mathrm{met}")
        symbol = {"m": "m", "m2": "m^2", "pt": r"p_T", "pz": r"p_z"}[self.kind]
        if self.kind == "m2":
            return f"m_{{{body}}}^2"
        return f"{symbol}_{{{body}}}"

    def L_degree(self, k: int = 4) -> int:
        """Degree in the momenta (masses count 1, or 2 when masses squared are stored: see BeamW)."""
        return _MOMENTUM_DEGREE[self.kind]


def beam_slots(spec: BeamSpec, wedges: bool = True) -> tuple:
    """
    The slots of w, grouped so that each orbit of the symmetry group is a consecutive run: masses (species with
    varying masses), pt's, transverse dots, then the bracket: pz's and wedges.
    """
    blocks = spec.group_spec().species_index_blocks[:-1]
    names = spec.object_labels
    groups: list[list] = [[] for _ in range(5)]
    for b, block in enumerate(blocks):
        if spec.mass_varies(b):
            groups[0] += [BeamSlot(spec.mass_variable, (i,), (names[i],)) for i in block]
        groups[1] += [BeamSlot("pt", (i,), (names[i],)) for i in block]
        groups[3] += [BeamSlot("pz", (i,), (names[i],)) for i in block]
    if spec.met:
        groups[1].append(BeamSlot("pt", (spec.n_visible,), ("met",)))
    transverse = list(blocks) + ([(spec.n_visible,)] if spec.met else [])
    for a in range(len(transverse)):
        for b in range(a, len(transverse)):
            pairs = (list(combinations(transverse[a], 2)) if a == b
                     else [(i, j) for i in transverse[a] for j in transverse[b]])
            groups[2] += [BeamSlot("dot", p, (names[p[0]], names[p[1]])) for p in pairs]
            if wedges:
                groups[4] += [BeamSlot("wedge", p, (names[p[0]], names[p[1]])) for p in pairs]
    return tuple(s for g in groups for s in g)


def beam_signed_permutation_table(spec: BeamSpec, slots: tuple, elements: tuple) -> tuple[list, list]:
    """
    For each group element g and slot s: (g.w)[s] = sign[g][s] * w[index[g][s]]. ``elements`` act on the visible
    labels and on (p, q); exchanging p and q is the swap (negates pz's and wedges); the parity flag is the reflection
    (negates wedges). The missing transverse momentum is fixed by every element.
    """
    where = {(s.kind, s.positions): t for t, s in enumerate(slots)}
    nv = spec.n_visible
    index, sign = [], []
    for g in elements:
        perm = list(g.perm[:nv]) + ([nv] if spec.met else [])
        swap = len(g.perm) > nv and g.perm[nv] != nv              # without (p, q) in the group: never
        row_index, row_sign = [], []
        for s in slots:
            image = tuple(perm[i] for i in s.positions)
            sgn = 1
            if len(image) == 2 and image[0] > image[1]:
                image = (image[1], image[0])
                sgn = -1 if s.kind == "wedge" else 1
            if s.kind in ("pz", "wedge") and swap:
                sgn = -sgn
            if s.kind == "wedge" and g.parity:
                sgn = -sgn
            row_index.append(where[(s.kind, image)])
            row_sign.append(sgn)
        index.append(row_index)
        sign.append(row_sign)
    return index, sign


class _SafeSqrt(torch.autograd.Function):
    """sqrt of a non-negative tensor, with a finite derivative at 0 (the argument is clamped in the backward pass)."""

    @staticmethod
    def forward(ctx, x):
        y = x.clamp_min(0).sqrt()
        ctx.save_for_backward(x)
        return y

    @staticmethod
    def backward(ctx, grad):
        (x,) = ctx.saved_tensors
        return grad * 0.5 * x.clamp_min(torch.finfo(x.dtype).tiny).rsqrt()


class BeamW(nn.Module):
    """
    Computes w from beam-frame inputs of shape (..., n_objects, 4), columns (m, px, py, pz). For particles of fixed
    mass, and for the missing transverse momentum, the mass column is ignored; so is pz for the latter.
    """

    def __init__(self, spec: BeamSpec, slots: tuple, dtype: torch.dtype):
        super().__init__()
        self.spec = spec
        kind_code = {kind: c for c, kind in enumerate(_KINDS)}
        first = [s.positions[0] for s in slots]
        second = [s.positions[-1] for s in slots]
        self.register_buffer("kind", torch.tensor([kind_code[s.kind] for s in slots], dtype=torch.int64))
        self.register_buffer("first", torch.tensor(first, dtype=torch.int64))
        self.register_buffer("second", torch.tensor(second, dtype=torch.int64))
        self.m = len(slots)

    def forward(self, x: Tensor) -> Tensor:
        mass, px, py, pz = x[..., 0], x[..., 1], x[..., 2], x[..., 3]
        a, b = self.first, self.second
        pxa, pya, pxb, pyb = px[..., a], py[..., a], px[..., b], py[..., b]
        ma = mass[..., a]
        values = torch.stack([
            ma,                                                                  # m
            _SafeSqrt.apply(pxa * pxa + pya * pya),                              # pt
            pxa * pxb + pya * pyb,                                               # dot
            pz[..., a],                                                          # pz
            pxa * pyb - pya * pxb,                                               # wedge
            ma * ma,                                                             # m2
        ], dim=-1)
        return torch.gather(values, -1, self.kind.expand(*values.shape[:-2], self.m).unsqueeze(-1)).squeeze(-1)


class BeamEmbedding(nn.Module):
    """
    The evaluator: beam-frame inputs (..., n_objects, 4), columns (m, px, py, pz), to outputs (..., D). It computes w
    (BeamW) and hands it to a multiset embedding ``core`` (any evaluator with ``from_w``). Conventions follow IWPC's
    Encoding (integer buffers ``input_shape`` and ``output_shape``) without depending on IWPC.
    """

    def __init__(self, beam_w: BeamW, core: nn.Module, n_objects: int, dtype: torch.dtype):
        super().__init__()
        self.beam_w, self.core = beam_w, core
        self.register_buffer("input_shape", torch.tensor([n_objects, 4], dtype=torch.int64))
        self.register_buffer("output_shape", core.output_shape.clone().reshape(1))
        self.register_buffer("anchor", torch.zeros((), dtype=dtype))      # dtype/device anchor

    def forward(self, x) -> Tensor:
        if not isinstance(x, Tensor):
            x = torch.as_tensor(x, dtype=self.anchor.dtype, device=self.anchor.device)
        if tuple(x.shape[-2:]) != tuple(self.input_shape.tolist()):
            raise ValueError(f"expected input of shape (..., {self.input_shape[0].item()}, 4), got {tuple(x.shape)}")
        batch_shape = x.shape[:-2]
        w = self.beam_w(x).reshape(-1, self.beam_w.m)
        out = self.core.from_w(w)
        return out.reshape(*batch_shape, out.shape[-1])


class ProjectiveCore(nn.Module):
    """
    Option A, plain (a baseline): each element of the multiset over the permutations is (first block of g.w,
    RP(bracket of g.w)), where RP stores the bracket b up to sign as the coefficients of -q(x)^2, q(x) = sum_a b_a x^a
    (the Bursarial embedding of {b, -b}, Echinocoder's real-projective embedder: 2 len(b) - 1 numbers). The multiset
    of these vectors is embedded by a plain Bursarial embedding.
    """

    def __init__(self, index: list, sign: list, n_even: int, dtype: torch.dtype):
        super().__init__()
        self.register_buffer("g_index", torch.tensor(index, dtype=torch.int64))
        self.register_buffer("g_sign", torch.tensor(sign, dtype=dtype))
        self.n_even = n_even
        m = len(index[0])
        self.n_odd = m - n_even
        self.width = n_even + max(0, 2 * self.n_odd - 1)
        self.output_shape = torch.tensor([bursar_size(len(index), self.width)], dtype=torch.int64)

    def from_w(self, w: Tensor) -> Tensor:
        gw = w[:, self.g_index] * self.g_sign                                # (B, N, m)
        even, odd = gw[..., :self.n_even], gw[..., self.n_even:]
        if self.n_odd:
            B, N, r = odd.shape
            square = nn.functional.conv1d(odd.reshape(1, B * N, r), odd.reshape(B * N, 1, r).flip(-1),
                                          padding=r - 1, groups=B * N).reshape(B, N, 2 * r - 1)
            even = torch.cat([even, -square], dim=-1)
        return bursar_coefficients(even)


# ----- front ends: other input conventions -> beam-frame (m, px, py, pz) ----------------------------------------------

def beam_inputs_from_lab(four_momenta: Tensor) -> Tensor:
    """
    Lab-frame four-momenta (..., n, 4), columns (E, px, py, pz), with the beams along +-z and the lab frame their
    centre-of-mass frame, to beam-frame inputs (m, px, py, pz) with m = sqrt(max(E^2 - |p|^2, 0)).
    """
    E, p = four_momenta[..., :1], four_momenta[..., 1:]
    m2 = E * E - (p * p).sum(-1, keepdim=True)
    return torch.cat([_SafeSqrt.apply(m2), p], dim=-1)


def beam_inputs_from_pt_eta_phi_m(pt: Tensor, eta: Tensor, phi: Tensor, m: Tensor) -> Tensor:
    """Detector coordinates to beam-frame inputs (m, px, py, pz); each argument has shape (..., n)."""
    return torch.stack([m, pt * torch.cos(phi), pt * torch.sin(phi), pt * torch.sinh(eta)], dim=-1)


def _minkowski(a: Tensor, b: Tensor) -> Tensor:
    return a[..., 0] * b[..., 0] - (a[..., 1:] * b[..., 1:]).sum(-1)


def beam_inputs_from_beams(p: Tensor, q: Tensor, v: Tensor) -> tuple[Tensor, Tensor]:
    """
    Four-momenta in any frame (for example with a beam crossing angle): beams p, q of shape (..., 4) and visible
    particles v of shape (..., n, 4), columns (E, px, py, pz), metric (+, -, -, -). Returns the beam-frame inputs
    (m, px, py, pz), shape (..., n, 4), and sqrt(s), shape (...).

    The beam frame has time axis e_t = (p + q)/sqrt(s) and z axis e_z along the part of p orthogonal to e_t; x and y
    complete a right-handed orthonormal frame. The choice of x is arbitrary (embeddings do not depend on it); it is
    made by Gram-Schmidt from whichever lab axis is least aligned with e_z.
    """
    P = p + q
    s = _minkowski(P, P)
    root_s = s.clamp_min(0).sqrt()
    e_t = P / root_s.unsqueeze(-1)
    pz_part = p - _minkowski(p, e_t).unsqueeze(-1) * e_t
    e_z = pz_part / (-_minkowski(pz_part, pz_part)).clamp_min(0).sqrt().unsqueeze(-1)
    eye = torch.eye(4, dtype=p.dtype, device=p.device)
    candidates = eye[1:].expand(*e_z.shape[:-1], 3, 4)
    overlap = _minkowski(candidates, e_z.unsqueeze(-2)).abs()
    ref = torch.gather(candidates, -2, overlap.argmin(-1)[..., None, None].expand(*e_z.shape[:-1], 1, 4)).squeeze(-2)
    x = ref - _minkowski(ref, e_t).unsqueeze(-1) * e_t + _minkowski(ref, e_z).unsqueeze(-1) * e_z
    e_x = x / (-_minkowski(x, x)).clamp_min(0).sqrt().unsqueeze(-1)
    # e_y: orthogonal to e_t, e_x, e_z. Lowering its index gives the four-dimensional "cross product" of the three,
    # c_mu = (-1)^mu det(rows e_t, e_x, e_z without column mu); then raise, normalise, and orient so that the frame
    # (e_t, e_x, e_y, e_z) has positive determinant, like the lab frame.
    M = torch.stack([e_t, e_x, e_z], dim=-2)
    c = torch.stack([(-1) ** mu * torch.linalg.det(torch.cat([M[..., :mu], M[..., mu + 1:]], dim=-1))
                     for mu in range(4)], dim=-1)
    y = c * torch.tensor([1.0, -1.0, -1.0, -1.0], dtype=p.dtype, device=p.device)
    y = y / (-_minkowski(y, y)).clamp_min(0).sqrt().unsqueeze(-1)
    e_y = y * torch.sign(torch.linalg.det(torch.stack([e_t, e_x, y, e_z], dim=-2))).unsqueeze(-1)
    comps = [_minkowski(v, e.unsqueeze(-2)) for e in (e_t, e_x, e_y, e_z)]
    E, px, py, pz = comps[0], -comps[1], -comps[2], -comps[3]
    m = _SafeSqrt.apply(_minkowski(v, v))
    return torch.stack([m, px, py, pz], dim=-1), root_s


# ----- the build step -------------------------------------------------------------------------------------------------

CERTIFICATION_BEAM = (
    "Beam frame: w = (masses, pt's, transverse dots | pz's, transverse wedges) determines the visible momenta up to "
    "rotations about the beam axis, for every input with positive energies: each E_i follows from (m_i, pt_i, pz_i), "
    "and the transverse dots and wedges determine the transverse momenta up to one common rotation (SO(2) invariants "
    "separate orbits, including degenerate configurations). The remaining symmetries act on w by signed "
    "permutations: exchanging identical particles permutes slots (negating wedges whose order is reversed), and the "
    "swap of the beams p, q (a rotation by pi about a transverse axis) negates exactly the bracket (pz's and wedges). "
    "In the statements above about the multiset MS, read 'L' as w."
)
CERTIFICATION_WEDGES_DROPPED = (
    "Parity in the group, wedges omitted (Level 0): the transverse dots (with pt_i^2 = i.i) form the transverse Gram "
    "matrix, which determines the transverse momenta up to O(2), that is, up to a rotation about the beam axis and "
    "possibly a reflection containing it; with parity in the group that reflection is a symmetry, so nothing is lost. "
    "Parity then acts trivially on w and the multiset is taken over the rest of the group."
)
CERTIFICATION_PROJECTIVE = (
    "Level 0 (by construction), bracket stored up to sign: each member of the multiset over the permutations is "
    "(first block of g.w, coefficients of -q(x)^2 with q(x) = sum_a b_a x^a for the bracket b of g.w). The second part "
    "is the Bursarial embedding of {b, -b}, which determines b up to sign (unique factorisation in R[x]); the swap "
    "negates exactly b, so this removes it. The multiset of these vectors is embedded by the Bursarial embedding "
    "(injective on multisets)."
)


class BeamMeta:
    """Metadata for a beam-frame embedding: the beam specification, plus the metadata of the multiset stage."""

    def __init__(self, spec: BeamSpec, swap: str, wedges_dropped: bool, slots: tuple, inner=None, options=None,
                 elements: tuple = (), D: int | None = None):
        self.spec, self.swap, self.wedges_dropped, self.slots = spec, swap, wedges_dropped, slots
        self.inner, self.options, self.elements = inner, options, elements
        self._D = D

    def __getattr__(self, name):                        # delegate everything else to the multiset stage's metadata
        inner = self.__dict__.get("inner")
        if inner is None:
            raise AttributeError(name)
        return getattr(inner, name)

    @property
    def D(self) -> int:
        return self._D if self.inner is None else self.inner.D

    @property
    def m(self) -> int:
        return len(self.slots)

    @property
    def N(self) -> int:
        return len(self.elements)

    @property
    def certification_level(self) -> int:
        return 0 if self.inner is None else self.inner.certification_level

    def slot_names(self) -> list[str]:
        return [s.name() for s in self.slots]

    def output_names(self) -> list[str]:
        if self.inner is not None and callable(getattr(self.inner, "output_names", None)):
            return self.inner.output_names()
        if self.inner is not None:
            return [self.inner.base_name(b) for b in self.inner.kept_base_indices]
        return [f"E_{{{c},{j}}}" for c, j in bursar_output_labels(self.N, self._width)]

    @property
    def certification(self) -> str:
        text = CERTIFICATION_PROJECTIVE if self.inner is None else self.inner.certification
        text += " " + CERTIFICATION_BEAM
        if self.wedges_dropped:
            text += " " + CERTIFICATION_WEDGES_DROPPED
        return text

    def to_dict(self) -> dict:
        d = self.inner.to_dict() if self.inner is not None else {"D": self.D, "multiset_size": self.N}
        d["certification"] = self.certification
        d["beam"] = {"spec": self.spec.to_dict(), "swap": self.swap, "wedges_dropped": self.wedges_dropped,
                     "input": "beam-frame (m, px, py, pz) per object, shape (..., n_objects, 4)",
                     "root_s": "not an output: the same for every event of a run, stored separately",
                     "slots": [{"name": s.name(), "kind": s.kind, "momentum_degree": s.L_degree()} for s in self.slots]}
        d["outputs"] = [{"index": i, "name": nm} for i, nm in enumerate(self.output_names())]
        return d

    def to_json(self, **kwargs) -> str:
        return json.dumps(self.to_dict(), **kwargs)


def build_beam_embedding(spec: BeamSpec, options: EmbeddingOptions | None = None, swap: str = "group"):
    """
    Build an embedding of the visible momenta of a collider event, invariant under rotations about the beam axis,
    the exchange of the beams, exchanges of identical particles (and reflections containing the beam axis if
    ``spec.parity``), and separating everything else, on every input with positive energies.

    Parameters
    ----------
    spec
        The visible species, masses, parity and missing transverse momentum (BeamSpec).
    options
        The multiset stage's options (embedder, reductions, dtype, ...), as for build_invariant_embedding.
        ``root_normalise``, ``parity_split`` and ``metric`` are not used here.
    swap
        How the exchange of the beams is handled. ``"group"`` (default): it is one more element of the group, so the
        multiset is twice as large and the graded or nested packings can use it. ``"projective"``: the bracket is
        stored up to sign by Echinocoder's real-projective embedder: with embedder ``"bursar"``, before a plain
        Bursarial embedding of the multiset over the permutations (a baseline); with ``"bursar_cascade"``, as the
        first stage of a cascade, so that the permutations are still exploited by graded stages.

    Returns
    -------
    (BeamEmbedding, BeamMeta)
        The evaluator, mapping beam-frame inputs (..., n_objects, 4), columns (m, px, py, pz), to (..., D); and its
        metadata.
    """
    from .build import _build_from_table
    options = EmbeddingOptions() if options is None else options
    if options.root_normalise or options.parity_split or options.metric is not None:
        raise NotImplementedError("root_normalise, parity_split and metric are not available for beam embeddings")
    if swap not in ("group", "projective"):
        raise ValueError(f"swap must be 'group' or 'projective', got {swap!r}")
    wedges_dropped = options.reductions == "safe" and spec.parity
    slots = beam_slots(spec, wedges=not wedges_dropped)
    parity_enumerated = spec.parity and not wedges_dropped
    group_spec = spec.group_spec(parity_enumerated) if swap == "group" else GroupSpec(spec.species, parity_enumerated)
    N = group_spec.order
    if options.embedder == "bursar" and bursar_size(N, len(slots)) > options.max_outputs:
        raise ValueError(f"the Bursarial embedding would have {bursar_size(N, len(slots))} outputs, more than "
                         f"max_outputs = {options.max_outputs}")
    elements = tuple(iter_group_elements(group_spec))
    index, sign = beam_signed_permutation_table(spec, slots, elements)
    beam_w = BeamW(spec, slots, options.dtype)
    if swap == "projective" and options.embedder == "bursar_cascade":
        # Option A in its graded form: a cascade whose first stage is the subgroup generated by the swap. That stage
        # stores the bracket up to sign (as -Q(x)^2, Echinocoder's real-projective embedder), with the bracket's
        # coordinates placed so that the permutations still act linearly on what follows.
        group_spec = spec.group_spec(parity_enumerated)
        elements = tuple(iter_group_elements(group_spec))
        index, sign = beam_signed_permutation_table(spec, slots, elements)
        nv = spec.n_visible
        swap_index = next(i for i, g in enumerate(elements) if not g.parity and g.perm[nv] == nv + 1
                          and all(g.perm[j] == j for j in range(nv)))
        core, inner = _build_from_table(group_spec, 4, options, slots, elements, index, sign, False,
                                        group_spec.species_index_blocks, vector_w=False, cascade_first=swap_index)
        meta = BeamMeta(spec, swap, wedges_dropped, slots, inner=inner, options=options, elements=elements)
    elif swap == "projective":
        if options.embedder != "bursar":
            raise ValueError("swap='projective' needs embedder 'bursar' (plain baseline) or 'bursar_cascade'")
        n_even = sum(1 for s in slots if s.kind not in ("pz", "wedge"))
        core = ProjectiveCore(index, sign, n_even, options.dtype)
        meta = BeamMeta(spec, swap, wedges_dropped, slots, options=options, elements=elements,
                        D=int(core.output_shape.item()))
        meta._width = core.width
    else:
        core, inner = _build_from_table(group_spec, 4, options, slots, elements, index, sign, False,
                                        group_spec.species_index_blocks, vector_w=False)
        meta = BeamMeta(spec, swap, wedges_dropped, slots, inner=inner, options=options, elements=elements)
    module = BeamEmbedding(beam_w, core, spec.n_objects, options.dtype)
    assert module.output_shape.item() == meta.D
    return module, meta
