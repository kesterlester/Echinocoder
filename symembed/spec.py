"""
User-facing specification objects: which group, and which build options.

See Section 5 of the invariant_multisets notes for the interface these
implement.
"""
from __future__ import annotations

from dataclasses import dataclass, asdict
from math import factorial, prod

import torch


@dataclass(frozen=True)
class GroupSpec:
    """
    The labelling symmetry group G, as a direct product of symmetric groups (one per species), optionally times the
    parity group P = {1, p}.

    Parameters
    ----------
    species
        A tuple of species, each a tuple of distinct labels, e.g. ``(("a", "b"), ("p",))`` for S(2) x S(1). The
        flattened label order defines the ``n`` axis of the input: ``X[..., i, :]`` is the vector carrying the i-th
        label of ``labels``.
    parity
        Whether G contains the parity element p, which reflects every vector (flips the sign of its last component).
        p leaves every dot product unchanged and negates every determinant.
    """
    species: tuple
    parity: bool = False

    def __post_init__(self):
        if not isinstance(self.species, tuple) or not all(isinstance(s, tuple) for s in self.species):
            raise TypeError(f"species must be a tuple of tuples of labels, got {self.species!r}")
        if any(len(s) == 0 for s in self.species):
            raise ValueError("every species must contain at least one label")
        labels = self.labels
        if len(labels) == 0:
            raise ValueError("at least one label is required")
        if len(set(labels)) != len(labels):
            raise ValueError(f"labels must be distinct across all species, got {labels!r}")
        if not isinstance(self.parity, bool):
            raise TypeError(f"parity must be a bool, got {self.parity!r}")

    @property
    def labels(self) -> tuple:
        """All labels, species by species, in the order of the input's n axis."""
        return tuple(label for s in self.species for label in s)

    @property
    def n(self) -> int:
        """The number of vectors."""
        return len(self.labels)

    @property
    def species_index_blocks(self) -> tuple:
        """For each species, the tuple of its positions along the n axis."""
        blocks, start = [], 0
        for s in self.species:
            blocks.append(tuple(range(start, start + len(s))))
            start += len(s)
        return tuple(blocks)

    @property
    def order(self) -> int:
        """|G|."""
        return prod(factorial(len(s)) for s in self.species) * (2 if self.parity else 1)

    def to_dict(self) -> dict:
        return {"species": [list(s) for s in self.species], "parity": self.parity}


_SUPPORTED = {
    "embedder": ("bursar", "bursar_involution", "bursar_graded", "bursar_nested"),
    "reductions": ("none", "safe"),
    "field": ("real",),
}


@dataclass(frozen=True)
class EmbeddingOptions:
    """
    Options controlling how the embedding is built.

    Parameters
    ----------
    embedder
        The multiset embedding used in layer 2. ``"bursar"``: the Bursarial embedding of {g.w}.
        ``"bursar_involution"``: involution packing (notes, Section 3.3). An involution g0 in G (the one with the
        smallest |trace| on w) is used to pack w into complex coordinates on which g0 acts as conjugation, so that
        the complex Bursarial polynomial has real coefficients; roughly halves the number of outputs. Level 0.
        ``"bursar_graded"``: graded packing (notes, Section 3.4). A dihedral subgroup <h, s> of G is used: coordinates
        placed at power x^j are scaled by exp(2 pi i j / r) under h (so every coefficient E_{c,j} with r not dividing
        j vanishes) and conjugated by s (so the rest are real); G-invariant coordinates may be output directly. The
        (h, s) with fewest outputs is chosen automatically. Subsumes "bursar_involution" (h = identity). Level 0.
        ``"bursar_nested"``: nested coset packing (notes, Section 3.5): graded packing on each H-orbit for a dihedral
        subgroup H, then a Bursarial embedding over the cosets. Used only if it beats "bursar_graded" (roughly by 2
        when there are many cosets, e.g. S(5)); otherwise the graded embedding is returned. Level 0.
    reductions
        ``"none"``: the unpruned baseline (Level 0). ``"safe"``: the reductions that need no judgement (Level 0/1):
        omit the determinants when parity is in G, and drop outputs that are exact linear combinations of kept
        outputs (checked for c <= ``linear_reduction_max_c``).
    homogenise
        Split outputs into parts homogeneous in L. Not implemented yet.
    root_normalise
        Replace each homogeneous output f of degree r by sign(f)|f|^(1/r). Not implemented yet.
    field
        Only ``"real"`` is supported.
    dtype
        Floating-point type used by the evaluator. The build itself is exact.
    max_outputs
        Refuse to build an embedding whose unpruned Bursarial stage has more outputs than this. A guard against
        accidentally requesting an embedding that could never be evaluated.
    metric
        The diagonal of the metric used for the dot products d, as a tuple of k entries +1 or -1; None (the default)
        means Euclidean. For example (1, -1, -1, -1) for Minkowski space. The determinants e do not depend on it, and
        nothing downstream of w does either (see the notes, Section 1, on indefinite metrics).
    linear_reduction_max_c
        With ``reductions="safe"``: the largest c for which the E_{c,j} are expanded exactly and checked for linear
        dependencies. ``None`` chooses it automatically from ``linear_reduction_budget``.
    linear_reduction_budget
        With ``reductions="safe"`` and ``linear_reduction_max_c=None``: the maximum number of integer coefficients
        in the exact expansions.
    """
    embedder: str = "bursar"
    reductions: str = "none"
    homogenise: bool = False
    root_normalise: bool = False
    field: str = "real"
    dtype: torch.dtype = torch.float32
    max_outputs: int = 10_000_000
    metric: tuple | None = None
    linear_reduction_max_c: int | None = None
    linear_reduction_budget: int = 2_000_000

    def __post_init__(self):
        for name, allowed in _SUPPORTED.items():
            value = getattr(self, name)
            if value not in allowed:
                raise NotImplementedError(f"{name}={value!r} is not implemented; supported: {allowed}")
        if self.homogenise:
            raise NotImplementedError("homogenise=True is not implemented yet")
        if self.root_normalise:
            raise NotImplementedError("root_normalise=True is not implemented yet")
        if self.linear_reduction_max_c is not None and self.linear_reduction_max_c < 0:
            raise ValueError(f"linear_reduction_max_c must be non-negative, got {self.linear_reduction_max_c}")
        if self.metric is not None and (not isinstance(self.metric, tuple) or any(x not in (1, -1) for x in self.metric)):
            raise ValueError(f"metric must be a tuple of +1/-1 entries, got {self.metric!r}")
        if self.dtype not in (torch.float32, torch.float64):
            raise ValueError(f"dtype must be torch.float32 or torch.float64, got {self.dtype!r}")

    def to_dict(self) -> dict:
        d = asdict(self)
        d["dtype"] = str(self.dtype)
        return d
