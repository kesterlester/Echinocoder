"""
The build step: specification in, (evaluator module, metadata) out. All symbolic and combinatorial work happens here,
once; evaluation then only does tensor arithmetic.
"""
from __future__ import annotations

from .bursar import bursar_size
from .meta import EmbeddingMeta
from .module import InvariantEmbedding
from .slots import build_slots, iter_group_elements, signed_permutation_table
from .spec import EmbeddingOptions, GroupSpec


def build_invariant_embedding(k: int, group_spec: GroupSpec, options: EmbeddingOptions | None = None):
    """
    Build an embedding of n vectors in k dimensions that is invariant under G (and rotations) and separates orbits on
    every input.

    Parameters
    ----------
    k
        The dimension of each vector.
    group_spec
        The labelling group G (species, and whether parity is included). Also fixes n and the order of the vectors
        along the input's n axis.
    options
        Build options. Defaults to ``EmbeddingOptions()``: the unpruned Bursarial baseline evaluated in float32.

    Returns
    -------
    (InvariantEmbedding, EmbeddingMeta)
        The evaluator, mapping tensors of shape (..., n, k) to (..., D), and its metadata.
    """
    if options is None:
        options = EmbeddingOptions()
    if not isinstance(k, int) or k < 1:
        raise ValueError(f"k must be a positive integer, got {k!r}")

    slots = build_slots(group_spec, k)
    N, m = group_spec.order, len(slots)
    D = bursar_size(N, m)
    if D > options.max_outputs:
        raise ValueError(f"the unpruned Bursarial embedding for |G| = {N}, m = {m} would have D = {D} outputs, more "
                         f"than max_outputs = {options.max_outputs}. Unpruned embeddings are only feasible for small "
                         f"groups.")

    elements = tuple(iter_group_elements(group_spec))
    assert len(elements) == N and elements[0].is_identity()
    index, sign = signed_permutation_table(slots, elements)

    module = InvariantEmbedding(group_spec.n, k, slots, index, sign, options.dtype)
    meta = EmbeddingMeta(group_spec=group_spec, k=k, options=options, slots=slots, elements=elements,
                         index=index, sign=sign)
    assert module.output_shape.item() == meta.D == D
    return module, meta
