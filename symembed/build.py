"""
The build step: specification in, (evaluator module, metadata) out. All symbolic and combinatorial work happens here,
once; evaluation then only does tensor arithmetic.
"""
from __future__ import annotations

from .bursar import bursar_output_labels, bursar_size
from .meta import EmbeddingMeta
from .module import InvariantEmbedding
from .reductions import LinearRelation, choose_c_max, exact_expansion, find_linear_relations
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
    safe = options.reductions == "safe"

    # Safe reduction (Level 0 for the smaller w): with parity in G, omit the determinants. Parity then acts trivially
    # on w, so the multiset is taken over the permutation part of G only.
    e_dropped = safe and group_spec.parity
    slots = build_slots(group_spec, k)
    enumeration_spec = group_spec
    if e_dropped:
        slots = tuple(s for s in slots if s.kind == "d")
        enumeration_spec = GroupSpec(group_spec.species, parity=False)

    N, m = enumeration_spec.order, len(slots)
    D_base = bursar_size(N, m)
    if D_base > options.max_outputs:
        raise ValueError(f"the Bursarial embedding for a multiset of {N} vectors of length m = {m} has {D_base} "
                         f"coefficients, more than max_outputs = {options.max_outputs}. The Bursarial stage is only "
                         f"feasible for small groups.")

    elements = tuple(iter_group_elements(enumeration_spec))
    assert len(elements) == N and elements[0].is_identity()
    index, sign = signed_permutation_table(slots, elements)

    # Safe reduction (Level 1): drop outputs that are exact linear combinations of kept ones, for c <= c_max.
    relations: list[LinearRelation] = []
    c_max = 0
    if safe:
        c_max = options.linear_reduction_max_c
        if c_max is None:
            c_max = choose_c_max(N, m, options.linear_reduction_budget)
        c_max = min(c_max, N)
        offsets = {}
        position = 0
        for c, j in bursar_output_labels(N, m):
            offsets.setdefault(c, position)
            position += 1
        for c, (_, coeffs) in enumerate(exact_expansion(index, sign, m, c_max), start=1):
            _, local_relations = find_linear_relations(coeffs)
            off = offsets[c]
            relations.extend(LinearRelation(off + r.dropped, tuple((off + t, f) for t, f in r.terms))
                             for r in local_relations)
    dropped = {r.dropped for r in relations}
    kept = [i for i in range(D_base) if i not in dropped]

    module = InvariantEmbedding(group_spec.n, k, slots, index, sign, options.dtype, keep=kept)
    meta = EmbeddingMeta(group_spec=group_spec, k=k, options=options, slots=slots, elements=elements,
                         index=index, sign=sign, kept=kept, relations=tuple(relations), e_dropped=e_dropped,
                         linear_reduction_max_c=c_max)
    assert module.output_shape.item() == meta.D == len(kept)
    return module, meta
