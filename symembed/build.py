"""
The build step: specification in, (evaluator module, metadata) out. All symbolic and combinatorial work happens here,
once; evaluation then only does tensor arithmetic.
"""
from __future__ import annotations

from .bursar import bursar_output_labels, bursar_size
from .involution import choose_involution, packed_forms, packing_from_signed_involution
from .meta import EmbeddingMeta
from .module import InvariantEmbedding
from .reductions import LinearRelation, choose_c_max, exact_expansion_forms, find_linear_relations
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

    # Guard before enumerating G: the packed embedding is at least half the size of the unpacked one; a graded one
    # can be much smaller, so for it only the size of G is limited.
    D_unpacked = bursar_size(N, m)
    estimate = {"bursar": D_unpacked, "bursar_involution": bursar_size(N, (m + 1) // 2), "bursar_graded": N}
    estimate = estimate[options.embedder]
    if estimate > options.max_outputs:
        raise ValueError(f"the Bursarial embedding for a multiset of {N} vectors of length m = {m} has about "
                         f"{estimate} coefficients, more than max_outputs = {options.max_outputs}. The Bursarial "
                         f"stage is only feasible for small groups.")

    elements = tuple(iter_group_elements(enumeration_spec))
    assert len(elements) == N and elements[0].is_identity()
    index, sign = signed_permutation_table(slots, elements)

    if options.embedder == "bursar_graded":
        from .graded import GradedEmbedding, GradedMeta, choose_packing
        graded = choose_packing(index, sign, elements, enumeration_spec.species_index_blocks)
        module = GradedEmbedding(group_spec.n, k, slots, index, sign, options.dtype, graded)
        meta = GradedMeta(group_spec, k, options, slots, elements, index, sign, graded, e_dropped)
        if meta.D > options.max_outputs:
            raise ValueError(f"the graded embedding has {meta.D} outputs, more than max_outputs = {options.max_outputs}")
        assert module.output_shape.item() == meta.D
        return module, meta

    # Layer 2: plain Bursarial, or involution packing (Level 0, real coefficients of a complex Bursarial embedding).
    packing, involution = None, None
    if options.embedder == "bursar_involution":
        involution = choose_involution(index, sign)
        if involution is None:
            raise ValueError("embedder='bursar_involution' needs an element of G acting on w as a non-trivial "
                             "involution, and there is none")
        packing = packing_from_signed_involution(index[involution], sign[involution])
    width = m if packing is None else packing.r
    D_base = bursar_size(N, width)
    if D_base > options.max_outputs:
        raise ValueError(f"the embedding has {D_base} coefficients, more than max_outputs = {options.max_outputs}")

    # Safe reduction (Level 1): drop outputs that are exact linear combinations of kept ones, for c <= c_max.
    relations: list[LinearRelation] = []
    c_max = 0
    if safe:
        c_max = options.linear_reduction_max_c
        if c_max is None:
            c_max = choose_c_max(N, m, options.linear_reduction_budget, width)
        c_max = min(c_max, N)
        offsets = {}
        for position, (c, _) in enumerate(bursar_output_labels(N, width)):
            offsets.setdefault(c, position)
        if packing is None:
            forms = [[[(index[g][jp], (sign[g][jp], 0))] for jp in range(m)] for g in range(N)]
        else:
            forms = packed_forms(index, sign, packing)
        for c, (_, re, im) in enumerate(exact_expansion_forms(forms, m, c_max), start=1):
            if im.any():
                raise AssertionError(f"involution packing produced a non-zero imaginary part at c = {c}")
            _, local_relations = find_linear_relations(re)
            off = offsets[c]
            relations.extend(LinearRelation(off + r.dropped, tuple((off + t, f) for t, f in r.terms))
                             for r in local_relations)
    dropped = {r.dropped for r in relations}
    kept = [i for i in range(D_base) if i not in dropped]

    module = InvariantEmbedding(group_spec.n, k, slots, index, sign, options.dtype, keep=kept, packing=packing)
    meta = EmbeddingMeta(group_spec=group_spec, k=k, options=options, slots=slots, elements=elements,
                         index=index, sign=sign, kept=kept, relations=tuple(relations), e_dropped=e_dropped,
                         linear_reduction_max_c=c_max, packing=packing, involution=involution)
    assert module.output_shape.item() == meta.D == len(kept)
    return module, meta
