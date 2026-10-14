"""
symembed: continuous, algebraic embeddings of n vectors in k dimensions that are invariant under a labelling group G
(a product of symmetric groups, optionally with parity) and under rotations, and that separate orbits on every input.

See the invariant_multisets notes for the mathematics and the interface. Only the unpruned Bursarial baseline
(certification Level 0) is implemented so far.

Typical use::

    from symembed import GroupSpec, EmbeddingOptions, build_invariant_embedding

    embed, meta = build_invariant_embedding(k=2, group_spec=GroupSpec(species=(("a", "b"), ("p",))))
    Y = embed(X)        # X.shape == (..., 3, 2)  ->  Y.shape == (..., meta.D)
"""
from .spec import GroupSpec, EmbeddingOptions
from .build import build_invariant_embedding
from .module import InvariantEmbedding, FlattenedInvariantEmbedding
from .meta import EmbeddingMeta

__all__ = ["GroupSpec", "EmbeddingOptions", "build_invariant_embedding", "InvariantEmbedding",
           "FlattenedInvariantEmbedding", "EmbeddingMeta"]
