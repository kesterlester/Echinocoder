# Echinocoder

This is a library contains functions which are able to perform:

  * Small [topological embeddings](https://en.wikipedia.org/wiki/Embedding) of real [symmetric product spaces](https://en.wikipedia.org/wiki/Symmetric_product_(topology)), $SP^n(\mathbb R^k)$.  These are continuous bijective mappings of multisets of size $n$ containing vectors in $\mathbb{R}^k$ into $\mathbb R^N$ for some $N$.

  * Small [topological embeddings](https://en.wikipedia.org/wiki/Embedding) of $\mathbb R^k/Z_2$. These are spaces in which $\vec x\in\mathbb R^k$ is identified with $-\vec x$.  I am not sure what these are really supposed to be called. This libarary currently calls them [real projective spaces](https://en.wikipedia.org/wiki/Real_projective_space) but that might be an abuse of terminology.

In this context, the 'order' of the embedding is deemed to be the number of real values which are generated by encoding one set.  An encoder is deemed to have 'small order' or to be 'small' if the order is at most $O(nk \log n \log k)$.

Most embedders work only reals inputs and generate only real embeddings as that's the whole purpose of the libarary. However, some embedders will accept complex numbers as inputs and can generate complex numbers as outputs.  Where this is the case it is not always documented. Some of the embedders which can process complex inputs and outputs are nonetheless used (in complex mode) as steps in the implementation of other embedders.  The capacity for some embedders to process complex numbers such routines should be considered private (unexposed) even if technically visible. This is to allow interface standardisation.

## Embedders for $SP^n(\mathbb R^k)$ -- for multisets of vectors:

All these are (or should be) instances of [MultisetEmbedder](MultisetEmbedder.py).

### Embedder summaries:

| Method  | Order (leading) | Exact order (for $n>1$ and $k>1$) | Piecewise Linear | Infinitely Differentiable | Notes | Source |
|---------|-----------------|---------------|------------------|---------------------------|-------|--------|
| Simplex 1 | $O(nk)$         | $2nk+1$       | Yes              |  No   |       | [link](C0HomDeg1_simplicialComplex_embedder_1_for_array_of_reals_as_multiset.py) |
| Simplex 2 | $O(nk)$         | $2nk+1-k$       | Yes              |  No   |       | [link](C0HomDeg1_simplicialComplex_embedder_2_for_array_of_reals_as_multiset.py) |
| Dotting Conjecture | $O(nk\cdot \log n)$   | $n((k-1)(\lfloor{ \log_2 n }\rfloor+1)+1)$  | Yes   |   No   | Conjectured (but not yet proved!) to be an embedding. | [link](C0HomDeg1_conjectured_dotting_embedder_for_array_of_reals_as_multiset.py) |
| Dotting Overkill | $O(nk\cdot nk)$ |          | Yes              |  No   | Is provably an embedding, but is probably wasteful of resources. | | 
| Polynomial | $O(nk\cdot k)$    | $nk(k-1)$     | No               |  Yes  |       | [link](Cinf_numpy_polynomial_embedder_for_array_of_reals_as_multiset.py) |
| Bursarial  | $O(nk\cdot n)$    | $n + (k-1) n (n+1)/2$  | No      |  Yes  |       | [link](Cinf_sympy_bursar_embedder_for_array_of_reals_as_multiset.py) |
| Hybrid | $O(nk\cdot \sqrt{nk})$ | Minimum of Polynomial and Bursarial orders | No      |  Yes  | This method uses whichever of Polynomial or Bursarial has smallest order. | [link](Cinf_hybrid_embedder_for_array_of_reals_as_multiset.py) |

The orders (i.e. embedding sizes) quoted in the table above are for $n>1$ and $k>1$ only. For $n\le 1$ or $k\le 1$ each algorithm should fall back to an optimal embedding, i.e. one for which the exact order is $nk$.

### Further details:

* The [Simplicial Complex](https://en.wikipedia.org/wiki/Simplicial_complex) embedder works for any $n$ and $k$ and embeds into $2 n k+1$ reals.  The [simplical complex embedder sources](C0HomDeg1_simplicialComplex_embedder_1_for_array_of_reals_as_multiset.py) may be browsed.
* The [conjectured dotting embedder](C0HomDeg1_conjectured_dotting_embedder_for_array_of_reals_as_multiset.py) is based on the [dotting encoder](C0HomDeg1_dotting_encoder_for_array_of_reals_as_multiset.py). It is CONJECTURED (but not proved) to be an embedder. It has order $O(n k \log n)$. 
* The [polynomial embedder](Cinf_numpy_polynomial_embedder_for_array_of_reals_as_multiset.py)
has order $O(n k^2)$ in general, but happens to be optimal (i.e. embeds into $nk$ reals) for $k=1$ or $k=2$.
* The [(vanilla) Busarial embedder](Cinf_sympy_bursar_embedder_for_array_of_reals_as_multiset.py)
has order $O(n^2 k)$.  Its exact order is $n + (k-1) n (n+1)/2$. 
* The ['even' Busarial embedder](Cinf_sympy_bursar_embedder_for_array_of_reals_as_multiset.py)
has order $Binom(n+k,n)-1$. Although this embedder is very inefficient, its one possible benefit is that it does not treat any components in the $k$-space differently than any other. It is `even handed' (hence the name) w.r.t. the axes of the vectors. 
* The [Hybrid embedder](Cinf_hybrid_embedder_for_array_of_reals_as_multiset.py)  uses whichever of the Busarial or Polynomial embedders has the smaller order. The Hybrid consequently has order $O((nk)^{\frac 3 2})$.

## Embedders which work on $SP^m(\mathbb R)$ only.

These are similar to the other encoders, but they have $k=1$ hard coded into them. They expect inputs to be one-dimensional lists not two-dimensional arrays.

* The [pure sorting embedder](C0_sorting_embedder_for_list_of_reals_as_multiset.py)
is optimal (i.e. it embeds into $n$ reals for any $n$). It is also just a trivial sort! It is piecewise linear.
* There is a [polynomial embedder for lists (not arrays!) of reals](Cinf_numpy_polynomial_embedder_for_list_of_reals_as_multiset.py)). It is also optimal.  It's encoding is infinitely differentiable.

## Obsolete/Retired/Historical embedders:
* An early (nonlinear) [Simplicial Complex](https://en.wikipedia.org/wiki/Simplicial_complex) embedder 
([embedder source](Historical/C0_simplicialComplex_embedder_1_for_array_of_reals_as_multiset.py))
worked for any $n$ and $k$ and embedded into $4 n k+1$ reals. 
In principle the method could be re-written to embed into just $2 n k + 1$ reals but the current (less efficient) implementation choice was selected as it was expected to make outputs more stable.

## Embedders for what this library calls $RP(\mathbb R^k)$ ([real projective space](https://en.wikipedia.org/wiki/Real_projective_space)s):

* By setting $n=2$ and embedding the multiset $\left\\{\vec x,-\vec x\right\\}$ with $\vec x$ in $R^k$ one can use either of the Bursarial embedders to embed something this library calls $RP^k$ (which is possibly an abuse of the notation for real projective space of order $k$).  This $RP^k$ embedding would (for the (vanilla) Bursarial embedder) naively therefore be of size $2+(k-1)2(2+1)/2 = 2+3(k-1)$.  However, since all $k$ terms of order 1 in the auxiliary variable $y$ always disappear for multisets of this sort, the coefficients of those terms do not need to be recorded. This leaves only $2k-1$ reals needing to be recorded in the embedding for $RP^k$.  A method named [Cinf_numpy_regular_embedder_for_list_of_realsOrComplex_as_realOrComplexprojectivespace](Cinf_numpy_regular_embedder_for_list_of_realsOrComplex_as_realOrComplexprojectivespace.py) implements this method. It is order $2n-1$ when $n>0$.
* A small optimisation of the above method (implemented as [Cinf_numpy_complexPacked_embedder_for_list_of_reals_as_realprojectivespace](Cinf_numpy_complexPacked_embedder_for_list_of_reals_as_realprojectivespace.py))  reduces the order by one when $n>0$ and $n$ is even.


## Testing/examples

[example.py](example.py) is a simple example showing how some of the embedders could be used.

Various unit tests can be run with [pytest](https://docs.pytest.org).  It runs functions begining test_ in files with names matching test_XXXX.py

## References:

Neither the [Don Davis papers](https://www.lehigh.edu/~dmd1/toppapers.html) nor [Don Davis immersion list](https://www.lehigh.edu/~dmd1/imms.html) has been used to create this library. Both may, however, be useful references and sources of other references, so some are cached in the [DOCS](DOCS) directory.
