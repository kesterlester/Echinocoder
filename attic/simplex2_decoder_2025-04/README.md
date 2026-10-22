# A decoder for the simplex2 embedder, by Alex Hodson (April 2025)

Alex Hodson, a student, wrote a decoder for the simplicial-complex ("simplex2") embedder,
`C0HomDeg1_simplicialComplex_embedder_2_for_array_of_reals_as_multiset.py`, and emailed two versions
of it in April 2025. They were never committed at the time; they are kept here as a record.

- `simplexClass_v1.py`: the first version, exactly as sent. Alex described it as a working decoder,
  not yet stress-tested, usable only for n*k <= 6 or n <= 2: removing the n! multiplicity from the
  list of subdivided simplices took at least O(#simplices^2). He pointed to `Decoder.fix_lin_comb()`
  for the conditions used in the row-wise canonicalisation, which may help with a proof of
  bijectivity. The decoder relies on `Eji_LinComb._index` being used in the encoder's hash.
- `simplexClass_v2.py`: the updated version, exactly as sent, with comments and much shorter
  runtimes (the expensive removal of equivalent simplices is replaced by a cutoff).
- `simplexClass_v1_lester.py`: v1 with small edits by Christopher Lester: a `#!` line, input read
  from the command line, `debug=False`.
- `simplexClass_v2_lester.py`: v2 with small edits by Christopher Lester: debug prints commented out,
  `debug=False`.

The imports assume the files sit in the top folder of the repository. The code predates the
April 2025 refactor of the embedders: `Eji_LinComb` and the `distinct_permutations` module no longer
exist, so none of these files runs against current code.
