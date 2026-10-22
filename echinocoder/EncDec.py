#!/usr/bin/env python

from collections.abc import Iterable
from fractions import Fraction
from itertools import pairwise
import numpy as np
import tools
from ndenumerate_slice import ndenumerate_fixed_axes_slice

from tools import numpy_array_of_frac_to_str
from tuple_ize import tuple_ize

class MonoLinComb:
    """A single term in a linear combination: one scalar coefficient
    times one basis matrix.

    Corresponds to a term such as α_s · L̂_(s) or δ_(s) · V_(s) in
    the Abel-summation notation of the simplicial embedder.  The basis
    vector is a k×n numpy array of the same shape as the input data.
    """

    def __init__(self, coeff, basis_vec):
        self.coeff = coeff
        self.basis_vec = basis_vec

    def __repr__(self):
        return f"MonoLinComb({repr(self.coeff)}, {repr(self.basis_vec)})"

    def __eq__(self, other):
        return self.coeff == other.coeff and np.array_equal(self.basis_vec, other.basis_vec)

    def __len__(self):
        return 1 # We are a mono (mono=1) lin comb after all.

class LinComb:
    """An ordered list of (coefficient, basis_matrix) terms.

    Terms are stored in the order they were added and are never
    auto-consolidated — two terms with the same basis vector remain
    separate entries.  This is intentional: the Abel-summation step
    relies on insertion order.  Call to_numpy_array() to evaluate the
    weighted sum back to a single matrix.
    """

    def __init__(self, initialiser=None):
        self.coeffs = []
        self.basis_vecs = []
        if initialiser is not None:
            self += initialiser

    def mlcs(self): # MLCS = Mono Lin CombinationS
        """Yield each term as a MonoLinComb (MLCS = Mono Lin CombinationS)."""
        return (MonoLinComb(c,b) for c,b in zip(self.coeffs, self.basis_vecs))

    def __len__(self):
        assert len(self.coeffs) == len(self.basis_vecs)
        return len(self.coeffs)

    def to_numpy_array(self):
        """Evaluate the linear combination and return the resulting matrix.

        Computes sum(coeff_i * basis_vec_i) over all terms.  Useful as a
        reconstruction check: for a correctly built LinComb this should
        reproduce the original input array (before canonicalisation).
        """
        ans = None
        first = True
        for c,b in zip(self.coeffs, self.basis_vecs):
            if first == True:
                ans = c*np.asarray(b)
                first = False
            else:
                ans += c*np.asarray(b)
        return ans

    def __add__(self, stuff):
        # Note: does not consolidate like terms — (3i+2j)+(5i) becomes (3i+2j+5i), not (8i+2j).
        # Consolidation is the caller's responsibility if needed.
        return LinComb((self, stuff))

    def __iadd__(self, stuff):
        #print(f"In iadd see stuff of type {stuff}")
        # Note that __add__ does not automatically consolidate. I.e. (3i+2j) + (5i) becomes (3i+2j+5i) not (8i+2j).
        # It is the user's responsibility to perform consolidation manually if they wish it to happen!!
        if isinstance(stuff, LinComb):
            self.coeffs.extend(stuff.coeffs)
            self.basis_vecs.extend(stuff.basis_vecs)
            return self

        if isinstance(stuff, MonoLinComb):
            self.coeffs.append(stuff.coeff)
            self.basis_vecs.append(stuff.basis_vec)
            return self

        if isinstance(stuff, Iterable):
            for elt in stuff:
                self += elt # Recurse!
            return self
        
        print(f"DDDD {stuff}")
        raise ValueError("LinComb.__iadd__ only knows how to add LimCombs and MonoLinCombs and iterables containing those.")

    def is_consolidated(self):
        """Return True if every basis vector in this combination is distinct.

        A consolidated LinComb has no repeated basis vectors; each Eji
        basis element appears at most once.  Used as a pre-condition check
        in barycentric_subdivide, because that function's sort-and-difference
        logic only makes sense on a deduplicated list.
        """
        # Need to convert 2D numpy arrays to tup(tup()) and 1D numpy arrays to tup() so that they are hashable:
        
        #print(f"Being asked for consolidation state of {self.basis_vecs}")
        basis_vecs_as_tuptups = [ tuple_ize(bv) for bv in self.basis_vecs ]
        #print(f"turned {self.basis_vecs} into {basis_vecs_as_tuptups}")
        return len(set(basis_vecs_as_tuptups)) == len(basis_vecs_as_tuptups)

    def __repr__(self):
        tmp = list(MonoLinComb(c,np.asarray(b)) for c,b in zip(self.coeffs, self.basis_vecs))
        return str(f"LinComb({tmp})")

    def __eq__(self, other):
        if len(self.coeffs) != len(other.coeffs):
            return False

        assert len(self.coeffs) == len(self.basis_vecs)
        assert len(other.coeffs) == len(other.basis_vecs)

        # Note that the order is required to match here, so eq means "same lin com in same order".
        for i,j in zip(self.mlcs(), other.mlcs()):
            if i != j:
               return False

        return True

def array_to_lin_comb(arr: np.array, fixed_axes=None, debug=False):
        """Expand a numpy array into a LinComb in the standard Eji basis.

        Each entry arr[i, j] becomes the term arr[i,j] · e^j_i, where e^j_i
        is the indicator matrix that is 1 at position (i, j) and 0 elsewhere.
        The result satisfies result.to_numpy_array() == arr.

        If fixed_axes is supplied (a dict mapping axis index to a coordinate
        value), only the entries along that coordinate slice are included.
        Algorithm 2 uses this to process one coordinate row at a time, e.g.
        fixed_axes={1: 0} extracts only the i=0 row contributions.
        """
        lin_comb = LinComb()
        if fixed_axes is None:
            locations = np.ndenumerate(arr)
        else:
            locations = ndenumerate_fixed_axes_slice(arr, fixed_axes=fixed_axes)

        for index, coeff in locations:
            #if debug:
            #    print(f"Considering pos {index} and coeff {coeff}.")
            basis_vec = np.zeros_like(arr)
            basis_vec[index] = 1
            lin_comb += MonoLinComb(coeff, basis_vec)
        return lin_comb

def barycentric_subdivide(lin_comb: LinComb, return_offset_separately=False, preserve_scale=True, debug=False, use_assertion_self_test=False):
    """Apply one Abel-summation (barycentric-subdivision) step to a LinComb.

    Sorts terms by coefficient descending, then replaces them with
    (gap, cumulative-prefix-sum) pairs.  In the document's notation this
    maps the (δ_(s), V_(s)) representation to the (Δ_(s), L_(s)) one:

        Δ_(s)  =  δ_(s) − δ_(s+1)      (consecutive gap differences)
        L_(s)  =  V_(0) + … + V_(s)     (cumulative prefix sums)

    If preserve_scale is True (default), prefix sums are normalised by
    their count, yielding L̂_(s) = L_(s)/(s+1), and the coefficients
    become α_s = (s+1)·Δ_(s).  This keeps every basis vector at
    one-norm 1 when the inputs have one-norm 1.

    If return_offset_separately is True, the minimum-coefficient term
    (the "anchor": min_coeff · sum_of_all_basis_vecs) is returned as a
    separate MonoLinComb rather than appended to the output LinComb.
    This separates the permutation-invariant offset (the row-minimum
    anchor) from the non-trivial difference terms.
    """
    
    if debug:
        print(f"lin_comb is\n{lin_comb}")

    if len(lin_comb) == 0:
        return LinComb()

    assert len(lin_comb) >= 1

    assert lin_comb.is_consolidated() # Note that LinComb does not consolidate elements, but our use case should only encounter conslolidated elements anyway! Or so we think. If that is not so, we need to know about it here!

    # Sort by coefficient in linear combination, big -> small
    # We need a list below as we will consume it twice when generating the diff_lin_comb
    sorted_lin_comb = sorted(zip(lin_comb.coeffs, lin_comb.basis_vecs), key=lambda x: x[0], reverse=True)

    if debug:
        print(f"sorted_lin_comb is\n{sorted_lin_comb}")

    coeffs = [ x for x, _ in sorted_lin_comb ]
    basis_vecs = [ x for _ , x in sorted_lin_comb ]

    if preserve_scale:
        diff_lin_comb = LinComb(MonoLinComb((fac := (i+1))*(x-y), sum(basis_vecs[:i+1], start=0*basis_vecs[0]+Fraction())/fac) for i, (x,y) in enumerate(pairwise(coeffs)))
        offset_mono_lin_comb = MonoLinComb((fac := len(basis_vecs))*coeffs[-1], sum(basis_vecs, start=0*basis_vecs[0]+Fraction())/fac)
    else:
        diff_lin_comb = LinComb(MonoLinComb((x-y), sum(basis_vecs[:i+1], start=0*basis_vecs[0])) for i, (x,y) in enumerate(pairwise(coeffs)))
        offset_mono_lin_comb = MonoLinComb(coeffs[-1], sum(basis_vecs, start=0*basis_vecs[0]))

    if debug:
        print(f"diff_lin_comb is\n{diff_lin_comb}")
        print(f"offset_mono_lin_comb is\n{offset_mono_lin_comb}")

    if return_offset_separately:
        ans = diff_lin_comb, offset_mono_lin_comb
    else:
        ans = diff_lin_comb + offset_mono_lin_comb

    """
    No need to do the following, but conceptually useful as documentation as it shows us one of the things
    that this subroute intends to achieve:
    """
    if use_assertion_self_test:
        assert np.allclose( lin_comb.to_numpy_array().astype(float), (diff_lin_comb + offset_mono_lin_comb).to_numpy_array().astype(float))

    if debug:
        print(f"About to return \n{ans}")

    return ans

def simplex_1_preprocess_steps(set_array : np.array,
                               preserve_scale_in_step_1=False,
                               preserve_scale_in_step_2=True,
                               canonicalise=False,
                               use_assertions=False,
                               debug=False):
    """Algorithm 1 preprocessing: global Eji expansion and two Abel summations.

    Implements the "flat" version of the simplicial embedder, treating the
    entire k×n input matrix as one list of Eji entries:

      Step 1 — Expand set_array into a LinComb in the Eji basis.
      Step 2 — First Abel summation: extract the global entry-minimum as
               the offset and compute first-difference (gap, prefix-sum)
               pairs.  The offset encodes the fully symmetric part of M.
      Step 3 — Second Abel summation: apply barycentric subdivision to
               the merged first-difference list, producing the final
               (α_s, L̂_(s)) coefficient/basis pairs.
      Step 4 — (if canonicalise=True) Sort the columns of every basis
               matrix lexicographically so that the encoding is
               column-permutation invariant.

    Returns (lin_comb, offset) where lin_comb holds the difference terms
    and offset holds the symmetric anchor.  If canonicalise=False, returns
    after Step 3 (useful for debugging without canonicalisation).

    preserve_scale_in_step_1 / preserve_scale_in_step_2 control whether
    each Abel summation normalises its prefix sums (see barycentric_subdivide).
    """

    """
    Step 1:

    Turn the array (which represents a set) into a linear combination of coefficients and Eji basis elements.
    Conceptually this step is turning:

       set_array = [[2,8],[4,5]]

    into

       lin_comb_0 = 2 * [[1,0],[0,0]] + 8 * [[0,1],[0,0]] + 4 * [[0,0],[1,0]] + 5 * [[0,0],[0,1]]

    """

    lin_comb_0 = array_to_lin_comb(set_array)

    if use_assertions:
        assert np.allclose(set_array.astype(float), lin_comb_0.to_numpy_array().astype(float))


    """
    Step 2:

    Re-write lin_comb_0 as a sum of the offset (which is the minimum coefficient times the perm-invariant
    basis all-ones element (e.g [[1,1],[1,1]]) and some (so called) differences.
    The latter are a set of non-negative coefficients times other basis vectors only
    containing zeros and ones.  Conceptually this step is turning:

       lin_comb_0 = 2 * [[1,0],[0,0]] + 8 * [[0,1],[0,0]] + 4 * [[0,0],[1,0]] + 5 * [[0,0],[0,1]]

    into

       lin_comb_1 + offset

    where the first differences are:
    
        lin_comb_1 = (8-5) * [[0,1],[0,0]]
                   + (5-4) * [[0,1],[0,1]]  
                   + (4-2) * [[0,1],[1,1]] 

    and where the offset is:

        offset = 2 * [[1,1],[1,1]]

    The above example assumed that preserve_scale=False is supplied to barycentric_subdivide, and that thus
    the one-norm of the basis vecs in the linear combination is growing as you go down the list, rather than
    constant as it would be if preserve_scale=True had been used instead.
    """

    lin_comb_1_first_diffs, offset = barycentric_subdivide(lin_comb_0, return_offset_separately=True, preserve_scale=preserve_scale_in_step_1, use_assertion_self_test=True)

    if use_assertions:
        assert np.allclose(set_array.astype(float), (lin_comb_1_first_diffs+offset).to_numpy_array().astype(float))

    """
    Step 3:

    Now we do a barycentric subdivision of lin_comb_1, storing the answer in lin_comb_2.  
    The purpose of this step is to make the resulting basis vectors sufficiently complicated that the
    process of canonicalise them will allow the set of all vertices to retain enough information that the
    canonicalisation process can be undone in all the materially important ways - according to PKH claim..
    [Provably the canonicalisation process would delete information if it were applied to lin_comb_1 directly.] 

    In principle this subdivision should be done with preserve_scale=True (as the vertices of mid-points of 
    simplex edges are things like (v1+v2)/2 not (v1+v2).  However, since this introduces a lot of fractions 
    into the output, a lot of debugging is done with preserve_scale=False.

    For our example input set above, this call turns:

        lin_comb_1 = 3 * [[0,1],[0,0]] +
                     1 * [[0,1],[0,1]] +
                     2 * [[0,1],[1,1]] 
    into

        lin_comb_2_second_diffs = (3-2) * [[0, 1], [0, 0]] +
                                  (2-1) * [[0, 2], [1, 1]] +
                                  (1-0) * [[0, 3], [1, 2]]
    """

    lin_comb_2_second_diffs = barycentric_subdivide(lin_comb_1_first_diffs, return_offset_separately=False, preserve_scale=preserve_scale_in_step_2, use_assertion_self_test=True)

    if use_assertions:
        assert np.allclose(set_array.astype(float), (lin_comb_2_second_diffs + offset).to_numpy_array().astype(float))


    """
    Step 3.9b:

    Option to return early, for debugging only:
    """

    if not canonicalise:
        return lin_comb_2_second_diffs, offset


    """
    Step 4:

    Canonicalise the basis vectors.

    We hope this step is a bijection (given the domain).  PKH claims it is, but I am suspicious. Claim is being tested!

    It re-orders the vectors within each set basis element so that the are listed in
    lexicographical order, and in doing so it ensures that we 
    our encoding is a set function.


    In our example event it turns

        lin_comb_2_second_diffs = (3-2) * [[0, 1], [0, 0]] +
                                  (2-1) * [[0, 2], [1, 1]] +
                                  (1-0) * [[0, 3], [1, 2]]

    into

        lin_comb_3_canonical = (3-2) * [[0, 0], [0, 1]] +    # ( because [0,1] is lexicographically AFTER  [0,0] )
                               (2-1) * [[0, 2], [1, 1]] +    # ( because [0,2] is lexicographically BEFORE [1,1] )
                               (1-0) * [[0, 3], [1, 2]]      # ( because [0,3] is lexicographically BEFORE [1,2] ).

    Note that after this step lin_comb_3_canonical+offset need not (in general) be the same as set_array so we do not have a line below saying

        assert np.allclose(set_array.astype(float), (lin_comb_3_canonical + offset).to_numpy_array().astype(float))

    even though we did have a line above saying

        assert np.allclose(set_array.astype(float), (lin_comb_2_second_diffs + offset).to_numpy_array().astype(float))
    .

    """

    lin_comb_3_canonical = LinComb(( MonoLinComb(coeff, tools.sort_np_array_rows_lexicographically(basis_vec)) for coeff, basis_vec in zip(lin_comb_2_second_diffs.coeffs, lin_comb_2_second_diffs.basis_vecs) ))

    """
    Step 5:

    We return both our canonicalised lin_comb_3, and our offsets (which are canonical by definition).
    """

    return lin_comb_3_canonical, offset

def simplex_2_preprocess_steps(set_array : np.array,
                               preserve_scale_in_step_1=False,
                               preserve_scale_in_step_2=True,
                               canonicalise=False,
                               use_assertions=False,
                               debug=False):
    """Algorithm 2 preprocessing: per-row Abel summation then global merge.

    Implements the "row-separated" version of the simplicial embedder.
    The key difference from Algorithm 1 is that the first Abel summation
    is applied independently to each coordinate row i, so the per-row
    minima m_i become separate anchor scalars:

      Step 1 — Decompose set_array into k per-row LinCombs, one for each
               coordinate row, each restricted to its own Eji slice.
      Step 2 — Apply a first Abel summation to each row separately,
               extracting per-row offsets (the m_i anchors) and per-row
               first-difference terms.  Merge all k per-row results into
               one flat LinComb.
      Step 3 — Second Abel summation on the merged list (same as
               Algorithm 1 from this point on).
      Step 4 — (if canonicalise=True) Lexicographic column sort.

    Returns (lin_comb, offsets) where offsets is a LinComb of k anchor
    terms (one per coordinate row).  Because each prefix element in
    Algorithm 2 contains Ejis from exactly one row, the resulting basis
    matrices have a cleaner structure than Algorithm 1.

    preserve_scale_in_step_1 / preserve_scale_in_step_2 control whether
    each Abel summation normalises its prefix sums (see barycentric_subdivide).
    """

    if debug:
        print(f"simplex_2_preprocess_steps was asked to encode {set_array}")

    n,k = set_array.shape

    """
    Step 1: 

    Turn the array (which represents a set) into a linear combination of coefficients and Eji basis elements, separated by component.
    Conceptually this step is turning:

       set_array = [[2,8],[4,5],[3,3]]

    into

       lin_comb_0[0] = 2 * [[1,0],[0,0],[0,0]] + 4 * [[0,0],[1,0],[0,0]] + 3 * [[0,0],[0,0],[1,0]]  # "x"-components
       lin_comb_0[1] = 8 * [[0,1],[0,0],[0,0]] + 5 * [[0,0],[0,1],[0,0]] + 3 * [[0,0],[0,0],[0,1]]  # "y"-components

    """

    lin_comb_0 = [ array_to_lin_comb(set_array, fixed_axes = {1:cpt_index}) for cpt_index in range(k) ]
    if debug:
        print(f"lin_comb_0 was {lin_comb_0}")
    if use_assertions:
        assert np.allclose(set_array.astype(float), LinComb(lin_comb_0).to_numpy_array().astype(float))

    """
    Step 2:

    Re-write each element of lin_comb_0 as a sum of an offset (which is the minimum coefficient times the
    perm-invariant all-ones element all (e.g. [[1,0],[1,0]] for the x-cpt) and some (so called) differences.
    The latter are a set of non-negative coefficients times other basis vectors only
    containing zeros and ones.  Conceptually this step is turning:

       lin_comb_0[0] = 2 * [[1,0],[0,0],[0,0]] + 4 * [[0,0],[1,0],[0,0]] + 3 * [[0,0],[0,0],[1,0]]  # "x"-components
       lin_comb_0[1] = 8 * [[0,1],[0,0],[0,0]] + 5 * [[0,0],[0,1],[0,0]] + 3 * [[0,0],[0,0],[0,1]]  # "y"-components

    into

       lin_comb_1 + offset

    where the first differences and offsets are:
    
        lin_comb_1[0]  = (4-3) * [[0,0],[1,0],[0,0]] + (3-2) * [[0,0],[1,0],[1,0]];   offsets[0]  = 2 * [[1,0],[1,0],[1,0]]   # and
        lin_comb_1[1]  = (8-5) * [[0,1],[0,0],[0,0]] + (5-3) * [[0,1],[0,1],[0,0]];   offsets[1]  = 3 * [[0,1],[0,1],[0,1]] 

    which are flattened into

        lin_comb_1  = (4-3) * [[0,0],[1,0],[0,0]] +
                      (3-2) * [[0,0],[1,0],[1,0]] +
                      (8-5) * [[0,1],[0,0],[0,0]] +
                      (5-3) * [[0,1],[0,1],[0,0]] 
    and

        offsets  = 2 * [[1,0],[1,0],[1,0]] +
                   3 * [[0,1],[0,1],[0,1]]

    The above example assumed that preserve_scale=False is supplied to barycentric_subdivide, and that thus
    the one-norm of the basis vecs in the linear combination is growing as you go down the list, rather than
    constant as it would be if preserve_scale=True had been used instead.
    """
    
    lin_comb_1_first_diffs = [None] * k
    offsets = [None] *k

    assert len(lin_comb_0) == k
    assert len(lin_comb_1_first_diffs) == k
    assert len(offsets) == k

    for i in range(k):
        lin_comb_1_first_diffs[i], offsets[i] = barycentric_subdivide(lin_comb_0[i], return_offset_separately=True, preserve_scale=preserve_scale_in_step_1, use_assertion_self_test=True)
        if debug:
            print(f"lin_comb_1[{i}] was {lin_comb_1_first_diffs[i]}")
            print(f"offsets[{i}] was {offsets[i]}")
            print()

    # Here is the merge phase:
    lin_comb_1_first_diffs = LinComb(lin_comb_1_first_diffs)
    offsets = LinComb(offsets)

    if debug:
        print(f"After the merge phase lin_comb_1 was {lin_comb_1_first_diffs}")
        print(f"offsets was {offsets}")
        print()

    if use_assertions:
        assert np.allclose(set_array.astype(float), (lin_comb_1_first_diffs+offsets).to_numpy_array().astype(float))
        if debug:
            print("Happy after step 2")

    """
    Step 3:

    Now we do a barycentric subdivision of lin_comb_1, storing the answer in lin_comb_2.  
    The purpose of this step is to make the resulting basis vectors sufficiently complicated that the
    process of canonicalise them will allow the set of all vertices to retain enough information that the
    canonicalisation process can be undone in all the materially important ways - according to CGL guess.
    [Provably the canonicalisation process would delete information if it were applied to lin_comb_1 directly.] 

    In principle this subdivision should be done with preserve_scale=True (as the vertices of mid-points of 
    simplex edges are things like (v1+v2)/2 not (v1+v2).  However, since this introduces a lot of fractions 
    into the output, a lot of debugging is done with preserve_scale=False.

    For our example input set above, this call turns:

        lin_comb_1  = 1 * [[0,0],[1,0],[0,0]] +
                      1 * [[0,0],[1,0],[1,0]] +
                      3 * [[0,1],[0,0],[0,0]] +
                      2 * [[0,1],[0,1],[0,0]]

    into

        lin_comb_2_second_diffs = (3-2) * [[0,1],[0,0],[0,0]]
                                  (2-1) * [[0,2],[0,1],[0,0]]
                                  (1-1) * [[0,2],[1,1],[0,0]]
                                  (1-0) * [[0,2],[2,1],[1,0]]

    without affecting the offsets:

        offsets  = 2 * [[1,0],[1,0],[1,0]] +
                   3 * [[0,1],[0,1],[0,1]]

    """

    lin_comb_2_second_diffs = barycentric_subdivide(lin_comb_1_first_diffs, return_offset_separately=False, preserve_scale=preserve_scale_in_step_2, use_assertion_self_test=True)

    if use_assertions:
        assert np.allclose(set_array.astype(float), (lin_comb_2_second_diffs + offsets).to_numpy_array().astype(float))
        if debug:
            print("Happy after step 3")

    """
    Step 3.9b:

    Option to return early, for debugging only:
    """

    if not canonicalise:
        return lin_comb_2_second_diffs, offsets


    """
    Step 4:

    Canonicalise the basis vectors.

    We believe this step is a bijection (given the domain) --- but this has not been proved.

    For the exmple above it changes

        lin_comb_2_second_diffs = 1 * [[0,1],[0,0],[0,0]] +
                                  1 * [[0,2],[0,1],[0,0]] +
                                  0 * [[0,2],[1,1],[0,0]] +
                                  1 * [[0,2],[2,1],[1,0]]
    to

        lin_comb_3_canonical =  1 * [[0,0],[0,0],[0,1]] +     # ( because [0,0] = [0,0] < [0,1] )
                                1 * [[0,0],[0,1],[0,2]] +     # ( because [0,0] < [0,1] < [0,2] )
                                0 * [[0,0],[0,2],[1,1]] +     # ( because [0,0] < [0,2] < [1,1] )
                                1 * [[0,2],[1,0],[2,1]]       # ( because [0,2] < [1,0] < [2,1] )

    and the offsets are not changed as they are already canonical by construction.

    """

    lin_comb_3_canonical = LinComb(( MonoLinComb(coeff, tools.sort_np_array_rows_lexicographically(basis_vec)) for coeff, basis_vec in zip(lin_comb_2_second_diffs.coeffs, lin_comb_2_second_diffs.basis_vecs) ))

    return lin_comb_3_canonical, offsets




    

