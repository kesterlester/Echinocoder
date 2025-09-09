import sympy as sp
from functools import partial
from vertex_matches import generate_viable_vertex_match_matrices
import sympy_tools as spt
import decider_functions.decider_function as df
from Match_Tracker import Match_Tracker
from collections import namedtuple

class Rational_Decider:

    """

    Note that for any match matrix, our "job" is to
    work out whether there exists non-collapsing solutions for it 
    (i.e. collections of blue and red dots that have not entirely 
    cancelled with each other).  Furthermore, if we can find multiple
    sets of not-fully cancelling blue and red dots, then we are
    interested in the sets among those that ARE SMALLEST.

    Note, therefore, the names of the methods:

        matrix_does_not_collapse(L_matrix), and
        matrix_collapses(L_matrix)

    defined below are perhaps badly named, because they may give the
    false impression that a certain L_matrix either always or never
    causes collapse. This is a false impression as whether collapse
    would be caused is actually a function of where the bats are.

    Any time that, given bats B, we find a non-collapsing solution
    with P points remaining, we know that the minimal confusable set
    size for this M and k is at most P.

    However:

        (1) given a P>0 for some B we don't immediately know if a
        better bat position B2 might lead to an even better bound,
        i.a. a P2(B2) for which 0 < P2(B2) < P(B).

        and

        (2) if we spot a collapse for some B (i.e. P(B)=0)
        we do not immediately know that our bat position was just
        poorly chosen for this match matrix.
        
    Both of threse issues are potentially addressable.

    Issue (1) would be addressed if it were possible to show that all
    nonzero P agree for a given match matrix.  We do not currently
    know if this is true, though.  AFTERTHOUGHT: We know that a non-
    collapsing soln MUST lead to AT LEAST the matches specified in
    the L-matrix .... and if MORE matches can be achieved, then they
    could be DEMANDED by an L-matrix deeper in the hierarchy. So it
    does not matter if (for a given L-matrix) we don't get the largest
    number of matches possible -- we only need to show that the matches
    which it REQUIRES can be achieved.

    Issue (2) could be stochastically addressed (i.e. to any level of
    desired precision) by trying multiple bat matrices simultaneously
    and demanding that all aways agree  -- as this would place bounds
    on the possion probability for any one bad-bat choice to have been
    conincidentally bad.  E.g. if always N independent bat matrices
    were used, for T tests, with any bat matrix having a per-test
    probability p<<1 of reporting a "false" collapse, then the
    probability of getting no disagreements over T tests would be:

        p(no disareemants in T tests) = p(no disagreements in one test)^T
            = ( 1 - P(all succeeed or all pass in one test) )^T
            = ( 1 - p^N - (1-p)^N )^T
            = ( 1 - p^N - (1 - Np + N(N-1)/2 p^2 + ... ) )^T
            = ( Np - N(N-1)/2 p^2 + ... + last_term_in_binom_exp - p^N )^T
            = (Np)^T ( 1 - (N-1)/2 p + (N-1)(N-2)/6 p^2 - ... + (1/(Np)) last_term_in_binom_exp - (1/N) p^(N-1) )^T
            = O( (Np)^T )

    ... hmm ... not a nice expansion ... have to think more about this.
    ... but something like this.
    
    """

    def __init__(self, M, k,
                 debug=False,
                 seed=0, starting_sigma=1,
                 unscaled_bad_bat_matrix = None,
                 ):
        self.M = M
        self.k = k
        self.debug = debug

        if unscaled_bad_bat_matrix is not None:
            self.unscaled_bad_bat_matrix = unscaled_bad_bat_matrix
            # Check that it is consistent with supplied M and k
            if unscaled_bad_bat_matrix.shape != (M, k):
                raise ValueError(f"In Rational_Decider constructor, bat matrix {unscaled_bad_bat_matrix} does not have expected shape {(M, k)}.")

            # The shapes were good, so we are done constructing:
            return

        # User did not supply a bat matrix, so it is now our job to make one given the supplied parameters.
        assert unscaled_bad_bat_matrix == None

        # Prepare bad bats
        self.unscaled_bad_bat_matrix = spt.general_position_integer_bat_matrix(M=M, k=k, seed=seed, starting_sigma=starting_sigma)

    def __repr__(self):
        return f"Rational_Decider(M={self.M}, k={self.k}, bat_matrix={self.unscaled_bad_bat_matrix})"

    def confusable_sets_or_None(self, L_matrix: sp.Matrix):
        return confusable_sets_or_None(L_matrix, self.unscaled_bad_bat_matrix, self.M)

    def vote_for_collapse_and_null_space(self, L_matrix : sp.Matrix) -> tuple:
        return vote_for_collapse_and_null_space(L_matrix, self.unscaled_bad_bat_matrix, self.M)

    def function_factory(self):
        return lambda mat : self.confusable_sets_or_None(mat)

#########################################################



def confusable_sets_or_None(L_matrix : sp.Matrix, unscaled_bad_bat_matrix:sp.Matrix, M:int):
    """
    If no scaling of the bad-bat lattice achieving the matches in L_matrix generates confusable sets,
    return None.
    Else, return a tuple containing a pair of confusable sets obtained from the bad-bat lattice (using the L_matrix matches)
    as well as the scaled_bad_bat_matrix that makes them.
    """

    vote_for_collapse, null_space = vote_for_collapse_and_null_space(L_matrix, unscaled_bad_bat_matrix, M)

    if vote_for_collapse:
        assert 0 <= len(null_space) <= M # Note <= not < in first inequality.
        return None

    assert not vote_for_collapse
    assert 0 < len(null_space) <= M # Note < not <= in first inequality.

    # OK - it is now our job to generate some confusable sets!
    import nonzero_lin_comb

    null_space_contribs, point_in_null_space =  nonzero_lin_comb.combine_many(null_space)

    assert len(null_space_contribs) == len(null_space)
    assert point_in_null_space.shape == (M, 1)
    assert all(alpha!=0 for alpha in point_in_null_space)

    import confusable_multisets

    scaled_bad_bat_matrix = confusable_multisets.scaled_bad_bat_matrix(unscaled_bad_bat_matrix, point_in_null_space)

    ##  E, O, C, EE, OO = confusable_multisets.analyze_B(scaled_bad_bat_matrix, plot_if_2d=False, show_C_if_plotting = False)
    _, _, _, EE, OO = confusable_multisets.mitm_compute_E_O_C_EE_OO(scaled_bad_bat_matrix)

    assert EE.total() == OO.total(), f"Must have {EE.total()}=={OO.total()} when scaled_bad_bat_matrix = {scaled_bad_bat_matrix}"

    return (EE, OO, scaled_bad_bat_matrix)

def vote_for_collapse_and_null_space(L_matrix: sp.Matrix, bat_matrix: sp.Matrix, M: int) -> tuple:
        """
        A key returned element is the null-space. It is probably only interesting when the
        null space is present.  The null space is a basis
        for the space of set of scales (alpha_1, alpha_2, ..., alpha_M) which, if applied
        to the unscaled bad-bat directions, would lead the bad-bat lattice to exhibit the
        matches specified in L_matrix.  A general point in the null space will always
        result in the L_matrix matches working, but might (sometimes) also lead to 
        bad bat lattice collapse, e.g. as would happen if at least one of the alphas was
        zero.  However, "vote for collapse" being False will assure us that there is 
        at least one point in the null space where non-collapse happens.

        In principle M could be found as the number of columns of L_matrix, or as the number
        of rows of bat_matrix, but this can fail for a few special cases, such as when
        L_matrix has no rows at all, etc. So to be on the safe side we pass in M
        """
        if L_matrix.rows != 0 and L_matrix.cols != M:
            raise ValueError()

        big_mat = alpha_attacking_matrix(L_matrix, bat_matrix)

        null_space = big_mat.nullspace()

        """
        Example output of

        print("---------------------------------")
        print(f"L_matrix   = {L_matrix}")
        print(f"big_mat    = {big_mat}")
        print(f"null space = {null_space}")
        print("---------------------------------")

        ---------------------------------
        L_matrix   = Matrix([[-1, -1, -1, -1, -1, -1, -1, 0], [-1, -1, -1, -1, 0, 0, 0, -1]])
        big_mat    = Matrix([[-13, -64, 54, -130, 70, 62, 233, 0], [13, -10, -36, -95, 127, -4, 22, 0], [-13, -64, 54, -130, 0, 0, 0, 125], [13, -10, -36, -95, 0, 0, 0, 73]])
        null space = [Matrix([
        [1422/481],
        [    9/37],
        [       1],
        [       0],
        [       0],
        [       0],
        [       0],
        [       0]]), Matrix([
        [2390/481],
        [ -225/74],
        [       0],
        [       1],
        [       0],
        [       0],
        [       0],
        [       0]]), Matrix([
        [          0],
        [          0],
        [          0],
        [          0],
        [ -1148/4077],
        [-28051/8154],
        [          1],
        [          0]]), Matrix([
        [ -1711/481],
        [     99/37],
        [         0],
        [         0],
        [ 2513/4077],
        [10765/8154],
        [         0],
        [         1]])]
        ---------------------------------
        """

        null_space_dimension = len(null_space)

        if null_space_dimension == 0:
            # The null-space is 0-dimensional, i.e. the only soln has all alphas equal to zero.
            # This is the strongest form of collapse we can have!
            return True, null_space

        if null_space_dimension == M:
            # The null-space is M-dimensional. There are only M alphas, so it is capable of spanning to
            # (alpha1, alpha2, .... ) = (1, 1, ... ).
            # So no collapse here!
            return False, null_space

        # OK - the shortcuts didn't work, so fall back to the default method that will always work.
        # Namely, for the reasons set out in "OneNote -> Research -> Symmetries -> Non collapsing null space"
        # it may be proved that 
        #
        #    (there exists a solution with all alphai != 0) iff (each coordinate direction is non-zero in at least one basis vector)
        #
        # or equivalently:
        #
        #    (there exists no solution with all alphai != 0) iff (at least one coordinate direction is zero in all basis vectors)
        #
        # so test if there is a non-zero value in each component:

        for cpt_index in range(M):
            if all( (basis_vec[cpt_index, 0] == 0) for basis_vec in null_space  ):
                # There are zeros in every basis vector at this cpt,
                # so there is no soln with every alphai nonzero, 
                # so this matrix collapses the bad bat lattice, and
                # so no coonfusable sets were found here
                return True, null_space
            # Try next cpt

        # OK - if we got here we tried all components but didn't find one that had zero in every basis vector,
        # so this matrix does admit solns with all alphai non-zero, i.e. the matrix does not collapse.
        # See proof in "OneNote -> Research -> Symmetries -> Non collapsing null space".
        return False, null_space


def alpha_attacking_matrix(
                L_matrix : sp.Matrix,  # Any number of rows, but M columns.
                bat_matrix : sp.Matrix, # M rows, and k columns, so that each row is a bat vector.
                ) -> sp.Matrix :
    """
    This method generates the matrix A for which the solns of A.(vec of alphas) are the same as the solutions to L.(alpha1 w1, alpha2 w2, ... , alphaM, wM) where w1 is the first bat (i.e. first row of bat matrix) and w2 the second, and so on.
    """

    M, k = bat_matrix.shape
    R, M_L = L_matrix.shape

    assert M == M_L or R==0, f"L and B must work with the same no. of vectors. Wanted M({M}) == M_L({M_L})"

    effective_row = 0
    ans = sp.zeros(R*k, M)
    for i in range(R):
        for kk in range(k):
            for j in range(M):
                ans[effective_row, j] = L_matrix[i,j]*bat_matrix[j,kk]
            effective_row += 1

    return ans


