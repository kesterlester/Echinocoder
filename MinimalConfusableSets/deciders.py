import sympy as sp
from functools import partial
from vertex_matches import generate_viable_vertex_match_matrices
import sympy_tools as spt
import decider_functions.decider_function as df
from Match_Tracker import Match_Tracker
from collections import namedtuple
import confusable_multisets

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
        return confusable_multisets.confusable_sets_or_None(L_matrix, self.unscaled_bad_bat_matrix, self.M)

    def vote_for_collapse_and_null_space(self, L_matrix : sp.Matrix) -> tuple:
        return confusable_multisets.vote_for_collapse_and_null_space(L_matrix, self.unscaled_bad_bat_matrix, self.M)

    def function_factory(self):
        return lambda mat : self.confusable_sets_or_None(mat)

#########################################################


