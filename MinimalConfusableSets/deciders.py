import sympy as sp
from functools import partial
from vertex_matches import generate_viable_vertex_match_matrices, alpha_attacking_matrix
import sympy_tools as spt
import decider_functions.decider_function as df
from Match_Tracker import Match_Tracker

class Old_Decider:
    def __init__(self, M, k, debug=False, method="electrostatic", iters = 500, learning_rate = 0.01, power = 2.0, sample = 'rr', spread = 'gauss'):
        self.M = M
        self.k = k
        self.debug = debug

        # prepare bad bats once. Arguments explained in decider_function
        self.B = df.prepare_B(k=k, M=M, method=method, iters = iters, learning_rate = learning_rate, power = power, sample = sample, spread = spread)

    def matrix_does_not_collapse(self, L_matrix : sp.Matrix):
         # The decide_collapse function below returns False on Collapse, True otherwise.
         ans = df.decide_collapse(L_matrix, self.B, num_trials=1000, tol=1e-12)
         if self.debug:
             if not ans:
                print(f"VETO by TOM (float) {L_matrix}")
         return ans

    def function_factory(self):
        return lambda mat : self.matrix_does_not_collapse(mat)


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
    know if this is true, though.

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
                 seed=0, starting_sigma=1, voting_copies=4,
                 bat_matrices = None,
                 ):
        self.M = M
        self.k = k
        self.debug = debug

        if bat_matrices is not None:
            self.bat_matrices = bat_matrices
            # Check that they are consistent with supplied M and k
            if any((bat_matrix.shape != (M,k)) for bat_matrix in bat_matrices):
                raise ValueError(f"In Rational_Decider constructor, bat matrix {bat_matrix} in {bat_matrices} does not have expected shape {(M,k)}.")

            # All the shapes were good, so we are done constructing:
            return

        # User did not supply bat matrices, so it is now our job to make them given the supplied parameters.
        assert bat_matrices == None
        assert voting_copies >= 1

        # Prepare bad bats
        self.bat_matrices = []
        for i in range(voting_copies):
            self.bat_matrices.append(spt.general_position_integer_bat_matrix(M=M, k=k, seed=seed*1234567*i, starting_sigma=starting_sigma))

    def __repr__(self):
        return f"Rational_Decider(M={self.M}, k={self.k}, bat_matrices={self.bat_matrices})"

    def __collapse_test_case(self, L_matrix : sp.Matrix, votes_for_collapse : tuple) -> str:
        return f"""
    
    
    ##################################
    import deciders
    from sympy import Matrix
    
    L_matrix = {repr(L_matrix)}
    bat_matrices = {self.bat_matrices}
    decider = deciders.Rational_Decider(M={self.M}, k={self.k}, bat_matrices=bat_matrices)
    
    votes_for_collapse = decider.votes_for_collapse(L_matrix)
    print("votes_for_collapse are: ",votes_for_collapse)
    
    # Expect votes_for_collapse={votes_for_collapse} 
    
    has_True = True in votes_for_collapse
    has_False = True in votes_for_collapse
    
    assert not ( has_True and has_False )
    ##################################
    
    
"""
        
    def votes_for_collapse(self, L_matrix : sp.Matrix) -> tuple:
        return tuple(self.__matrix_collapses(L_matrix, bat_matrix) for bat_matrix in self.bat_matrices)
        
    def matrix_does_not_collapse(self, L_matrix : sp.Matrix) -> bool:
        return not self.matrix_collapses(L_matrix)

    def matrix_collapses(self, L_matrix : sp.Matrix) -> bool:
        votes_for_collapse = self.votes_for_collapse(L_matrix)

        if not False in votes_for_collapse:
            # All votes are True so all agree on collapse:
            print(".", end="")
            return True

        if not True in votes_for_collapse:
            # All votes are False, so all agree on non collapse:
            print(".", end="")
            return False

        # Oh dear, some votes are True and some are False. This should not happen if
        # the principle on which method is based is fully understood.
        print(f"\n\nFor L_matrix = {L_matrix},")
        print(f"votes_for_collapse are {votes_for_collapse}.")
        print(f"Decider can be reconstructed like this:\n{repr(self)}\n")
        print("Consider this unit test case:")
        print(self.__collapse_test_case(L_matrix, votes_for_collapse))
       
        raise RuntimeError()
        assert False, "Implementation of matrix_collapse function seems to be broken."

    def __matrix_does_not_collapse(self, L_matrix : sp.Matrix, bat_matrix : sp.Matrix) -> bool:
        return not self.matrix_collapses(L_matrix, bat_matrix)

    def __matrix_collapses(self, L_matrix : sp.Matrix, bat_matrix : sp.Matrix) -> bool:

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
            return True

        if null_space_dimension == self.M:
            # The null-space is M-dimensional. There are only M alphas, so it is capable of spanning to 
            # (alpha1, alpha2, .... ) = (1, 1, ... ).
            # So no collapse here!
            return False

        # OK - the shortcuts didn't work, so fall back to the default method that would always work.
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

        for cpt_index in range(self.M):
            if all( (basis_vec[cpt_index, 0] == 0) for basis_vec in null_space  ):
                # There are zeros in every basis vector at this cpt,
                # so there is no soln with every alphai nonzero, 
                # so this matrix collapses:
                return True
            # Try next cpt
        # OK - we tried all components but didn't find one that had zero in every basis vector,
        # so this matrix does admit solns with all alphai non-zero, i.e. the matrix does not collapse.
        # See proof in "OneNote -> Research -> Symmetries -> Non collapsing null space".
        return False

    def function_factory(self):
        return lambda mat : self.matrix_does_not_collapse(mat)

#########################################################
