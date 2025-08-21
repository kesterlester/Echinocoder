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

    def __collapse_test_case(self, L_matrix : sp.Matrix) -> str:
        return f"""

import deciders
from sympy import Matrix
L_matrix = {repr(L_matrix)}
bat_matrices = {self.bat_matrices}
decider = deciders.Rational_Decider(M={self.M}, k={self.k}, bat_matrices=bat_matrices)
votes_for_collapse = decider.votes_for_collapse(L_matrix)
print("votes_for_collapse are: ",votes_for_collapse)
has_True = True in votes_for_collapse
has_False = True in votes_for_collapse
assert not ( has_True and has_False )


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
        print("Consider this unit test case")
        print(self.__collapse_test_case(L_matrix))
        #print("Bat matrix list was:")
        #print(f"{self.bat_matrices}")
        #for i, bm in enumerate(self.bat_matrices):
        #    print(f"  ({i}): {repr(bm)}")

       
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
        return False

    def function_factory(self):
        return lambda mat : self.matrix_does_not_collapse(mat)

#########################################################
