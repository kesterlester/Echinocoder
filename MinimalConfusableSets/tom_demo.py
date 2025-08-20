import sympy as sp
from functools import partial
from vertex_matches import generate_viable_vertex_match_matrices
from sympy_tools import strip_zero_rows, some_row_causes_collapse
import decider_functions.decider_function as df
from Match_Tracker import Match_Tracker


#########################################################

### def old_demo():
###     print("== Test of Matrix Generation =========")
### 
###     def max_row_requirement(mat, max_rows):
###         return sp.shape(mat)[0] <= max_rows
### 
###     def f(mat):
###         return sp.shape(mat)[0] <= 4 # True if mat has 4 or fewer rows.
### 
###     mat_gen = generate_viable_vertex_match_matrices(
###         M=M,
###         k=k,
###         # All of the next three lines have the same effect, but different pros/cons.
###         # Try changing which one(s) is(are) commented out.
###         #yield_matrix = partial(max_row_requirement, max_rows=4),
###         #go_deeper = partial(max_row_requirement, max_rows=3), # fastest option, where possible
###         yield_matrix = f,
###         )
### 
###     return(mat_gen)


class Decider:
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

#########################################################



def demo(M_and_k_tuple=None):

    def max_row_requirement(mat, max_rows):
        return sp.shape(mat)[0] <= max_rows

    def f(mat: sp.Matrix):
        return sp.shape(mat)[0] <= 5 # True if mat has 4 or fewer rows.



    for M,k in (() if M_and_k_tuple is None else (M_and_k_tuple,)) + (
            #(5,4),
            #(7,4),
            #(7,3),
            ##(7,2),
            #(9,4),
            #(11,4),
            #(13,4),
            #(15,4),
        ):
        print()
        print()
        size_of_smallest_confusable_set_constructed_so_far = 2**(M-1)
        smallest_set_mat = sp.Matrix()
        print( "====================================================================")
        print(f"For M={M} and k={k} the initial E and O sets have size {size_of_smallest_confusable_set_constructed_so_far}.")
        print( "--------------------------------------------------------------------")
        print(f"For M={M} and k={k} the not obviously bad vertex match matrices are:")
        print( "--------------------------------------------------------------------")

        debug = False

        collapse_checker = Decider(M=M, k=k, debug=debug)
        collapse_checking_function = collapse_checker.function_factory()

        mat_gen = generate_viable_vertex_match_matrices(
            M=M,
            k=k,
            return_mat = True,
            return_hashable_rre = True,
            remove_duplicates_via_hash = True,
            yield_matrix = collapse_checking_function,
            # yield_matrix = partial(max_row_requirement, max_rows=4),
            # go_deeper = partial(max_row_requirement, max_rows=3), # fastest option, where possible
            # yield_matrix = partial(matrix_is_not_definitely_bad, k=k),
            debug = debug,
            debug_test_max_rows=True,
            )

        mat_gen_slow = generate_viable_vertex_match_matrices(
            M=M,
            k=k,
            return_mat = True,
            return_hashable_rre = True,
            remove_duplicates_via_hash = True,
            yield_matrix = collapse_checking_function,
            # yield_matrix = partial(max_row_requirement, max_rows=4),
            # go_deeper = partial(max_row_requirement, max_rows=3), # fastest option, where possible
            # yield_matrix = partial(matrix_is_not_definitely_bad, k=k),
            debug = debug,
            debug_test_max_rows=False,
            )
       
        number_enumerated = 0
        ((last_mat,last_rre), (last_mat_slow,last_rre_slow)) = ((None,None),(None,None))
        #for i, (mat,rre) in enumerate(mat_gen):
        for i, ((mat,rre), (mat_slow,rre_slow)) in enumerate(zip(mat_gen, mat_gen_slow)):
            #tracker = Match_Tracker(M, mat)
            #e_vertices = tracker.number_of_even_vertices_present()
            #print(f"    {i} ev:{e_vertices},  raw={mat}, rre={repr(rre)}")
            pr = False
            mismatch = ((mat != mat_slow) or (rre != rre_slow))
            if mismatch:
                print(f"(Pre){i-1} FAST raw={mat     }, rre={repr(rre)     }")
                print(f"(pre){i-1} SLOW raw={mat_slow}, rre={repr(rre_slow)}")
                pr = True
            if i % 500 == 0 or pr:
                print(f"     {i  } FAST raw={mat     }, rre={repr(rre)     }")
                print(f"     {i  } SLOW raw={mat_slow}, rre={repr(rre_slow)}")
            ((last_mat,last_rre), (last_mat_slow,last_rre_slow)) = ((mat,rre), (mat_slow,rre_slow)) 
            if mismatch:
                 assert False, "Stopping after mismatch"
            ### if (
            ###    e_vertices < size_of_smallest_confusable_set_constructed_so_far or 
            ###    e_vertices <= size_of_smallest_confusable_set_constructed_so_far and (sp.shape(mat)[0])<(sp.shape(smallest_set_mat)[0])
            ###    ): # better, or as good but with fewer rows
            ###     size_of_smallest_confusable_set_constructed_so_far = e_vertices
            ###     smallest_set_mat = mat
            ###     print(f"SET ABOVE WAS SMALLEST CONFUSABLE SET SO FAR FOUND: ev:{e_vertices}, M={M}, k={k}")

            number_enumerated += 1

        print(f"There were {number_enumerated} found for M={M} and k={k}.")
        #print(f"The smallest confusable set was found for ev:{size_of_smallest_confusable_set_constructed_so_far}, M={M}, k={k}, smallest_set_mat = {smallest_set_mat}")
        print("====================================================================")

if __name__ == "__main__":
    import sys
    print(f"sys.argv = {sys.argv}")
    if len(sys.argv) == 3:
        M, k = (int(n) for n in sys.argv[1:3])
        demo( (M,k) )
    else:
        demo()
