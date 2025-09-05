import sympy as sp
from functools import partial
from vertex_matches import generate_viable_vertex_match_matrices, alpha_attacking_matrix
import sympy_tools as spt
import decider_functions.decider_function as df
from Match_Tracker import Match_Tracker

import deciders

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

        collapse_checker_1 = deciders.Rational_Decider(M=M, k=k, debug=debug, seed=10, starting_sigma=1)
        collapse_checker_2 = deciders.Rational_Decider(M=M, k=k, debug=debug, seed=0, starting_sigma=100)
        collapse_checker_3 = deciders.Rational_Decider(M=M, k=k, debug=debug, seed=100, starting_sigma=10000)
        collapse_checking_function_1 = collapse_checker_1.function_factory()
        collapse_checking_function_2 = collapse_checker_2.function_factory()
        collapse_checking_function_3 = collapse_checker_3.function_factory()

        mat_gen_fast = generate_viable_vertex_match_matrices(
            M=M,
            k=k,
            return_mat = True,
            return_hashable_rre = True,
            return_confusable_sets = True,
            remove_duplicates_via_hash = True,
            confusable_sets_or_None_function = collapse_checking_function_1,
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
            return_confusable_sets = True,
            remove_duplicates_via_hash = True,
            confusable_sets_or_None_function = collapse_checking_function_2,
            # yield_matrix = partial(max_row_requirement, max_rows=4),
            # go_deeper = partial(max_row_requirement, max_rows=3), # fastest option, where possible
            # yield_matrix = partial(matrix_is_not_definitely_bad, k=k),
            debug = debug,
            debug_test_max_rows=True,
            )
      
        mat_gen_very_slow = generate_viable_vertex_match_matrices(
            M=M,
            k=k,
            return_mat = True,
            return_hashable_rre = True,
            return_confusable_sets = True,
            remove_duplicates_via_hash = True,
            confusable_sets_or_None_function = collapse_checking_function_3,
            # yield_matrix = partial(max_row_requirement, max_rows=4),
            # go_deeper = partial(max_row_requirement, max_rows=3), # fastest option, where possible
            # yield_matrix = partial(matrix_is_not_definitely_bad, k=k),
            debug = debug,
            debug_test_max_rows=True,
            )
      

        for name, mat_gen in (
                ("VSLW", mat_gen_very_slow),
                ("SLOW", mat_gen_slow),
                ("FAST",mat_gen_fast),
                ):

            print("")
            print("=================================")
            print(f"STARTING {name}")
            print("=================================")

            number_enumerated = 0
            smallest_siz_so_far = None

            for i, (mat,rre,(EE,OO)) in enumerate(mat_gen):

                pr = False
                siz = EE.total()

                if smallest_siz_so_far == None or siz < smallest_siz_so_far:
                    smallest_siz_so_far, smallest_EE, smallest_OO = siz, EE, OO
                    pr = True

                if i % 10000 == 0 or pr:
                    print(f"{name}:     {i}:  raw={mat}, rre={repr(rre)}, EE.total()={EE.total()}, OO.total()={OO.total()}     ")
                    print(f"{name}: The smallest confusable sets so far have {smallest_siz_so_far}=={smallest_EE.total()} points and are {smallest_EE} and {smallest_OO}.")
                    print()

                number_enumerated += 1

            print(f"{name}:  M={M} and k={k} smallest confusable set was size {smallest_siz_so_far} and was found after checking {number_enumerated} match matrices.")
            print("====================================================================")

def ddd():
    import sys
    print(f"sys.argv = {sys.argv}")
    if len(sys.argv) == 3:
        M, k = (int(n) for n in sys.argv[1:3])
        demo( (M,k) )
    else:
        demo()

if __name__ == "__main__":
    ddd()
