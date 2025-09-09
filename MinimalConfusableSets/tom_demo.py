import sympy as sp
from functools import partial
from vertex_matches import generate_viable_vertex_match_matrices
import sympy_tools as spt
import decider_functions.decider_function as df
from Match_Tracker import Match_Tracker

import deciders

def demo(M_and_k_tuple=None, show_cs=False, scale_cols=False):

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

        debug = False

        print("Creating fast decider")
        collapse_checker_1 = deciders.Rational_Decider(M=M, k=k, debug=debug, seed=10, starting_sigma=1)
        print("Creating slow decider")
        collapse_checker_2 = deciders.Rational_Decider(M=M, k=k, debug=debug, seed=0, starting_sigma=100)
        print("Creating vslow decider")
        collapse_checker_3 = deciders.Rational_Decider(M=M, k=k, debug=debug, seed=100, starting_sigma=10000)
        print("done")
        collapse_checking_function_1 = collapse_checker_1.function_factory()
        collapse_checking_function_2 = collapse_checker_2.function_factory()
        collapse_checking_function_3 = collapse_checker_3.function_factory()

        for sort_cols in (
               #True, 
               False,
               ):

            mat_gen_1= generate_viable_vertex_match_matrices(
                M=M,
                k=k,
                sort_cols = sort_cols,
                return_mat = True,
                return_rre = True,
                return_confusable_sets = True,
                remove_duplicates_via_hash = True,
                confusable_sets_or_None_function = collapse_checking_function_1,
                # yield_matrix = partial(max_row_requirement, max_rows=4),
                # go_deeper = partial(max_row_requirement, max_rows=3), # fastest option, where possible
                # yield_matrix = partial(matrix_is_not_definitely_bad, k=k),
                debug = debug,
                debug_test_max_rows=True,
                )

            mat_gen_2 = generate_viable_vertex_match_matrices(
                M=M,
                k=k,
                sort_cols = sort_cols,
                return_mat = True,
                return_rre = True,
                return_confusable_sets = True,
                remove_duplicates_via_hash = True,
                confusable_sets_or_None_function = collapse_checking_function_2,
                # yield_matrix = partial(max_row_requirement, max_rows=4),
                # go_deeper = partial(max_row_requirement, max_rows=3), # fastest option, where possible
                # yield_matrix = partial(matrix_is_not_definitely_bad, k=k),
                debug = debug,
                debug_test_max_rows=True,
                )
      
            mat_gen_3 = generate_viable_vertex_match_matrices(
                M=M,
                k=k,
                sort_cols = sort_cols,
                return_mat = True,
                return_rre = True,
                return_confusable_sets = True,
                remove_duplicates_via_hash = True,
                confusable_sets_or_None_function = collapse_checking_function_3,
                # yield_matrix = partial(max_row_requirement, max_rows=4),
                # go_deeper = partial(max_row_requirement, max_rows=3), # fastest option, where possible
                # yield_matrix = partial(matrix_is_not_definitely_bad, k=k),
                debug = debug,
                debug_test_max_rows=True,
                scale_cols_in_hash=scale_cols,
                )
      

            for name, mat_gen, decider in (
                #    (f"SLOW sort={sort_cols}", mat_gen_2, collapse_checker_2),
                    (f"VSLW sort={sort_cols} scale_cols={scale_cols}", mat_gen_3, collapse_checker_3),
                #    (f"FAST sort={sort_cols}", mat_gen_1, collapse_checker_1),
                    ):

                print("")
                print("=================================")
                print(f"STARTING {name}")
                print("=================================")

                number_enumerated = 0
                smallest_siz_so_far = None

                for i, (mat,rre,(EE,OO,scalings, scaled_bad_bats)) in enumerate(mat_gen):

                    new_best = False
                    siz = EE.total()

                    if smallest_siz_so_far == None or siz < smallest_siz_so_far:
                        smallest_siz_so_far = siz
                        # smallest_EE, smallest_OO = EE, OO
                        best_scalings, best_scaled_bad_bats = scalings, scaled_bad_bats
                        best_mat, best_rre = mat, rre
                        new_best = True

                    if i % 10000 == 0:
                        print(f"{name}: CURRENT {i}, M={M}, k={k}, raw={mat}, rre={repr(rre)}, EE.total()={EE.total()}, OO.total()={OO.total()}\n")

                    if new_best:
                        prefix = f"{name} SO FAR: "
                        mes = f"\n\nfor M={M}, k={k} the smallest confusable sets have size {smallest_siz_so_far},\n"\
                              f"raw=\n{repr(best_mat)},\nrre=\n{repr(best_rre)},\n"\
                              f"unscaled_bad_bats=\n{repr(decider.unscaled_bad_bat_matrix)},\n"\
                              f"scalings={sp.srepr(best_scalings)}\n" \
                              f"scalingsSREPR={sp.srepr(best_scalings)}\n" \
                              f"scaled_bad_bats=\n{sp.srepr(best_scaled_bad_bats)}.\n"\
                              f"{number_enumerated} matrices have been scanned.\n\n"
                        for line in mes.split("\n"):
                            print(prefix, line)

                    number_enumerated += 1

                prefix = f"{name} AT END: "
                mes = f"\n\nfor M={M}, k={k} the smallest confusable sets have size {smallest_siz_so_far},\n" \
                      f"raw=\n{repr(best_mat)},\nrre=\n{repr(best_rre)},\n" \
                      f"unscaled_bad_bats=\n{repr(decider.unscaled_bad_bat_matrix)},\n" \
                      f"scalings={repr(best_scalings)}\n" \
                      f"scalingsSREPR={sp.srepr(best_scalings)}\n" \
                      f"scaled_bad_bats=\n{sp.srepr(best_scaled_bad_bats)}.\n" \
                      f"{number_enumerated} matrices have been scanned.\n\n"
                for line in mes.split("\n"):
                   print(prefix, line)
                    
                import confusable_multisets as cs
                #cs.analyze_B(best_scaled_bad_bats, plot_if_2d=True, show_C_if_plotting=True)
                #cs.analyze_B(best_scaled_bad_bats, plot_if_2d=True)
                print("====================================================================")

def ddd():
    import sys
    print(f"sys.argv = {sys.argv}")
    if len(sys.argv) == 3:
        M, k = (int(n) for n in sys.argv[1:3])
        demo( (M,k) )
    elif len(sys.argv) == 4:
        M, k = (int(n) for n in sys.argv[1:3])
        demo( (M,k), scale_cols = True if sys.argv[3]=="1" else False  )
    else:
        demo()

if __name__ == "__main__":
    ddd()
