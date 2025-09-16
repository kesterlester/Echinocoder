#!/usr/bin/env python3

import sympy as sp
from functools import partial
from vertex_matches import generate_viable_vertex_match_matrices
import sympy_tools as spt
import decider_functions.decider_function as df
from Match_Tracker import Match_Tracker

import deciders
import config
import sys

def demo(M_and_k_tuple=None, show_cs=False, 
              scale_cols=False, # Warning - breaks meaning of bats
              use_hash=False,  # Warning - may not be rigorous
              prune_max_rows = False, # Warning - buggy
              prune_short_rows = True, # safe
              prune_pivots = False, # Warning - buggy
              print_period = 10000,
              ):

    cfg = f"config: hash={use_hash}, prune_short_rows={prune_short_rows}, prune_max_rows={prune_max_rows}, prune_pivots={prune_pivots}, odd_ones={config.only_output_odd_ones}"

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
        collapse_checker_1 = deciders.Rational_Decider(M=M, k=k, debug=debug, seed=10, starting_sigma=10) # 1
        print("Creating slow decider")
        collapse_checker_2 = deciders.Rational_Decider(M=M, k=k, debug=debug, seed=0, starting_sigma=100) # 100
        print("Creating vslow decider")
        collapse_checker_3 = deciders.Rational_Decider(M=M, k=k, debug=debug, seed=100, starting_sigma=10000) # 10000
        print("done")
        collapse_checking_function_1 = collapse_checker_1.function_factory()
        collapse_checking_function_2 = collapse_checker_2.function_factory()
        collapse_checking_function_3 = collapse_checker_3.function_factory()

        mat_gen_1= generate_viable_vertex_match_matrices(
            M=M,
            k=k,
            sort_cols = False,
            return_mat = True,
            return_hashable_rre = True,
            return_confusable_sets = True,
            remove_duplicates_via_hash = use_hash,
            prune_max_rows = prune_max_rows,
            prune_short_rows = prune_short_rows,
            prune_pivots = prune_pivots,
            confusable_sets_or_None_function = collapse_checking_function_1,
            # yield_matrix = partial(max_row_requirement, max_rows=4),
            # go_deeper = partial(max_row_requirement, max_rows=3), # fastest option, where possible
            # yield_matrix = partial(matrix_is_not_definitely_bad, k=k),
            debug = debug,
            )

        mat_gen_2 = generate_viable_vertex_match_matrices(
            M=M,
            k=k,
            sort_cols = False,
            return_mat = True,
            return_hashable_rre = True,
            return_confusable_sets = True,
            remove_duplicates_via_hash = use_hash,
            prune_max_rows = prune_max_rows,
            prune_short_rows = prune_short_rows,
            prune_pivots = prune_pivots,
            confusable_sets_or_None_function = collapse_checking_function_2,
            # yield_matrix = partial(max_row_requirement, max_rows=4),
            # go_deeper = partial(max_row_requirement, max_rows=3), # fastest option, where possible
            # yield_matrix = partial(matrix_is_not_definitely_bad, k=k),
            debug = debug,
            )
      
        mat_gen_3 = generate_viable_vertex_match_matrices(
            M=M,
            k=k,
            sort_cols = False,
            return_mat = True,
            return_hashable_rre = True,
            return_confusable_sets = True,
            remove_duplicates_via_hash = use_hash,
            prune_max_rows = prune_max_rows,
            prune_short_rows = prune_short_rows,
            prune_pivots = prune_pivots,
            confusable_sets_or_None_function = collapse_checking_function_3,
            # yield_matrix = partial(max_row_requirement, max_rows=4),
            # go_deeper = partial(max_row_requirement, max_rows=3), # fastest option, where possible
            # yield_matrix = partial(matrix_is_not_definitely_bad, k=k),
            debug = debug,
            )
      
        for name, mat_gen, decider in (
                (f"VSLW", mat_gen_3, collapse_checker_3),
                #(f"SLOW", mat_gen_2, collapse_checker_2),
                #(f"FAST", mat_gen_1, collapse_checker_1),
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

                if i % print_period == 0:
                    print(f"{name}: CURRENT {i}, M={M}, k={k}, raw={mat}, rre={repr(rre)}, EE.total()={EE.total()}, OO.total()={OO.total()}\n")
                    mes = f"\n\nbest_siz={smallest_siz_so_far} for M={M}, k={k},\n" \
                          f"{cfg}\n" \
                          f"CMD: {sys.argv}\n"
                    print(mes)

                if new_best:
                    prefix = f"{name} SO FAR: "
                    mes = f"\n\nbest_siz={smallest_siz_so_far} for M={M}, k={k},\n" \
                          f"{cfg}\n" \
                          f"CMD: {sys.argv}\n" \
                          f"raw=\n{repr(best_mat)},\nrre=\n{repr(best_rre)},\n"\
                          f"unscaled_bad_bats=\n{repr(decider.unscaled_bad_bat_matrix)},\n"\
                          f"scalings={sp.srepr(best_scalings)}\n" \
                          f"scalingsSREPR={sp.srepr(best_scalings)}\n" \
                          f"scaled_bad_bats=\n{sp.srepr(best_scaled_bad_bats)}.\n"\
                          f"{number_enumerated} matrices have been scanned.\n\n"
                    for line in mes.split("\n"):
                        print(prefix, line)

                number_enumerated += 1

            print("---------------------- ENDING -------------------")
            prefix = f"{name} AT END: "
            mes = f"\n\nbest_siz={smallest_siz_so_far} for M={M}, k={k},\n" \
                  f"{cfg}\n" \
                  f"CMD: {sys.argv}\n" \
                  f"raw=\n{repr(best_mat)},\nrre=\n{repr(best_rre)},\n" \
                  f"unscaled_bad_bats=\n{repr(decider.unscaled_bad_bat_matrix)},\n" \
                  f"scalings={repr(best_scalings)}\n" \
                  f"scalingsSREPR={sp.srepr(best_scalings)}\n" \
                  f"scaled_bad_bats=\n{sp.srepr(best_scaled_bad_bats)}.\n" \
                  f"{number_enumerated} matrices have been scanned.\n\n"
            for line in mes.split("\n"):
               print(prefix, line)
            print("---------------------- END -------------------")
                
            import confusable_multisets as cs
            #cs.analyze_B(best_scaled_bad_bats, plot_if_2d=True, show_C_if_plotting=True)
            #cs.analyze_B(best_scaled_bad_bats, plot_if_2d=True)
            print("====================================================================")

def str_to_bool(value):
    if isinstance(value, bool):
        return value
    if value.lower() in {'false', 'f', '0', 'no', 'n'}:
        return False
    elif value.lower() in {'true', 't', '1', 'yes', 'y'}:
        return True
    raise ValueError(f'{value} is not a valid boolean value')

def ddd():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("M", type=int, help="number of bad bats")
    parser.add_argument("k", type=int, help="dimension of space")
    parser.add_argument("-f", "--frequency", type=int, help="period of status updates.", default=10000)
    parser.add_argument("-s", "--prune-short-rows",  
                        type=str_to_bool, default=True, nargs='?', const=True,
                        help="if true, discount L-matrix if any row has fewer than k+1 non-zero elements")
    parser.add_argument("-o", "--only-odd-ones",
                        type=str_to_bool, default=False, nargs='?', const=True,
                        help="if true, require all vertex match deltas to have an odd number of ones")
    parser.add_argument("-H", "--prune-hash",
                        type=str_to_bool, default=False, nargs='?', const=True,
                        help="if true, prune tree if RREF of L-matrix was seen before")
    parser.add_argument("-m", "--prune-max-rows",
                        type=str_to_bool, default=False, nargs='?', const=True,
                        help="if true, apply max row constraint to RREF of L-matrix")
    parser.add_argument("-p", "--prune-pivots",
                        type=str_to_bool, default=False, nargs='?', const=True,
                        help="if true, apply constraint to RREF L_matrix pivot position")
    args = parser.parse_args()

    config.only_output_odd_ones=args.only_odd_ones
    demo( (args.M,args.k), 
            use_hash = args.prune_hash,
            prune_pivots = args.prune_pivots,
            prune_max_rows = args.prune_max_rows,
            prune_short_rows = args.prune_short_rows,
            print_period = args.frequency,
            )

if __name__ == "__main__":
    ddd()
