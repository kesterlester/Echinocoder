import sympy as sp
from functools import partial
from vertex_matches import generate_viable_vertex_match_matrices
from decider_function import prepare_B, decide_collapse

k = 2
M = 5

def demo():
    print("== Test of Matrix Generation =========")

    def max_row_requirement(mat, max_rows):
        return sp.shape(mat)[0] <= max_rows

    def f(mat):
        return sp.shape(mat)[0] <= 4 # True if mat has 4 or fewer rows.

    mat_gen = generate_viable_vertex_match_matrices(
        M=M,
        k=k,
        # All of the next three lines have the same effect, but different pros/cons.
        # Try changing which one(s) is(are) commented out.
        #yield_matrix = partial(max_row_requirement, max_rows=4),
        #go_deeper = partial(max_row_requirement, max_rows=3), # fastest option, where possible
        yield_matrix = f,
        ) 

    return(mat_gen)
# prepare bad bats once. Arguments explained in decider_function
B = prepare_B(k=k, M=M, method="electrostatic", iters = 500, learning_rate = 0.01, power = 2.0, sample = 'rr', spread = 'gauss')

# loop over all L's without regenerating bad bats.
# ignore num_trials on trust that nullspace basis is dense with totally non-zero vectors
# tol is tolerance when checking if eg alpha = 0

for L in demo():
    tom_answer = decide_collapse(L, B, num_trials=1000, tol=1e-12)
    print(tom_answer)

