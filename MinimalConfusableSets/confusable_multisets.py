import sympy as sp

def scaled_bad_bat_matrix(B_unscaled, # B_unscaled is an UNSCALED bad bat matrix. This means that it is a matrix with m rows, each being a k-dim direction vector which is perp to lots of good bats. The edges of the "even-odd hypercube construction" from which the alternate red and blue points are extracted to make the confusable multisets are scaled versions of these rows.  The lengths of these vectors are yet to be scaled by the alphas which are calculated from the null space of an L-vertex-match-matrix which has been appropriatley transformed to act on the alphas -- hence "unscaled"
                          vector_of_alphas, # vector_of_alphas is an m-dimensional vector, whose first component tells us how much to scale the first row of B, and whose second component tells us how much to scale the second row of B, and so on
                          ):

    return scale_rows(B_unscaled, vector_of_alphas)

def confusable_multisets_given(B_scaled, # B_scaled is a SCALED bad bat matrix. This means that it is a matrix with m rows, each being a k-dim direction vector which is perp to lots of good bats. The edges of the "even-odd hypercube construction" from which the alternate red and blue points are extracted to make the confusable multisets are these rows
                              ):
    pass




def scale_rows(M: sp.Matrix, a) -> sp.Matrix:
    """
    Returns a matrix S that is like M except that
    row i of S is row i of M multiplied by a[i].

    Parameters:
    - M: sympy.Matrix (m x n)
    - a: sympy.Matrix, list, or tuple of length m

    Returns:
    - sympy.Matrix (m x n)
    """
    # Convert a to sympy column vector if it isn't already
    if isinstance(a, (list, tuple)):
        a = sp.Matrix(a)
    elif not isinstance(a, sp.Matrix):
        raise TypeError("Argument 'a' must be a sympy.Matrix, list, or tuple.")

    # Ensure column vector
    a = a.reshape(len(a), 1)

    # Input validation
    assert M.rows == a.rows, (
        f"Dimension mismatch: M has {M.rows} rows, "
        f"but vector 'a' has length {a.rows}."
    )

    # Construct diagonal matrix from 'a' and scale rows
    return sp.diag(*a) * M


from sympy import Matrix
from collections import Counter
import matplotlib.pyplot as plt

# --- Helper functions ---

def _row_tuple_rows(rows):
    """Convert sympy row objects to tuples of rationals/ints."""
    return [tuple(r) for r in rows]

def subset_sums_with_parity(rows):
    """Return (even_counter, odd_counter) for all subset sums of given rows."""
    n = len(rows)
    even = Counter()
    odd = Counter()
    for mask in range(1 << n):
        parity = (bin(mask).count("1") % 2)
        s = [0]*len(rows[0])
        for i in range(n):
            if (mask >> i) & 1:
                row = rows[i]
                for j in range(len(s)):
                    s[j] += row[j]
        s_t = tuple(s)
        if parity == 0:
            even[s_t] += 1
        else:
            odd[s_t] += 1
    return even, odd

def convolve_counters(A: Counter, B: Counter):
    """Convolution of counters of k-vectors."""
    out = Counter()
    for a, ca in A.items():
        for b, cb in B.items():
            s = tuple(x+y for x,y in zip(a,b))
            out[s] += ca * cb
    return out

def mitm_compute_E_O_C_EE_OO(B: Matrix):
    """Meet-in-the-middle computation of E, O, C, EE, OO for matrix B."""
    m, k = B.shape
    rows = _row_tuple_rows([B.row(i) for i in range(m)])
    mid = m // 2
    L_rows = rows[:mid]
    R_rows = rows[mid:]
    L_even, L_odd = subset_sums_with_parity(L_rows)
    R_even, R_odd = subset_sums_with_parity(R_rows)
    E1 = convolve_counters(L_even, R_even)
    E2 = convolve_counters(L_odd, R_odd)
    O1 = convolve_counters(L_even, R_odd)
    O2 = convolve_counters(L_odd, R_even)
    E = E1 + E2
    O = O1 + O2
    C = Counter()
    for key in set(E.keys()) & set(O.keys()):
        C[key] = min(E[key], O[key])
    EE = Counter()
    OO = Counter()
    for key, ce in E.items():
        cc = C.get(key, 0)
        if ce > cc:
            EE[key] = ce - cc
    for key, co in O.items():
        cc = C.get(key, 0)
        if co > cc:
            OO[key] = co - cc
    return E, O, C, EE, OO

def expand_points(counter):
    """Expand Counter of k-vectors into list of points with multiplicities."""
    pts = []
    for pt, c in counter.items():
        pts.extend([pt]*c)
    return pts

def analyze_B(B: Matrix, plot_if_2d=True):
    """Compute and summarize E,O,C,EE,OO. If k=2, also make scatter plot."""
    E, O, C, EE, OO = mitm_compute_E_O_C_EE_OO(B)
    
    print(f"Matrix B: {B.shape[0]}x{B.shape[1]}")
    print(f"  Distinct |E|={len(E)}, |O|={len(O)}")
    print(f"  Multiset sizes: sum(E)={sum(E.values())}, sum(O)={sum(O.values())}")
    print(f"  |C|={sum(C.values())}, |EE|={sum(EE.values())}, |OO|={sum(OO.values())}")
    print(f"  Sanity: |EE|+|OO|+|C| = {sum(EE.values())+sum(OO.values())+sum(C.values())} "
          f"should equal sum(E)+sum(O) = {sum(E.values())+sum(O.values())}")
    
    # Plot only if 2D
    if plot_if_2d and B.shape[1] == 2:
        EE_pts = expand_points(EE)
        OO_pts = expand_points(OO)
        C_pts  = expand_points(C)
        plt.figure(figsize=(6,6))
        if EE_pts:
            xs, ys = zip(*EE_pts)
            plt.scatter(xs, ys, color="red", label="EE", alpha=0.6)
        if OO_pts:
            xs, ys = zip(*OO_pts)
            plt.scatter(xs, ys, color="blue", label="OO", alpha=0.6)
        if C_pts:
            xs, ys = zip(*C_pts)
            plt.scatter(xs, ys, color="green", label="C", alpha=0.6)
        plt.legend()
        plt.gca().set_aspect("equal", adjustable="box")
        plt.title("Subset sum partition (EE=red, OO=blue, C=green)")
        plt.show()
    
    return E, O, C, EE, OO


def demo():
    B = Matrix([[0, -4], [2, -4], [3, -3], [4, -2], [4, 0], [4, 2], [3, 3]])
    E, O, C, EE, OO = analyze_B(B)

if __name__ == "__main__":
    demo()
    
