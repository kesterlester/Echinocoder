import sympy as sp

def size_of_confusable_multiset_formed_from(B_scaled, # B_scaled is a SCALED bad bat matrix. This means that it is a matrix with m rows, each being a k-dim direction vector which is perp to lots of good bats. The edges of the "even-odd hypercube construction" from which the alternate red and blue points are extracted to make the confusable multisets are these rows
                              ):
    """
    This algorithm is basically the same as "size_of_EE_multiset_formed_from" except that
    it throws rather than return 0. The reason is that the EE/OO construction always gives
    a multiset, but CONFUSABLE multisets are non-empty by definition, so if the caller
    asked us for a CONFUSABLE mutltiset size and we returned the size of EE we cannot be
    giving the caller what he/she wanted. Presumably this means that the caller gave us the
    wrong kind of B_scaled. Hence rather than return 0 we raise an exception.
    """

    ans = size_of_EE_multiset_formed_from(B_scaled)

    if ans>0:
        return ans
    else:
        raise ValueError(f"The Even/Odd construction resulted in total cancelation (and so no confusable multiset) when supplied with scaled bat matrix {B_scaled}.")

def size_of_EE_multiset_formed_from(B_scaled, # B_scaled is a SCALED bad bat matrix. This means that it is a matrix with m rows, each being a k-dim direction vector which is perp to lots of good bats. The edges of the "even-odd hypercube construction" from which the alternate red and blue points are extracted to make the confusable multisets are these rows
                              ):
    """
    E are the points with even numbers of edges from B.
    O are the points with odd numbers of edged from B.
    C are the points in their intersection.
    EE are the points just in E.
    OO are the points just in O.
    """

    E, O, C, EE, OO = analyze_B(B_scaled, plot_if_2d=False)

    assert len(E) == len(O)
    assert len(EE) == len(OO)
    assert len(C) + len(EE) == len(E)

    return len(EE)

def scaled_bad_bat_matrix(B_unscaled, # B_unscaled is an UNSCALED bad bat matrix. This means that it is a matrix with m rows, each being a k-dim direction vector which is perp to lots of good bats. The edges of the "even-odd hypercube construction" from which the alternate red and blue points are extracted to make the confusable multisets are scaled versions of these rows.  The lengths of these vectors are yet to be scaled by the alphas which are calculated from the null space of an L-vertex-match-matrix which has been appropriatley transformed to act on the alphas -- hence "unscaled"
                          vector_of_alphas, # vector_of_alphas is an m-dimensional vector, whose first component tells us how much to scale the first row of B, and whose second component tells us how much to scale the second row of B, and so on
                          ):

    return scale_rows(B_unscaled, vector_of_alphas)

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

    assert E.total() == O.total()
    assert EE.total() == OO.total()
    assert E.total() == EE.total() + C.total()

    return E, O, C, EE, OO

def plot_with_rings(counter, color, label, double_count = False):
    """Plot points with concentric rings for multiplicity."""
    for (x, y), mult in counter.items():
        if double_count:
            mult *= 2
        # base filled dot
        plt.scatter([x], [y], color=color, alpha=0.7, s=20, edgecolor='black')
        # add rings if mult > 1
        for r in range(1, mult):
            circle = plt.Circle((x, y), radius=0.45*r, fill=False,
                                edgecolor=color, linewidth=1.2, alpha=0.8)
            plt.gca().add_patch(circle)

def analyze_B(B: Matrix, debug=False, plot_if_2d = True, show_C_if_plotting = False):
    """Compute and summarize E,O,C,EE,OO. If k=2, also make scatter plot."""
    E, O, C, EE, OO = mitm_compute_E_O_C_EE_OO(B)

    if debug:
        print(f"Matrix B: {B.shape[0]}x{B.shape[1]}")
        print(f"  Distinct |E|={len(E)}, |O|={len(O)}")
        print(f"  Multiset sizes: sum(E)={sum(E.values())}, sum(O)={sum(O.values())}")
        print(f"  |C|={sum(C.values())}, |EE|={sum(EE.values())}, |OO|={sum(OO.values())}")
        print(f"  Sanity: |EE|+|OO|+2|C| = {sum(EE.values())+sum(OO.values())+2*sum(C.values())} "
              f"should equal sum(E)+sum(O) = {sum(E.values())+sum(O.values())}")

    # Plot only if 2D
    if plot_if_2d and B.shape[1] == 2:
        plt.figure(figsize=(6,6))
        plot_with_rings(EE, color="red", label="EE")
        plot_with_rings(OO, color="blue", label="OO")
        if show_C_if_plotting:
            plot_with_rings(C,  color="green", label="C", double_count=True)
        #plt.legend(["EE","OO","C"])
        plt.gca().set_aspect("equal", adjustable="box")
        plt.title(f"Confusable sets of size {EE.total()}\n(EE=red, OO=blue, C=green)")
        plt.show()

    return E, O, C, EE, OO

def demo():
    B = Matrix([[0, -4], [2, -4], [3, -3], [4, -2], [4, 0], [4, 2], [3, 3]])
    E, O, C, EE, OO = analyze_B(B, show_C_if_plotting = True)

    B = Matrix([[0, -4], [2, -4], [3, -3], [4, -2], [4, 0], [4, 2], [3, 3]])
    E, O, C, EE, OO = analyze_B(B, show_C_if_plotting = False)

    B = Matrix([[0, -4], [2, -4], [3, -3], [4, -2], [4, 0], [4, 2], [3, 3], [2,2]])
    E, O, C, EE, OO = analyze_B(B, show_C_if_plotting = True)

    B = Matrix([[0, -4], [2, -4], [3, -3], [4, -2], [4, 0], [4, 2], [3, 3], [2,2]])
    E, O, C, EE, OO = analyze_B(B, show_C_if_plotting = False)

if __name__ == "__main__":
    demo()

