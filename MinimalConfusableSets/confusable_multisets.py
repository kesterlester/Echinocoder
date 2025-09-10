import sympy as sp
import sympy_tools as spt
from math import sqrt

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

def size_of_EE_multiset_formed_from(scaled_bad_bat_matrix,  # scaled_bat_matrix is a SCALED bad bat matrix. This means that it is a matrix with m rows, each being a k-dim direction vector which is perp to lots of good bats. The edges of the "even-odd hypercube construction" from which the alternate red and blue points are extracted to make the confusable multisets are these rows
                                    ):
    """
    E are the points with even numbers of edges from B.
    O are the points with odd numbers of edged from B.
    C are the points in their intersection.
    EE are the points just in E.
    OO are the points just in O.
    """

    E, O, C, EE, OO = mitm_compute_E_O_C_EE_OO(scaled_bad_bat_matrix)

    assert E.total() == O.total()
    assert EE.total() == OO.total()
    assert C.total() + EE.total() == E.total()

    return EE.total()

def scaled_bad_bat_matrix(B_unscaled, vector_of_alphas,):
    """
    Args:
        B_unscaled: this is an UNSCALED bad bat matrix. This means that it is a matrix with m rows,
                     each being a k-dim direction vector which is perp to lots of good bats.
                     The edges of the "even-odd hypercube construction" from which the alternate
                     red and blue points are extracted to make the confusable multisets are scaled
                     versions of these rows.  The lengths of these vectors are yet to be scaled by
                     the alphas which are calculated from the null space of an
                     L-vertex-match-matrix which has been appropriatley transformed to act on the
                     alphas -- hence "unscaled"

        vector_of_alphas: this is an m-dimensional vector, whose first component tells us how
                     much to scale the first row of B_unscaled, and whose second component tells us
                     how much to scale the second row of B_unscaled, and so on.

    Returns: B_unscaled with its rows scaled by the coefficients in vector_of_alphas.
    """

    return spt.scale_rows(B_unscaled, vector_of_alphas)

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

def mitm_compute_E_O_C_EE_OO(scaled_bad_bat_matrix: Matrix):
    """Meet-in-the-middle computation of E, O, C, EE, OO for scaled_bad_bat_matrix."""
    m, k = scaled_bad_bat_matrix.shape
    rows = _row_tuple_rows([scaled_bad_bat_matrix.row(i) for i in range(m)])
    if m>=2: #2: # could change 2 to a number bigger than 2, but not less than 2 as otherwise split leads ot empty row
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
    else:
        E, O = subset_sums_with_parity(rows)
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


def confusable_sets_or_None(L_matrix : sp.Matrix, unscaled_bad_bat_matrix:sp.Matrix, M:int):
    """
    If no scaling of the bad-bat lattice achieving the matches in L_matrix generates confusable sets,
    return None.
    Else, return a tuple containing:
        * EE and OO (a pair of confusable sets obtained from the bad-bat lattice using the L_matrix matches)
        * the length-M vector of scalings that would scale the M unscaled_bad_bats to the M scaled_bad_bats that made EE and OO
        * the scaled_bad_bats (this is really a convenience output, as it is would be easy to recompute via:

                import confusable_multisets
                scaled_bad_bat_matrix = confusable_multisets.scaled_bad_bat_matrix(unscaled_bad_bat_matrix, point_in_null_space)

            so I may remove it in future. However, since we already have it to hand, we may as well output it for now.
    """

    vote_for_collapse, null_space = vote_for_collapse_and_null_space(L_matrix, unscaled_bad_bat_matrix, M)

    if vote_for_collapse:
        assert 0 <= len(null_space) <= M # Note <= not < in first inequality.
        return None

    assert not vote_for_collapse
    assert 0 < len(null_space) <= M # Note < not <= in first inequality.

    # OK - it is now our job to generate some confusable sets!
    import nonzero_lin_comb

    null_space_contribs, point_in_null_space =  nonzero_lin_comb.combine_many(null_space)

    assert len(null_space_contribs) == len(null_space)
    assert point_in_null_space.shape == (M, 1)
    assert all(alpha!=0 for alpha in point_in_null_space)

    import confusable_multisets

    scaled_bad_bat_matrix = confusable_multisets.scaled_bad_bat_matrix(unscaled_bad_bat_matrix, point_in_null_space)

    ##  E, O, C, EE, OO = confusable_multisets.analyze_B(scaled_bad_bat_matrix, plot_if_2d=False, show_C_if_plotting = False)
    _, _, _, EE, OO = confusable_multisets.mitm_compute_E_O_C_EE_OO(scaled_bad_bat_matrix)

    assert EE.total() == OO.total(), f"Must have {EE.total()}=={OO.total()} when scaled_bad_bat_matrix = {scaled_bad_bat_matrix}"

    return (EE, OO, point_in_null_space, scaled_bad_bat_matrix)

def vote_for_collapse_and_null_space(L_matrix: sp.Matrix, bat_matrix: sp.Matrix, M: int) -> tuple:
        """
        A key returned element is the null-space. It is probably only interesting when the
        null space is present.  The null space is a basis
        for the space of set of scales (alpha_1, alpha_2, ..., alpha_M) which, if applied
        to the unscaled bad-bat directions, would lead the bad-bat lattice to exhibit the
        matches specified in L_matrix.  A general point in the null space will always
        result in the L_matrix matches working, but might (sometimes) also lead to
        bad bat lattice collapse, e.g. as would happen if at least one of the alphas was
        zero.  However, "vote for collapse" being False will assure us that there is
        at least one point in the null space where non-collapse happens.

        In principle M could be found as the number of columns of L_matrix, or as the number
        of rows of bat_matrix, but this can fail for a few special cases, such as when
        L_matrix has no rows at all, etc. So to be on the safe side we pass in M
        """
        if L_matrix.rows != 0 and L_matrix.cols != M:
            raise ValueError()

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
            return True, null_space

        if null_space_dimension == M:
            # The null-space is M-dimensional. There are only M alphas, so it is capable of spanning to
            # (alpha1, alpha2, .... ) = (1, 1, ... ).
            # So no collapse here!
            return False, null_space

        # OK - the shortcuts didn't work, so fall back to the default method that will always work.
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

        for cpt_index in range(M):
            if all( (basis_vec[cpt_index, 0] == 0) for basis_vec in null_space  ):
                # There are zeros in every basis vector at this cpt,
                # so there is no soln with every alphai nonzero,
                # so this matrix collapses the bad bat lattice, and
                # so no coonfusable sets were found here
                return True, null_space
            # Try next cpt

        # OK - if we got here we tried all components but didn't find one that had zero in every basis vector,
        # so this matrix does admit solns with all alphai non-zero, i.e. the matrix does not collapse.
        # See proof in "OneNote -> Research -> Symmetries -> Non collapsing null space".
        return False, null_space


def alpha_attacking_matrix(
                L_matrix : sp.Matrix,  # Any number of rows, but M columns.
                bat_matrix : sp.Matrix, # M rows, and k columns, so that each row is a bat vector.
                ) -> sp.Matrix :
    """
    This method generates the matrix A for which the solns of A.(vec of alphas) are the same as the solutions to L.(alpha1 w1, alpha2 w2, ... , alphaM, wM) where w1 is the first bat (i.e. first row of bat matrix) and w2 the second, and so on.
    """

    M, k = bat_matrix.shape
    R, M_L = L_matrix.shape

    assert M == M_L or R==0, f"L and B must work with the same no. of vectors. Wanted M({M}) == M_L({M_L})"

    effective_row = 0
    ans = sp.zeros(R*k, M)
    for i in range(R):
        for kk in range(k):
            for j in range(M):
                ans[effective_row, j] = L_matrix[i,j]*bat_matrix[j,kk]
            effective_row += 1

    return ans


# Plotting and cosmetics appear after this line.



def old_plot_with_rings(counter, color, label, double_count = False):
    """Plot points with concentric rings for multiplicity."""
    if counter.total()==0:
        # nothing to plot
        return
    max_x = max(x for (x,_),_ in counter.items())
    min_x = min(x for (x,_),_ in counter.items())
    max_y = max(y for (_,y),_ in counter.items())
    min_y = min(y for (_,y),_ in counter.items())

    scale = max( (max_x - min_x, max_y-min_y) )

    for (x, y), mult in counter.items():
        if double_count:
            mult *= 2
        # base filled dot
        plt.scatter([x], [y], color=color, alpha=0.7, s=20, edgecolor='black')
        # add rings if mult > 1
        for r in range(1, mult):
            circle = plt.Circle((x, y), radius=0.015*sqrt(r)*scale, fill=False,
                                edgecolor=color, linewidth=0.8, alpha=0.8)
            plt.gca().add_patch(circle)

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D  # registers 3D projection

def plot_with_rings(counter, color, label, double_count=False):
    """
    Plot multiset points with multiplicities.
    - For 2D: solid dot + concentric rings.
    - For 3D: marker size encodes multiplicity.
    - For higher dims: project to first 2 components.
    """
    if not counter:
        return

    if counter.total()==0:
        # nothing to plot
        return

    dim = len(next(iter(counter.keys())))
    ax = plt.gca()

    if dim==2:
        max_x = max(x for (x, _), _ in counter.items())
        min_x = min(x for (x, _), _ in counter.items())
        max_y = max(y for (_, y), _ in counter.items())
        min_y = min(y for (_, y), _ in counter.items())
        scale = max((max_x - min_x, max_y - min_y))
    elif dim==3:
        max_x = max(x for (x, _, _), _ in counter.items())
        min_x = min(x for (x, _, _), _ in counter.items())
        max_y = max(y for (_, y, _), _ in counter.items())
        min_y = min(y for (_, y, _), _ in counter.items())
        max_z = max(z for (_, _, z), _ in counter.items())
        min_z = min(z for (_, _, z), _ in counter.items())
        scale = max((max_x - min_x, max_y - min_y, max_z-min_z))
    else:
        raise NotImplementedError("Ranging for dims greater than 3 is not yet implemented.")

    if dim == 2:
        # --- 2D mode with rings ---
        for (x, y), mult in counter.items():
            if double_count: mult*=2
            x, y = float(x), float(y)
            ax.scatter([x], [y], color=color, alpha=0.7, s=20,
                       edgecolor='black', label=label)
            label = None  # only use label once
            for r in range(1, mult):
                circ = plt.Circle((x, y), radius=0.015*sqrt(r)*scale, fill=False,
                                  edgecolor=color, linewidth=1.2, alpha=0.8)
                ax.add_patch(circ)

    elif dim == 3:
        # --- 3D mode ---
        if ax.name != "3d":
            # upgrade axes to 3D
            ax = plt.gcf().add_subplot(111, projection="3d")

        xs, ys, zs, sizes = [], [], [], []
        for (x, y, z), mult in counter.items():
            if double_count: mult*=2
            xs.append(float(x))
            ys.append(float(y))
            zs.append(float(z))
            sizes.append(float(40 + 30*(mult-1)))
        ax.scatter(xs, ys, zs, color=color, alpha=0.7, s=sizes,
                   edgecolor='black', label=label)

    else:
        # --- fallback: project higher dims to first 2 coords ---
        print(f"[plot_with_rings] Warning: dim={dim}>3, projecting to 2D.")
        for pt, mult in counter.items():
            if double_count: mult*=2
            x, y = float(pt[0]), float(pt[1])
            ax.scatter([x], [y], color=color, alpha=0.7, s=20,
                       edgecolor='black', label=label)
            label = None
            for r in range(1, mult):
                circ = plt.Circle((x, y), radius=0.015*sqrt(r)*scale, fill=False,
                                  edgecolor=color, linewidth=1.2, alpha=0.8)
                ax.add_patch(circ)



def analyze_B(scaled_bad_bat_matrix: Matrix, title_add ="",
              debug=False, try_to_plot = True,
              show_C_if_plotting = False):
    """Compute and summarize E,O,C,EE,OO. If k=2, also make scatter plot."""
    E, O, C, EE, OO = mitm_compute_E_O_C_EE_OO(scaled_bad_bat_matrix)

    if debug:
        print(f"Scaled bad bat matrix: {scaled_bad_bat_matrix.shape[0]}x{scaled_bad_bat_matrix.shape[1]}")
        print(f"  Distinct |E|={len(E)}, |O|={len(O)}")
        print(f"  Multiset sizes: sum(E)={sum(E.values())}, sum(O)={sum(O.values())}")
        print(f"  |C|={sum(C.values())}, |EE|={sum(EE.values())}, |OO|={sum(OO.values())}")
        print(f"  Sanity: |EE|+|OO|+2|C| = {sum(EE.values())+sum(OO.values())+2*sum(C.values())} "
              f"should equal sum(E)+sum(O) = {sum(E.values())+sum(O.values())}")

    # Plot only if 2D
    if try_to_plot and scaled_bad_bat_matrix.shape[1] >= 1:
        plt.figure(figsize=(6,6))
        #print("Plotting EE")
        plot_with_rings(EE, color="red", label="EE")
        #print("Plotting OO")
        plot_with_rings(OO, color="blue", label="OO")
        if show_C_if_plotting:
            #print("Plotting C")
            plot_with_rings(C,  color="green", label="C", double_count=True)
        #plt.legend(["EE","OO","C"])
        plt.gca().set_aspect("equal", adjustable="box")
        plt.title(title_add + f"\nConfusable sets of size {EE.total()}\n(EE=red, OO=blue, C=green)")
        plt.show()

    return E, O, C, EE, OO

def demo():
    ImmutableDenseMatrix = sp.ImmutableDenseMatrix
    Rational = sp.Rational
    Integer = sp.Integer

    #scaled_bad_bats_M7k3 = sp.ImmutableDenseMatrix([[Rational(846339635568568320615535800, 39552309291267152811577), Rational(-211895659946238530725168272, 39552309291267152811577), Rational(-570977642691572355566887176, 39552309291267152811577)], [Rational(-1660721099925955406255386400, 118656927873801458434731), Rational(2934958208582377808040309715, 118656927873801458434731), Rational(-1089848221826408235355097325, 39552309291267152811577)], [Rational(10168732570317317945339997515, 237313855747602916869462), Rational(1703244585904896448435423675, 39552309291267152811577), Rational(10800745029427560372163695151, 237313855747602916869462)], [Rational(-4064330700823362155565441489, 79104618582534305623154), Rational(-12382903778792177825228586335, 118656927873801458434731), Rational(6751801748311561961771984435, 237313855747602916869462)], [Rational(133831959296634705086942476, 118656927873801458434731), Rational(4973898792333826264057510411, 118656927873801458434731), Rational(-1264598598438539798067295430, 39552309291267152811577)], [Rational(-60360917979435175285391912, 39552309291267152811577), Rational(127747876473061361453205057, 39552309291267152811577), Rational(3670596363614301199787346, 39552309291267152811577)], [Integer(4748), Integer(-19311), Integer(-9925)]])
    #E, O, C, EE, OO = analyze_B(scaled_bad_bats_M7k3, show_C_if_plotting = True)

    B = Matrix([[0, -4], [2, -4], [3, -3], [4, -2], [4, 0], [4, 2], [3, 3]])
    E, O, C, EE, OO = analyze_B(B, show_C_if_plotting = True)

    B = Matrix([[0, -4, 1], [2, -4, 2], [3, -3,3], [4, -2,2],])
    E, O, C, EE, OO = analyze_B(B, show_C_if_plotting = True)

    
    if False:
        B = Matrix([[0, -4], [2, -4], [3, -3], [4, -2], [4, 0], [4, 2], [3, 3]])
        E, O, C, EE, OO = analyze_B(B, show_C_if_plotting = False)

        B = Matrix([[0, -4], [2, -4], [3, -3], [4, -2], [4, 0], [4, 2], [3, 3], [2,2]])
        E, O, C, EE, OO = analyze_B(B, show_C_if_plotting = True)

        B = Matrix([[0, -4], [2, -4], [3, -3], [4, -2], [4, 0], [4, 2], [3, 3], [2,2]])
        E, O, C, EE, OO = analyze_B(B, show_C_if_plotting = False)



        # M=9, k=2 update:
        #VSLW sort=False SO FAR:
        #VSLW sort=False SO FAR:  for M=9, k=2 the smallest confusable sets have size 99,
        #VSLW sort=False SO FAR:  raw=
        #VSLW sort=False SO FAR:  Matrix([
        #VSLW sort=False SO FAR:  [-1, -1, -1, -1, -1, -1, -1, -1, -1],
        #VSLW sort=False SO FAR:  [-1, -1, -1, -1, -1,  0,  0,  0,  0],
        #VSLW sort=False SO FAR:  [-1, -1,  0,  0,  0, -1, -1, -1,  0],
        #VSLW sort=False SO FAR:  [ 0,  0, -1, -1,  0, -1, -1,  0, -1]]),
        #VSLW sort=False SO FAR:  rre=
        #VSLW sort=False SO FAR:  Matrix([
        #VSLW sort=False SO FAR:  [1, 1, 0, 0, 0, 0, 0,  0, -1],
        #VSLW sort=False SO FAR:  [0, 0, 1, 1, 0, 0, 0, -1,  0],
        #VSLW sort=False SO FAR:  [0, 0, 0, 0, 1, 0, 0,  1,  1],
        #VSLW sort=False SO FAR:  [0, 0, 0, 0, 0, 1, 1,  1,  1]]),
        #VSLW sort=False SO FAR:  unscaled_bad_bats=
        #VSLW sort=False SO FAR:  Matrix([
        #VSLW sort=False SO FAR:  [-11575,  2898],
        #VSLW sort=False SO FAR:  [  7809,  5440],
        #VSLW sort=False SO FAR:  [ -9614, 10710],
        #VSLW sort=False SO FAR:  [  7015,  7050],
        #VSLW sort=False SO FAR:  [  7451, 11043],
        #VSLW sort=False SO FAR:  [ 22430, -6115],
        #VSLW sort=False SO FAR:  [   472, 17542],
        #VSLW sort=False SO FAR:  [-13380,  3256],
        #VSLW sort=False SO FAR:  [ -6891,  -198]]),
        #VSLW sort=False SO FAR:  scalings=MutableDenseMatrix([[Rational(17970429, 42799241)], [Rational(-11130984, 42799241)], [Rational(-14572415930046, 40971109326821)], [Rational(69645109659177, 204855546634105)], [Rational(6271584, 43003949)], [Rational(65586733880272, 1420396568278305)], [Rational(30653135591672, 284079313655661)], [Rational(-74622015, 172015796)], [Integer(1)]])
        #VSLW sort=False SO FAR:  scalingsSREPR=MutableDenseMatrix([[Rational(17970429, 42799241)], [Rational(-11130984, 42799241)], [Rational(-14572415930046, 40971109326821)], [Rational(69645109659177, 204855546634105)], [Rational(6271584, 43003949)], [Rational(65586733880272, 1420396568278305)], [Rational(30653135591672, 284079313655661)], [Rational(-74622015, 172015796)], [Integer(1)]])
        #VSLW sort=False SO FAR:  scaled_bad_bats=
        #VSLW sort=False SO FAR:  ImmutableDenseMatrix([[Rational(-208007715675, 42799241), Rational(52078303242, 42799241)], [Rational(-86921854056, 42799241), Rational(-60552552960, 42799241)], [Rational(6091269858759228, 1781352579427), Rational(-156070574610792660, 40971109326821)], [Rational(4248351689209797, 1781352579427), Rational(98199604619439570, 40971109326821)], [Rational(46729572384, 43003949), Rational(69257102112, 43003949)], [Rational(294222088186900192, 284079313655661), Rational(-80212575535572656, 284079313655661)], [Rational(14468279999269184, 284079313655661), Rational(537717304549110224, 284079313655661)], [Rational(249610640175, 43003949), Rational(-60742320210, 43003949)], [Integer(-6891), Integer(-198)]]).
        #VSLW sort=False SO FAR:  212 matrices have been scanned.
        #VSLW sort=False SO FAR:
        #VSLW sort=False SO FAR:


        # M=9, k=2
        B = ImmutableDenseMatrix([[Rational(-208007715675, 42799241), Rational(52078303242, 42799241)], [Rational(-86921854056, 42799241), Rational(-60552552960, 42799241)], [Rational(6091269858759228, 1781352579427), Rational(-156070574610792660, 40971109326821)], [Rational(4248351689209797, 1781352579427), Rational(98199604619439570, 40971109326821)], [Rational(46729572384, 43003949), Rational(69257102112, 43003949)], [Rational(294222088186900192, 284079313655661), Rational(-80212575535572656, 284079313655661)], [Rational(14468279999269184, 284079313655661), Rational(537717304549110224, 284079313655661)], [Rational(249610640175, 43003949), Rational(-60742320210, 43003949)], [Integer(-6891), Integer(-198)]])
        for x in True, False:
            E, O, C, EE, OO = analyze_B(B, show_C_if_plotting = x, title_add = "M=9, k=2")

        # M=13, k=5
        B = ImmutableDenseMatrix([[Rational(169287, 94868), Integer(0), Rational(56429, 94868), Rational(56429, 94868), Rational(-56429, 94868)], [Rational(43929, 47434), Integer(0), Rational(-131787, 94868), Rational(43929, 94868), Rational(-43929, 94868)], [Rational(121505, 94868), Rational(-121505, 189736), Rational(121505, 94868), Rational(-364515, 189736), Rational(121505, 47434)], [Rational(-94233, 47434), Rational(-94233, 94868), Rational(94233, 94868), Rational(471165, 94868), Integer(0)], [Integer(0), Rational(117515, 189736), Rational(-23503, 94868), Rational(23503, 189736), Rational(-23503, 94868)], [Rational(-71035, 23717), Rational(71035, 23717), Rational(71035, 94868), Rational(-213105, 94868), Rational(-213105, 94868)], [Rational(23489, 23717), Rational(-46978, 23717), Rational(-46978, 23717), Rational(-46978, 23717), Rational(23489, 23717)], [Rational(-29055, 23717), Integer(0), Integer(0), Rational(-29055, 23717), Rational(58110, 23717)], [Integer(0), Rational(48045, 23717), Rational(19218, 23717), Rational(-19218, 23717), Rational(-9609, 23717)], [Rational(-141837, 94868), Rational(47279, 23717), Rational(-47279, 47434), Rational(47279, 47434), Rational(-47279, 47434)], [Integer(-5), Integer(-2), Integer(-1), Integer(-4), Integer(-2)], [Integer(4), Integer(2), Integer(1), Integer(-1), Integer(-1)], [Integer(-1), Integer(6), Integer(1), Integer(2), Integer(-1)]])
        for x in True, False:
            E, O, C, EE, OO = analyze_B(B, show_C_if_plotting = x, title_add = "M=13, k=5")

        # M=11, k=3
        #     VSLW:     478445:  raw=Matrix([
        #                                 [-1, -1, -1, -1, -1,  0,  0,  0,  0, 0, 0],
        #                                 [-1, -1,  0,  0,  0, -1, -1, -1,  0, 0, 0], 
        #                                 [ 0,  0, -1, -1,  0,  0,  0,  1, -1, 0, 1],
        # ]), rre=Matrix([
        # [1, 1, 0, 0, 0,  1,  1,  1,  0, 0,  0],
        # [0, 0, 1, 1, 0,  0,  0, -1,  1, 0, -1],
        # [0, 0, 0, 0, 1, -1, -1,  0, -1, 0,  1]]), EE.total()=824, OO.total()=824
        # VSLW: The smallest confusable sets so far have scaled bad bats ImmutableDenseMatrix([[Rational(-2504245519491118249040122727900630072700, 341859150003571766011271429186537353), Rational(89568695038386432406273072082147805624, 48837021429081680858753061312362479), Rational(1689473283948694808358904395868338681444, 341859150003571766011271429186537353)], [Rational(551555026084548715703548960439332083520, 48837021429081680858753061312362479), Rational(-974751842054568263377558769423481369662, 48837021429081680858753061312362479), Rational(1085873957603955284041362015864935039430, 48837021429081680858753061312362479)], [Rational(-525905350930859629933716378717151187755, 48837021429081680858753061312362479), Rational(-528529255033864631651133352809111314850, 48837021429081680858753061312362479), Rational(-558591699185436222756396398834140199567, 48837021429081680858753061312362479)], [Rational(311061213857843453602078447695296853447, 48837021429081680858753061312362479), Rational(631812281701659754079020155918274782470, 48837021429081680858753061312362479), Rational(-172248421872744065813339645717354003335, 48837021429081680858753061312362479)], [Rational(3980251254334877660452581675088806168, 9239436486583020703007335923960469), Rational(21132435685091532663335105249518110714, 1319919498083288671858190846280067), Rational(-112830003777543777747575302569254717220, 9239436486583020703007335923960469)], [Rational(-403503010695697194068345911624061177344, 9239436486583020703007335923960469), Rational(853973970117951278969585895884952571584, 9239436486583020703007335923960469), Rational(24537345245008613152804818950111828352, 9239436486583020703007335923960469)], [Rational(176486952942893149836009762137099177064, 9239436486583020703007335923960469), Rational(-717805296604930416276997581429975191298, 9239436486583020703007335923960469), Rational(-368920178592715777616343068494252176150, 9239436486583020703007335923960469)], [Rational(1006136353393289553673102505120414214060, 48837021429081680858753061312362479), Rational(165434444161642985310461749507881696812, 48837021429081680858753061312362479), Rational(493081978098397612643210960887186987296, 48837021429081680858753061312362479)], [Rational(33151263457585138467662082113329465769, 1319919498083288671858190846280067), Rational(-13198364542020563344791609749033287408, 1319919498083288671858190846280067), Rational(37496746212965251001275500044625902903, 1319919498083288671858190846280067)], [Integer(-4062), Integer(2714), Integer(398)], [Integer(115), Integer(-11272), Integer(3347)]]).
        # VSLW: The smallest confusable sets so far have 824==824 points.

        #########################################################

        # VSLW:  M=11 and k=3 smallest confusable set was size 824 and was found after checking 1986006 match matrices. It has scaled bad bats ImmutableDenseMatrix([[Rational(-2504245519491118249040122727900630072700, 341859150003571766011271429186537353), Rational(89568695038386432406273072082147805624, 48837021429081680858753061312362479), Rational(1689473283948694808358904395868338681444, 341859150003571766011271429186537353)], [Rational(551555026084548715703548960439332083520, 48837021429081680858753061312362479), Rational(-974751842054568263377558769423481369662, 48837021429081680858753061312362479), Rational(1085873957603955284041362015864935039430, 48837021429081680858753061312362479)], [Rational(-525905350930859629933716378717151187755, 48837021429081680858753061312362479), Rational(-528529255033864631651133352809111314850, 48837021429081680858753061312362479), Rational(-558591699185436222756396398834140199567, 48837021429081680858753061312362479)], [Rational(311061213857843453602078447695296853447, 48837021429081680858753061312362479), Rational(631812281701659754079020155918274782470, 48837021429081680858753061312362479), Rational(-172248421872744065813339645717354003335, 48837021429081680858753061312362479)], [Rational(3980251254334877660452581675088806168, 9239436486583020703007335923960469), Rational(21132435685091532663335105249518110714, 1319919498083288671858190846280067), Rational(-112830003777543777747575302569254717220, 9239436486583020703007335923960469)], [Rational(-403503010695697194068345911624061177344, 9239436486583020703007335923960469), Rational(853973970117951278969585895884952571584, 9239436486583020703007335923960469), Rational(24537345245008613152804818950111828352, 9239436486583020703007335923960469)], [Rational(176486952942893149836009762137099177064, 9239436486583020703007335923960469), Rational(-717805296604930416276997581429975191298, 9239436486583020703007335923960469), Rational(-368920178592715777616343068494252176150, 9239436486583020703007335923960469)], [Rational(1006136353393289553673102505120414214060, 48837021429081680858753061312362479), Rational(165434444161642985310461749507881696812, 48837021429081680858753061312362479), Rational(493081978098397612643210960887186987296, 48837021429081680858753061312362479)], [Rational(33151263457585138467662082113329465769, 1319919498083288671858190846280067), Rational(-13198364542020563344791609749033287408, 1319919498083288671858190846280067), Rational(37496746212965251001275500044625902903, 1319919498083288671858190846280067)], [Integer(-4062), Integer(2714), Integer(398)], [Integer(115), Integer(-11272), Integer(3347)]]).


if __name__ == "__main__":
    demo()

