import sympy as sp
import sympy_tools as spt

def canonical_nullspace_key(M: sp.Matrix):
    """
    Compute a canonical key representing the nullspace of M
    up to coordinate permutations.

    M: sympy.Matrix (rectangular allowed), exact rationals.

    Returns: ImmutableMatrix (hashable canonical matrix).
    """
    M = spt.lex_sort_sympy_matrix_by_cols(M)
    rref = M.rref()[0]
    stripped = spt.strip_zero_rows(rref)
    return sp.ImmutableMatrix(stripped)

### from sage.all import Graph, QQ, Matrix as SageMatrix
### 
### def canonical_nullspace_key_sage(M: sp.Matrix):
###     """
###     Canonical key for Null(M) up to permutation of coordinates.
###     Input:
###         M : sympy.Matrix with exact rational entries (SymPy Rationals)
###     Output:
###         key : tuple-of-tuples of sage.rings.rational_field.QQ entries
###               This is the canonical weighted adjacency matrix (P) under
###               simultaneous row/column permutations.
###     """
### 
###     # 1) compute exact rational nullspace basis in SymPy
###     basis = M.nullspace()
###     if not basis:  # trivial {0} nullspace
###         return (("zero", M.shape[1]),)
### 
###     # 2) build exact projection P = B * (B^T B)^{-1} * B^T in SymPy
###     B = sp.Matrix.hstack(*basis)    # n x k (columns are basis vectors)
###     Gram = B.T * B                  # k x k
###     P_sym = B * Gram.inv() * B.T    # n x n, entries are exact SymPy Rationals
### 
###     n = P_sym.rows
### 
###     # 3) convert P_sym entries to Sage QQ rationals and build labelled graph
###     #    We'll include all pairs (i,j) with their rational label (including i==j).
###     G = Graph(multiedges=False, loops=True)
###     G.add_vertices(range(n))
### 
###     # Add an edge for each pair (i,j) with a nonzero label (store label as QQ)
###     for i in range(n):
###         for j in range(i, n):
###             val_sym = P_sym[i, j]
###             # Convert SymPy Rational to Sage QQ via string (safe for arbitrary precision)
###             val_q = QQ(str(val_sym))
###             # Add edge (i,j) with label even for diagonal entries.
###             # Use label only when non-zero to keep graph sparse and canonicalizer efficient.
###             if val_q != QQ(0):
###                 G.add_edge(i, j, label=val_q)
### 
###     # 4) canonicalize the labelled graph
###     #    Use Sage's canonical_label; preserve edge labels in the canonicalization.
###     #    canonical_label returns a new graph (relabelled). We request edge_labels=True.
###     Gcanon = G.canonical_label(edge_labels=True)
### 
###     # 5) extract canonical weighted adjacency by reading edge labels from Gcanon
###     #    For each ordered pair (i,j) we query if there's an edge; if not, 0.
###     #    Because the graph is undirected, treat (i,j) same as (j,i).
###     canon_matrix = [[QQ(0) for _ in range(n)] for __ in range(n)]
###     for i in range(n):
###         for j in range(i, n):
###             if Gcanon.has_edge(i, j):
###                 # try to get label. Sage edge label retrieval can be:
###                 #   Gcanon.edge_label((i,j))  or  Gcanon.edge_label(i,j)
###                 # We'll attempt both (the environment may accept one).
###                 try:
###                     lab = Gcanon.edge_label((i, j))
###                 except Exception:
###                     lab = Gcanon.edge_label(i, j)
###                 # lab should already be a QQ; if not, coerce:
###                 lab_q = QQ(lab)
###                 canon_matrix[i][j] = lab_q
###                 canon_matrix[j][i] = lab_q
### 
###     # 6) turn into a hashable tuple-of-tuples and return
###     key = tuple(tuple(canon_matrix[i][j] for j in range(n)) for i in range(n))
###     return key


## import sympy as sp
## 
## def old_wrong_canonical_nullspace_key(M: sp.Matrix):
##     """
##     Compute a canonical key representing the nullspace of M
##     up to coordinate permutations.
## 
##     M: sympy.Matrix (rectangular allowed), exact rationals.
## 
##     Returns: tuple of tuples (hashable canonical matrix).
##     """
## 
##     # Step 1: Compute a basis for the nullspace (exact rational vectors).
##     # SymPy returns a list of column vectors (as sympy.Matrix objects).
##     basis = M.nullspace()
## 
##     if not basis:  # trivial nullspace {0}
##         # Encode by a special marker + ambient dimension
##         return (("zero", M.shape[1]),)
## 
##     # Step 2: Build projection matrix P onto the nullspace.
##     # If B is the n×k matrix whose columns are a basis of the nullspace,
##     # then the projection is P = B * (B^T * B)^(-1) * B^T.
##     # This works over rationals without needing orthonormalization.
##     B = sp.Matrix.hstack(*basis)      # n×k basis matrix
##     Gram = B.T * B                    # k×k Gram matrix
##     P = B * Gram.inv() * B.T          # n×n projection matrix, exact
## 
##     # Step 3: Canonicalize under coordinate permutation.
##     #   - Represent each row as a tuple of rationals.
##     #   - Sort rows lexicographically to choose a canonical row order.
##     #   - Apply the same permutation to columns (to preserve symmetry).
##     rows = [tuple(P.row(i)) for i in range(P.rows)]
##     idx = sorted(range(len(rows)), key=lambda i: rows[i])
##     P_perm = P.extract(idx, idx)
## 
##     # Step 4: Serialize into a hashable key (tuple of tuples).
##     key = tuple(tuple(P_perm[i, j] for j in range(P_perm.cols))
##                 for i in range(P_perm.rows))
## 
##     return key

