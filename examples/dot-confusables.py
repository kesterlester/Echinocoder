#!/usr/bin/env python

this is slow/broken

"""
We work with (multu)sets of $n$ vectors in $k$-dimensions.

Encoding dots are done with $m$ directions: $D=\{\vec d_1, \vec d_2, ... , \vec d_m\}$

There are $M$ vectors ($M=ceil(m/(k-1)$) which are orthogonal to the vectors in $D$ (these are the $M$ bad-bats): $B=\{\vec q_1, \vec q_2, \cdots, \vec q_M\}$. E.g. $\vec q_1 \cdot m_i = 0$ for all $i\in \{ 1, 2, \cdots, k-1 \}$.

The (multi)sets $E$ and $O$ have even and odd sums of the elements of $B$:

$E=\{ \vec 0, \vec q_1 + \vec q_2, \cdots \}$
$O=\{  \vec q_1, \vec q_2, \cdots, \vec q_2 + \vec q_5 + \vec q_9, \cdots \}$

We can represent elements of E and O (and their later trimmed version) by a set of dimension-M vectors (one for each free parameter) whose components show us how much of each cpt is to be used.

E.g. if the coefficients are lambda_1, lambda_2, ... \lambda_M which are buried inside each of the q_i .... i.e.:

     \vec q_1 = \lambda_1(parameter_1) \vec Q_1
     \vec q_2 = \lambda_2(parameter_4) \vec Q_2
     \vec q_3 = \lambda_3(parameter_1) \vec Q_3
     \vec q_4 = \lambda_4(parameter_5) \vec Q_4
     \vec q_5 = \lambda_5(parameter_2) \vec Q_5
     \vdots

then the vector

    5 q_1 + 2 q_2 + 3 q_3 - 8 q_4 + 9 q_5

would be represented by

    v_{param 1} = [  5,  0,  3,  0,  0, ... ]
    v_{param 2} = [  0,  0,  0,  0,  9, ... ]
    v_{param 3} = [  0,  0,  0,  0,  0, ... ]
    v_{param 4} = [  0,  2,  0,  0,  0, ... ]
    v_{param 5} = [  0,  0,  0, -8,  0, ... ]
    \vdots

which we could put into a matrix V_{ij}, with index i being the parameter, and index j being the occurrence count of vec q_j.

With such a format in mind, we could identify whether two different linear combinations(A and B) could be made equal by checking whether their difference (A-B) has at least k+1 non-zero rows ... and in that case 

"""

import os, sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))  # the repository root, for its modules

import numpy as np
from itertools import combinations
from tuple_ize import tuple_ize

# Set small dimensions for illustration
t, M, k = 4, 3, 4
t, M, k = 4, 3, 1

# Create dummy data for n
# Shape: (2, t, M) => (2, 4, 3)
# A = n[0], B = n[1]
n = np.array([
    [  # A: Two example matrices of shape (4,3)
        [ [1, 2, 3],
          [4, 5, 6],
          [7, 8, 9],
          [10, 11, 12] ],

     #   [ [1, 0, 0],
     #     [0, 1, 0],
     #     [0, 0, 1],
     #     [1, 1, 1] ],

     #   [ [1, 2, 3],
     #     [4, 5, 6],
     #     [7, 8, 9],
     #     [0, 0, 0] ]
    ],
    [  # B: Two example matrices of shape (4,3)
        [ [0, 0, 0],
          [0, 0, 0],
          [0, 0, 0],
          [0, 0, 0] ],

     #    [ [1, 0, 0],
     #      [0, 0, 0],
     #      [0, 0, 1],
     #      [1, 1, 0] ],

     #    [ [1, 1, 1],
     #      [1, 1, 1],
     #      [1, 1, 1],
     #      [0, 0, 0] ]
    ]
])

# Generator function:
def generate_variants(n, k):
    """
    n is numpy array with shape (2, something, t, M)

    This generator treats n[0] as equivalent to a list of txM matrices (A), 
    and n[1] as equivalent to another list of txM matrices (B).
    Then for each a in A and each b in B it constructs the matrix d=a-b ...
    ... i.e. it loops over the differences of all pairs of matrices, with one
    from A and one from B.  And then it yields successively all "variants" of d ... these
    being zero matrices into which k+1 non-zero rows of d have been copied, provided that
    there are at least k+1 rows in d.

    The purpose of doing so, is that the 

    """
    A = n[0]  # shape (something, t, M)
    B = n[1]  # shape (something, t, M)

    seen_differences = set()
    for a in A:
        for b in B:
            d = a - b  # shape (t, M)
            d_tup = tuple_ize(d) 

            if d_tup in seen_differences:
                # Don't need to use it again:
                continue
            else:
                # Record it (and its negative) so we don't use it again:
                seen_differences.add(d_tup)
                seen_differences.add(tuple_ize(-d))

            non_zero_row_indices = np.where(~np.all(d == 0, axis=1))[0]

            if len(non_zero_row_indices) < k + 1:
                # The only way of setting the difference represented by d to zero would 
                # necessarily set every member zero, and that
                # would "collapse" E into O. So ignore and continue:
                continue

            for rows_to_keep in combinations(non_zero_row_indices, k + 1):
                variant = np.zeros_like(d)
                variant[list(rows_to_keep)] = d[list(rows_to_keep)]
                yield variant

# Run and print results
for i, variant in enumerate(generate_variants(n, k)):
    print(f"Variant {i+1}:\n{variant}\n")

