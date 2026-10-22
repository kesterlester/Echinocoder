#!/bin/sh

# See further down for expected output of these examples.

# Work on a random n=4, k=4 set uaing Simplex1 (a default) and overruding the lower bound for array entries (which are normally in [0,99]).
./play_simplex_encoders_via_EncDec.py n 4 k 4 random method 2 low -10

# Use Simplex_1 method:
./play_simplex_encoders_via_EncDec.py array "[[5,4,2 ], [6,3,1] ]" method 1

# Use Simplex_2 method:
./play_simplex_encoders_via_EncDec.py array "[[5,4,2 ], [6,3,1] ]" method 2

# ARCHIVE: ./play_simplex_encoders_via_EncDec.py array "[[17, 4, 28, ],  [10, 14, -7, ],  [15, 28, -8, ],  [26, 23, -2, ],  ]"
# ARCHIVE: ./play_simplex_encoders_via_EncDec.py array "[[12, 29, 27, ],  [6, -8, 0, ],  [29, 9, 18, ],  [29, 15, 6, ],  ]"
# ARCHIVE: ./play_simplex_encoders_via_EncDec.py array "[[709, 537, 355, ],  [725, 534, 223, ],  ]"

############
# When run this script should produce output resembling:
# 
# XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
# =======================
# <function simplex_2_preprocess_steps at 0x10580f2e0> encoded
# [[29, 5, 0, 25, ],  [20, -4, 21, 15, ],  [7, 8, 22, 27, ],  [-6, 28, 13, -5, ],  ]
# without canonicalisation to second diff(s)
#  +  0.0 * [[0, 0, 0, 0, ],  [0, 0, 0, 0, ],  [0, 0, 0, 0, ],  [0, 1, 0, 0, ],  ]  (one-norm 1)
#  +  7.0 * [[0, 0, 0, 1, ],  [0, 0, 0, 1, ],  [0, 0, 0, 1, ],  [0, 1, 0, 0, ],  ]  (one-norm 4)
#  +  0.0 * [[1, 0, 0, 1, ],  [1, 0, 0, 1, ],  [0, 0, 0, 1, ],  [0, 1, 0, 0, ],  ]  (one-norm 6)
#  +  0.0 * [[2, 0, 0, 1, ],  [2, 0, 0, 1, ],  [1, 0, 0, 1, ],  [0, 1, 0, 0, ],  ]  (one-norm 9)
#  +  3.0 * [[2, 0, 0, 1, ],  [2, 0, 1, 1, ],  [1, 0, 1, 1, ],  [0, 1, 1, 0, ],  ]  (one-norm 12)
#  +  1.0 * [[2, 0, 0, 2, ],  [2, 0, 1, 1, ],  [1, 0, 1, 2, ],  [0, 1, 1, 0, ],  ]  (one-norm 14)
#  +  0.0 * [[3, 0, 0, 2, ],  [2, 0, 1, 1, ],  [1, 0, 1, 2, ],  [0, 1, 1, 0, ],  ]  (one-norm 15)
#  +  1.0 * [[3, 1, 0, 2, ],  [2, 0, 1, 1, ],  [1, 1, 1, 2, ],  [0, 2, 1, 0, ],  ]  (one-norm 18)
#  +  5.0 * [[3, 1, 0, 2, ],  [2, 0, 2, 1, ],  [1, 1, 2, 2, ],  [0, 2, 1, 0, ],  ]  (one-norm 20)
#  +  1.0 * [[3, 1, 0, 2, ],  [2, 0, 2, 1, ],  [1, 2, 2, 2, ],  [0, 3, 1, 0, ],  ]  (one-norm 22)
#  +  1.0 * [[3, 1, 0, 2, ],  [2, 0, 2, 1, ],  [1, 2, 2, 3, ],  [0, 3, 1, 0, ],  ]  (one-norm 23)
#  +  1.0 * [[3, 1, 0, 2, ],  [2, 0, 2, 1, ],  [1, 2, 3, 3, ],  [0, 3, 1, 0, ],  ]  (one-norm 24)
# and offset(s)
#  +  -6.0 * [[1, 0, 0, 0, ],  [1, 0, 0, 0, ],  [1, 0, 0, 0, ],  [1, 0, 0, 0, ],  ]  (one-norm 4)
#  +  -4.0 * [[0, 1, 0, 0, ],  [0, 1, 0, 0, ],  [0, 1, 0, 0, ],  [0, 1, 0, 0, ],  ]  (one-norm 4)
#  +  0.0 * [[0, 0, 1, 0, ],  [0, 0, 1, 0, ],  [0, 0, 1, 0, ],  [0, 0, 1, 0, ],  ]  (one-norm 4)
#  +  -5.0 * [[0, 0, 0, 1, ],  [0, 0, 0, 1, ],  [0, 0, 0, 1, ],  [0, 0, 0, 1, ],  ]  (one-norm 4)
# ========================
# The first diff basis vector in the above is
# [[0, 0, 0, 0, ],  [0, 0, 0, 0, ],  [0, 0, 0, 0, ],  [0, 1, 0, 0, ],  ]  with one-norm  1
# and the differences between the subsequent (non-offset) basis vectors are:
# [[0, 0, 0, 1, ],  [0, 0, 0, 1, ],  [0, 0, 0, 1, ],  [0, 0, 0, 0, ],  ]  with one-norm  3
# [[1, 0, 0, 0, ],  [1, 0, 0, 0, ],  [0, 0, 0, 0, ],  [0, 0, 0, 0, ],  ]  with one-norm  2
# [[1, 0, 0, 0, ],  [1, 0, 0, 0, ],  [1, 0, 0, 0, ],  [0, 0, 0, 0, ],  ]  with one-norm  3
# [[0, 0, 0, 0, ],  [0, 0, 1, 0, ],  [0, 0, 1, 0, ],  [0, 0, 1, 0, ],  ]  with one-norm  3
# [[0, 0, 0, 1, ],  [0, 0, 0, 0, ],  [0, 0, 0, 1, ],  [0, 0, 0, 0, ],  ]  with one-norm  2
# [[1, 0, 0, 0, ],  [0, 0, 0, 0, ],  [0, 0, 0, 0, ],  [0, 0, 0, 0, ],  ]  with one-norm  1
# [[0, 1, 0, 0, ],  [0, 0, 0, 0, ],  [0, 1, 0, 0, ],  [0, 1, 0, 0, ],  ]  with one-norm  3
# [[0, 0, 0, 0, ],  [0, 0, 1, 0, ],  [0, 0, 1, 0, ],  [0, 0, 0, 0, ],  ]  with one-norm  2
# [[0, 0, 0, 0, ],  [0, 0, 0, 0, ],  [0, 1, 0, 0, ],  [0, 1, 0, 0, ],  ]  with one-norm  2
# [[0, 0, 0, 0, ],  [0, 0, 0, 0, ],  [0, 0, 0, 1, ],  [0, 0, 0, 0, ],  ]  with one-norm  1
# [[0, 0, 0, 0, ],  [0, 0, 0, 0, ],  [0, 0, 1, 0, ],  [0, 0, 0, 0, ],  ]  with one-norm  1
# XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
# 
# XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
# =======================
# <function simplex_2_preprocess_steps at 0x10580f2e0> encoded
# [[29, 5, 0, 25, ],  [20, -4, 21, 15, ],  [7, 8, 22, 27, ],  [-6, 28, 13, -5, ],  ]
# with canonicalisation to second diff(s)
#  +  0.0 * [[0, 0, 0, 0, ],  [0, 0, 0, 0, ],  [0, 0, 0, 0, ],  [0, 1, 0, 0, ],  ]  (one-norm 1)
#  +  7.0 * [[0, 0, 0, 1, ],  [0, 0, 0, 1, ],  [0, 0, 0, 1, ],  [0, 1, 0, 0, ],  ]  (one-norm 4)
#  +  0.0 * [[0, 0, 0, 1, ],  [0, 1, 0, 0, ],  [1, 0, 0, 1, ],  [1, 0, 0, 1, ],  ]  (one-norm 6)
#  +  0.0 * [[0, 1, 0, 0, ],  [1, 0, 0, 1, ],  [2, 0, 0, 1, ],  [2, 0, 0, 1, ],  ]  (one-norm 9)
#  +  3.0 * [[0, 1, 1, 0, ],  [1, 0, 1, 1, ],  [2, 0, 0, 1, ],  [2, 0, 1, 1, ],  ]  (one-norm 12)
#  +  1.0 * [[0, 1, 1, 0, ],  [1, 0, 1, 2, ],  [2, 0, 0, 2, ],  [2, 0, 1, 1, ],  ]  (one-norm 14)
#  +  0.0 * [[0, 1, 1, 0, ],  [1, 0, 1, 2, ],  [2, 0, 1, 1, ],  [3, 0, 0, 2, ],  ]  (one-norm 15)
#  +  1.0 * [[0, 2, 1, 0, ],  [1, 1, 1, 2, ],  [2, 0, 1, 1, ],  [3, 1, 0, 2, ],  ]  (one-norm 18)
#  +  5.0 * [[0, 2, 1, 0, ],  [1, 1, 2, 2, ],  [2, 0, 2, 1, ],  [3, 1, 0, 2, ],  ]  (one-norm 20)
#  +  1.0 * [[0, 3, 1, 0, ],  [1, 2, 2, 2, ],  [2, 0, 2, 1, ],  [3, 1, 0, 2, ],  ]  (one-norm 22)
#  +  1.0 * [[0, 3, 1, 0, ],  [1, 2, 2, 3, ],  [2, 0, 2, 1, ],  [3, 1, 0, 2, ],  ]  (one-norm 23)
#  +  1.0 * [[0, 3, 1, 0, ],  [1, 2, 3, 3, ],  [2, 0, 2, 1, ],  [3, 1, 0, 2, ],  ]  (one-norm 24)
# and offset(s)
#  +  -6.0 * [[1, 0, 0, 0, ],  [1, 0, 0, 0, ],  [1, 0, 0, 0, ],  [1, 0, 0, 0, ],  ]  (one-norm 4)
#  +  -4.0 * [[0, 1, 0, 0, ],  [0, 1, 0, 0, ],  [0, 1, 0, 0, ],  [0, 1, 0, 0, ],  ]  (one-norm 4)
#  +  0.0 * [[0, 0, 1, 0, ],  [0, 0, 1, 0, ],  [0, 0, 1, 0, ],  [0, 0, 1, 0, ],  ]  (one-norm 4)
#  +  -5.0 * [[0, 0, 0, 1, ],  [0, 0, 0, 1, ],  [0, 0, 0, 1, ],  [0, 0, 0, 1, ],  ]  (one-norm 4)
# ========================
# The first diff basis vector in the above is
# [[0, 0, 0, 0, ],  [0, 0, 0, 0, ],  [0, 0, 0, 0, ],  [0, 1, 0, 0, ],  ]  with one-norm  1
# and the differences between the subsequent (non-offset) basis vectors are:
# [[0, 0, 0, 1, ],  [0, 0, 0, 1, ],  [0, 0, 0, 1, ],  [0, 0, 0, 0, ],  ]  with one-norm  3
# [[0, 0, 0, 0, ],  [0, 1, 0, -1, ],  [1, 0, 0, 0, ],  [1, -1, 0, 1, ],  ]  with one-norm  2
# [[0, 1, 0, -1, ],  [1, -1, 0, 1, ],  [1, 0, 0, 0, ],  [1, 0, 0, 0, ],  ]  with one-norm  3
# [[0, 0, 1, 0, ],  [0, 0, 1, 0, ],  [0, 0, 0, 0, ],  [0, 0, 1, 0, ],  ]  with one-norm  3
# [[0, 0, 0, 0, ],  [0, 0, 0, 1, ],  [0, 0, 0, 1, ],  [0, 0, 0, 0, ],  ]  with one-norm  2
# [[0, 0, 0, 0, ],  [0, 0, 0, 0, ],  [0, 0, 1, -1, ],  [1, 0, -1, 1, ],  ]  with one-norm  1
# [[0, 1, 0, 0, ],  [0, 1, 0, 0, ],  [0, 0, 0, 0, ],  [0, 1, 0, 0, ],  ]  with one-norm  3
# [[0, 0, 0, 0, ],  [0, 0, 1, 0, ],  [0, 0, 1, 0, ],  [0, 0, 0, 0, ],  ]  with one-norm  2
# [[0, 1, 0, 0, ],  [0, 1, 0, 0, ],  [0, 0, 0, 0, ],  [0, 0, 0, 0, ],  ]  with one-norm  2
# [[0, 0, 0, 0, ],  [0, 0, 0, 1, ],  [0, 0, 0, 0, ],  [0, 0, 0, 0, ],  ]  with one-norm  1
# [[0, 0, 0, 0, ],  [0, 0, 1, 0, ],  [0, 0, 0, 0, ],  [0, 0, 0, 0, ],  ]  with one-norm  1
# XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
# 
# 
# XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
# =======================
# <function simplex_1_preprocess_steps at 0x105763240> encoded
# [[5, 4, 2, ],  [6, 3, 1, ],  ]
# without canonicalisation to second diff(s)
#  +  0.0 * [[0, 0, 0, ],  [1, 0, 0, ],  ]  (one-norm 1)
#  +  0.0 * [[1, 0, 0, ],  [2, 0, 0, ],  ]  (one-norm 3)
#  +  0.0 * [[2, 1, 0, ],  [3, 0, 0, ],  ]  (one-norm 6)
#  +  0.0 * [[3, 2, 0, ],  [4, 1, 0, ],  ]  (one-norm 10)
#  +  1.0 * [[4, 3, 1, ],  [5, 2, 0, ],  ]  (one-norm 15)
# and offset(s)
#  +  1.0 * [[1, 1, 1, ],  [1, 1, 1, ],  ]  (one-norm 6)
# ========================
# The first diff basis vector in the above is
# [[0, 0, 0, ],  [1, 0, 0, ],  ]  with one-norm  1
# and the differences between the subsequent (non-offset) basis vectors are:
# [[1, 0, 0, ],  [1, 0, 0, ],  ]  with one-norm  2
# [[1, 1, 0, ],  [1, 0, 0, ],  ]  with one-norm  3
# [[1, 1, 0, ],  [1, 1, 0, ],  ]  with one-norm  4
# [[1, 1, 1, ],  [1, 1, 0, ],  ]  with one-norm  5
# XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
# 
# XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
# =======================
# <function simplex_1_preprocess_steps at 0x105763240> encoded
# [[5, 4, 2, ],  [6, 3, 1, ],  ]
# with canonicalisation to second diff(s)
#  +  0.0 * [[0, 0, 0, ],  [1, 0, 0, ],  ]  (one-norm 1)
#  +  0.0 * [[1, 0, 0, ],  [2, 0, 0, ],  ]  (one-norm 3)
#  +  0.0 * [[2, 1, 0, ],  [3, 0, 0, ],  ]  (one-norm 6)
#  +  0.0 * [[3, 2, 0, ],  [4, 1, 0, ],  ]  (one-norm 10)
#  +  1.0 * [[4, 3, 1, ],  [5, 2, 0, ],  ]  (one-norm 15)
# and offset(s)
#  +  1.0 * [[1, 1, 1, ],  [1, 1, 1, ],  ]  (one-norm 6)
# ========================
# The first diff basis vector in the above is
# [[0, 0, 0, ],  [1, 0, 0, ],  ]  with one-norm  1
# and the differences between the subsequent (non-offset) basis vectors are:
# [[1, 0, 0, ],  [1, 0, 0, ],  ]  with one-norm  2
# [[1, 1, 0, ],  [1, 0, 0, ],  ]  with one-norm  3
# [[1, 1, 0, ],  [1, 1, 0, ],  ]  with one-norm  4
# [[1, 1, 1, ],  [1, 1, 0, ],  ]  with one-norm  5
# XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
# 
# 
# XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
# =======================
# <function simplex_2_preprocess_steps at 0x101aa32e0> encoded
# [[5, 4, 2, ],  [6, 3, 1, ],  ]
# without canonicalisation to second diff(s)
#  +  0.0 * [[0, 0, 0, ],  [1, 0, 0, ],  ]  (one-norm 1)
#  +  0.0 * [[0, 1, 0, ],  [1, 0, 0, ],  ]  (one-norm 2)
#  +  1.0 * [[0, 1, 1, ],  [1, 0, 0, ],  ]  (one-norm 3)
# and offset(s)
#  +  5.0 * [[1, 0, 0, ],  [1, 0, 0, ],  ]  (one-norm 2)
#  +  3.0 * [[0, 1, 0, ],  [0, 1, 0, ],  ]  (one-norm 2)
#  +  1.0 * [[0, 0, 1, ],  [0, 0, 1, ],  ]  (one-norm 2)
# ========================
# The first diff basis vector in the above is
# [[0, 0, 0, ],  [1, 0, 0, ],  ]  with one-norm  1
# and the differences between the subsequent (non-offset) basis vectors are:
# [[0, 1, 0, ],  [0, 0, 0, ],  ]  with one-norm  1
# [[0, 0, 1, ],  [0, 0, 0, ],  ]  with one-norm  1
# XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
# 
# XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
# =======================
# <function simplex_2_preprocess_steps at 0x101aa32e0> encoded
# [[5, 4, 2, ],  [6, 3, 1, ],  ]
# with canonicalisation to second diff(s)
#  +  0.0 * [[0, 0, 0, ],  [1, 0, 0, ],  ]  (one-norm 1)
#  +  0.0 * [[0, 1, 0, ],  [1, 0, 0, ],  ]  (one-norm 2)
#  +  1.0 * [[0, 1, 1, ],  [1, 0, 0, ],  ]  (one-norm 3)
# and offset(s)
#  +  5.0 * [[1, 0, 0, ],  [1, 0, 0, ],  ]  (one-norm 2)
#  +  3.0 * [[0, 1, 0, ],  [0, 1, 0, ],  ]  (one-norm 2)
#  +  1.0 * [[0, 0, 1, ],  [0, 0, 1, ],  ]  (one-norm 2)
# ========================
# The first diff basis vector in the above is
# [[0, 0, 0, ],  [1, 0, 0, ],  ]  with one-norm  1
# and the differences between the subsequent (non-offset) basis vectors are:
# [[0, 1, 0, ],  [0, 0, 0, ],  ]  with one-norm  1
# [[0, 0, 1, ],  [0, 0, 0, ],  ]  with one-norm  1
# XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX


