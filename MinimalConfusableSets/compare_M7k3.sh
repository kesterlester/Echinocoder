#!/bin/sh

python M7k3.py '[[-1, -1, -1, -1, -1,  0, 0], 
                 [-1, -1,  0,  0,  1, -1, 1]]' &

# Changing a sign on row 1:
python M7k3.py '[[-1, -1, -1, -1, +1,  0, 0], 
                 [-1, -1,  0,  0,  1, -1, 1]]' &

# # Changing a sign on row 2:
# python M7k3.py '[[-1, -1, -1, -1, -1,  0, 0], 
#                  [-1, -1,  0,  0, -1, -1, 1]]' &
# 
# # Changing a sign on col 5:
# python M7k3.py '[[-1, -1, -1, -1, +1,  0, 0], 
#                  [-1, -1,  0,  0, -1, -1, 1]]' &
#
