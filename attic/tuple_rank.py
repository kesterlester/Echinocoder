#!/opt/local/bin/python3

# Christopher Lester

from math import comb
from itertools import pairwise

def tuple_rank(tup, k):
    """
    The purpose of this function is to map a supplied tuple of non-decreasing integers in [0,k-1] to its position in this list:

    (), (0), (1), (2), ... , (k-1), (0,0), (0,1), (0,k-1), (1,1), (1,2), ... , (1,k-1), (2,2), (2,3), ... , (2,k-1), (3,3), ... , (3,k-1), ... , (k-1,k-1), (0,0,0), (0,0,1), (0,0,2), ... , (0,0,k-1), (0,1,1), (0,1,2), ... , (1,1,1), (1,1,2), ... , (k-1,k-1,k-1), (0,0,0,0), ... , (k-1,k-1,k-1,k-1), ... .

    This list orders the tuples primarily by length, and secondarily by lexicographical order when tuples have the same length.

    Behaviour is formally undefined if the input is not a tuple of the type specified above. However, at time of writing the method will throw an exception if inputs seem to be inconsistent.

    Expected output for k=3:

    print(tuple_rank((), 3))        # Expected output: 0
    print(tuple_rank((0,), 3))      # Expected output: 1
    print(tuple_rank((1,), 3))      # Expected output: 2
    print(tuple_rank((2,), 3))      # Expected output: 3
    print(tuple_rank((0,0,), 3))    # Expected output: 4
    print(tuple_rank((0,1,), 3))    # Expected output: 5
    print(tuple_rank((0,2,), 3))    # Expected output: 6
    print(tuple_rank((1,1,), 3))    # Expected output: 7
    print(tuple_rank((1,2,), 3))    # Expected output: 8
    print(tuple_rank((2,2,), 3))    # Expected output: 9
    print(tuple_rank((0,0,0), 3))   # Expected output: 10

    There is a unit-test function which provides more test-cases.
    """

    assert isinstance(k,int), "k should be an integer."
    assert k>=0, "k should be a non-negative integer"
    assert all(isinstance(x, int) and x>=0 and x<k for x in tup), "values in tup should be integers in [0,"+str(k-1)+"]"
    assert len(tup)<=1 or min([b-a for a,b in list(pairwise(tup))])>=0, "tup should be non-decreasing"
    
    rank = 0
    n = len(tup)

    # Count all tuples of lengths less than the current tuple's length
    for length in range(n):
        rank += comb(k + length -1 , length)

    # Now, for the current tuple of length n, compute the lexicographical rank
    prev = 0
    for i in range(n):
        for val in range(prev, tup[i]):
            rank += comb(n - i + k - 1 - (val + 1), n - i - 1)
        prev = tup[i]

    return rank

def unit_test_tuple_rank():
    fails = 0

    fails += 0 != tuple_rank((), 1)
    fails += 0 != tuple_rank((), 10)
    fails += 0 != tuple_rank((), 100)
    fails += 0 != tuple_rank((), 1000)

    fails += 10 != tuple_rank((0,0,0), 3)
    fails += 15 != tuple_rank((0,0,0), 4)
    fails += 21 != tuple_rank((0,0,0), 5)

    fails += 0 != tuple_rank((),2)
    fails += 1 != tuple_rank((0,),2)
    fails += 2 != tuple_rank((1,),2)
    fails += 3 != tuple_rank((0,0,),2)
    fails += 4 != tuple_rank((0,1,),2)
    fails += 5 != tuple_rank((1,1,),2)
    fails += 6 != tuple_rank((0,0,0,),2)
    fails += 7 != tuple_rank((0,0,1,),2)
    fails += 8 != tuple_rank((0,1,1,),2)
    fails += 9 != tuple_rank((1,1,1,),2)
    fails += 10 != tuple_rank((0,0,0,0,),2)

    fails += 0 != tuple_rank((),3)
    fails += 1 != tuple_rank((0,),3)
    fails += 2 != tuple_rank((1,),3)
    fails += 3 != tuple_rank((2,),3)
    fails += 4 != tuple_rank((0,0,),3)
    fails += 5 != tuple_rank((0,1,),3)
    fails += 6 != tuple_rank((0,2,),3)
    fails += 7 != tuple_rank((1,1,),3)
    fails += 8 != tuple_rank((1,2,),3)
    fails += 9 != tuple_rank((2,2,),3)
    fails += 10 != tuple_rank((0,0,0,),3)
    fails += 11 != tuple_rank((0,0,1,),3)
    fails += 12 != tuple_rank((0,0,2,),3)
    fails += 13 != tuple_rank((0,1,1,),3)
    fails += 14 != tuple_rank((0,1,2,),3)
    fails += 15 != tuple_rank((0,2,2,),3)
    fails += 16 != tuple_rank((1,1,1,),3)
    fails += 17 != tuple_rank((1,1,2,),3)
    fails += 18 != tuple_rank((1,2,2,),3)
    fails += 19 != tuple_rank((2,2,2,),3)
    fails += 20 != tuple_rank((0,0,0,0,),3)
    fails += 21 != tuple_rank((0,0,0,1,),3)
    fails += 22 != tuple_rank((0,0,0,2,),3)
    fails += 23 != tuple_rank((0,0,1,1,),3)
    fails += 24 != tuple_rank((0,0,1,2,),3)
    fails += 25 != tuple_rank((0,0,2,2,),3)
    fails += 26 != tuple_rank((0,1,1,1,),3)
    fails += 27 != tuple_rank((0,1,1,2,),3)
    fails += 28 != tuple_rank((0,1,2,2,),3)
    fails += 29 != tuple_rank((0,2,2,2,),3)
    fails += 30 != tuple_rank((1,1,1,1,),3)
    fails += 31 != tuple_rank((1,1,1,2,),3)
    fails += 32 != tuple_rank((1,1,2,2,),3)
    fails += 33 != tuple_rank((1,2,2,2,),3)
    fails += 34 != tuple_rank((2,2,2,2,),3)
    fails += 35 != tuple_rank((0,0,0,0,0,),3)

    if fails>0:
        raise Exception("tuple_rank fails unit test")

    if fails==0:
        print("tuple_rank passes Unit Test")

if __name__ == "__main__":
    unit_test_tuple_rank()
