from numpy import array
from fractions import Fraction
from numbers import Number

assert isinstance(3, Number)
assert isinstance(3.5, Number)
assert isinstance(Fraction(3,4), Number)

def tuple_ize(a):
    """
    to_tuple turns numpy arrays (of any dimension or shape) into tuple arrays to make them hashable.

    Some examples of what should convert to what are given in the tuple_ize.unit_test_input_output_pairs attribute defined below.
    """

    if isinstance(a, Number):
        return a
    else:
        return tuple(map(tuple_ize, a))


"""
Use these below in unit tests if you want to check tuple_ize is working properly.
If all is working the first element of each pair shoudl convert to the second element of each pair.
"""
tuple_ize.unit_test_input_output_pairs = [
     (   array([1,2,5]),    (1,2,5)   ),
     (   array([[1,2],[5,6]]),    ((1,2),(5,6))   ),
     (   array([[1.2,2.2],[5.2,6.2]]),    ((1.2,2.2),(5.2,6.2))   ),
     (   array([[[1],[2]],[[5],[6]]]),    (((1,),(2,)),((5,),(6,)))   ),
     (   array([0,1,4])+Fraction(1),    (Fraction(1),Fraction(2),Fraction(5),)   ),
     (   array([0,1.5,4])+Fraction(1),    (Fraction(1),Fraction(5,2),Fraction(5),)   ),
     (   array([]),    tuple()   ),
     (   array([[],[]]),    ((),(),)   ),
     (   array([[[],],[[],]]),    (((),),((),))   ),
   ]
