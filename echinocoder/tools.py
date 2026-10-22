import numpy as np
import inspect
from itertools import permutations, product
#from distinct_permutations_with_leftovers import distinct_permutations_with_leftovers as distinct_permutations # Use either this or the next line!
from more_itertools import distinct_permutations # Use either this or the last line!


class LineNo:
    def __str__(self):
        return str(inspect.currentframe().f_back.f_lineno)

__line__ = LineNo()

def numpy_array_of_frac_to_str(arr : np.array):
    """
    Numpy arrays of fractions print in a really ugly way. This method is designed to output a human-readable
    array representation that looks more readable. Its output is only for human consumption, not archival, so
    don't try to read it back in. If the array happens not to conain fractions, it should still cope ... it will just
    print the array as usual.

    Args:
        arr: array to print

    Returns:

    """

    # This next line is heuristic but good enough for now as it will distinguish float and integer from Fraction
    # (which is all we need it to do).
    # is_array_of_fractions = ( arr.dtype == "O" )

    # if not is_array_of_fractions:
    #  Fall back to normal str
    #   return str(arr)

    ans = "["
    for row in arr:
        ans += "["
        for elt in row:
            ans += str(elt) + ", "
        ans += "],  "
    ans += "]"
    return ans

def ascending_data(np_array_representing_set : np.ndarray) -> np.ndarray:
        # The input is assumed to be an np_array with n elements, each being a vector, and with each vector holding k components.
        # It represents a multiset of those n vectors. 
        # The following "ascending data" has the x-components in ascending order, the y-components in asceding order,
        # and so on. This has broken up the vectors.  I.e. the j=1 vector in ascending_data is not likely to
        # be any of the vectors in the input (unless the data was already sorted appropriately).
        # You can think of "ascending data" as representing all the things we want to encode EXCEPT the associations
        # which link every element of each vector up in the right way.
        return sort_each_np_array_column(np_array_representing_set)

def expand_complex_to_real_pairs(np_array):
    """
    This function expands numpy arrays of complex numbers into arrays which are twice as long but which have real parts separated from complex parts.
    E.g.:   [3+4i, 4+2i] -> [3, 4, 4, 2]
    Its inverse is "real_pairs_to_complex_zip".
    """

    a = np_array.real
    b = np_array.imag
    c = np.empty((a.size+b.size,), dtype=a.dtype)
    c[0::2] = a
    c[1::2] = b
    return c

def real_pairs_to_complex_zip(np_array):
    """
    This function zips adjacent reals in an array of length n into complex numbers in an array of length floor(n/2).
    E.g.:  [3, 4, 4, 2] -> [3+4i, 4+2i]
    or:    [3, 4, 4,  ] -> [3+4i, 4   ]
    Its inverse is "expand_complex_to_real_pairs".
    """
    reals = np_array[0::2]
    imags = np_array[1::2]
    if len(imags)<len(reals):
        imags = np.append(imags, [0])
    print("input ",np_array)
    print("reals ",reals)
    print("images ",imags)
    return reals+1j*imags

def invert_perm(perm):
    """given [3, 0, 1, 2] should generate [1, 2, 3, 0]"""
    return sorted(range(len(perm)), key=perm.__getitem__)

def sort_each_np_array_row(array: np.ndarray) -> np.ndarray:
    """
    Sort the elements of each row into numerical order.
    For example:
                  np.array([[1, 0, 2],
                            [0, 5, 2],
                            [3, 0, 8]])
    sorts to:
                  np.array([[0, 1, 2],
                            [0, 2, 5],
                            [0, 3, 8]])
    """
    return np.sort(array)

def sort_each_np_array_column(array: np.ndarray) -> np.ndarray:
    """
    Sort the elements of each column into numerical order.
    For example:
                  np.array([[1, 0, 3],
                            [0, 5, 2],
                            [3, 0, 8]])
    sorts to:     
                  np.array([[0, 0, 2],
                            [1, 0, 3],
                            [3, 5, 8]])
    """
    return np.sort(array.T).T

def sort_np_array_rows_lexicographically(array: np.ndarray) -> np.ndarray:
    """
    Permutes rows of a numpy array (individual rows are preserved) so that the rows end up in lexicographical order.
    E.g.:
                  np.array([[1, 0, 2],
                            [0, 5, 2],
                            [3, 0, 8]])
    sorts to:     
                  np.array([[0, 5, 2],
                            [1, 0, 2],
                            [3, 0, 8]])
    """
    return array[np.lexsort(array.T[::-1])]

def permute_columns_except_first(arr):
    # Ensure the input is a numpy array
    arr = np.array(arr)
    rows, cols = arr.shape

    # Separate the fixed first column from the other columns
    fixed_col = arr[:, [0]]
    variable_cols = arr[:, 1:]

    # Generate permutations for each column independently as generators
    # print("Making all perms ... ")
    all_perms = (distinct_permutations(variable_cols[:, i]) for i in range(variable_cols.shape[1]))
    # print("                 ... done.")

    # Create the product of all column permutations (as a generator)
    #print("Making all products ... ")
    for perm_set in product(*all_perms):
        #print("                         ... got a product")
        combined = np.hstack([fixed_col] + [np.array(col).reshape(-1, 1) for col in perm_set])
        yield combined

