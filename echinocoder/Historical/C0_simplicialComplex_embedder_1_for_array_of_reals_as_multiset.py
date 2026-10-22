#!/usr/bin/env python3

# USE WITH CAUTION!  No known bugs, but not tested to destruction.
# This is a re-implementation of C0_bug2_numpy_simplicialComplex_embedder_for_array_of_reals_as_multiset.py but aiming to fix the bug in the part of the code which mods out S(n).
#
# Patrick Kennedy-Hunt
# Christopher Lester
#
# Where used in this file the expression Deltak (or $\Delta^k$ in TeX)
# refers to the space inside a unit k-simplex.
# A point in Deltak could be parameterised by  k real numbers:
#
#      x_0, x_1, ... , x_(k-1)
#
# with the property that:
#
#      0 <= x_i <= 1 for all i
#
# and for which
#
#      sum_i x_i <= 1.
#
# The implementation below  also handles points that live inside spaces which are
# products of n copies of Deltak. These would be in $(\Delta^k)^n$ in TeX, but
# in the sourcecode we term this space Deltakn.
#
# For practical reasons, the implementation below mostly stores points in Deltak
# or in Deltakn not as lists, i.e. not as
#
#  x = [x_0, x_1, ... , x_(k-1)]
#
# but instead stores them as dictionaries mapping key index pairs (j,i) to the
# simplex coordinates. In these key index pairs (j,i) the value j in [0,n-1]
# indexes the copy of Deltak within Dektakn,  and the i value in [0,k-1]
# indexes the coordinate position within Deltak.
#
# For example:, for n=2 and k=2 ...
#
# a point x within the first Deltak within Deltakn might be stored as
#
#          x = { (0,0):0.25,  (0,1):0.50, (0,2):0.10, }
#
# and a point y within the second Deltak within Deltakn might be stored as
#
#          y = { (1,0):0.05, (1,1):0.00,  (1,2):0.90, }
#
# or equivalentlay (since zero coefficients may be omitted) as
#
#          y = { (1,0):0.05,              (1,2):0.90, }
#
# and a point delta = x*y  within Detlakn would be stored as
#
#          delta = { **x, **y }
#                = {
#                    (0,0):0.25,  (0,1):0.50, (0,2):0.10,
#                    (1,0):0.05,              (1,2):0.90,
#                  }.

from MultisetEmbedder import MultisetEmbedder
from sys import version_info
import numpy
import tools
from typing import Self

if not version_info >= (3, 7):
    assert False, "We need at least python 3.7 as we rely on dictionaries being ordered! See DIC_ORDER_USED"

from itertools import pairwise
import numpy as np
import unittest
from dataclasses import dataclass, field
from typing import Union, Any
import hashlib

class Embedder(MultisetEmbedder):
    def embed_generic(self, data: np.ndarray, debug=False) -> (np.ndarray, Any):
        metadata = None
        return self.embed_internal(data), metadata

    def embed_kOne(self, data: np.ndarray, debug=False) -> (np.ndarray, Any):
        metadata = None
        return MultisetEmbedder.embed_kOne_sorting(data), metadata
 
    def embed_internal(self, data: Union[np.ndarray, 'Position_within_Simplex_Product'],
               use_n2k2_optimisation=False, input_is_in_DeltakToN=False) -> np.ndarray:
        """
        By default, this function takes as input a length-n list of k-vectors (i.e. points in (R^k)^n)
        to embed in a perutation invariant way in R^(2nk+1).
    
        E.g. to embed the three two-vectors [1,2], [4,2] and [2,0] one might supply it with
    
            data = np.array([[1,2], [4,2], [2,0]]).
    
        However, if the flag "input_is_in_DkToN" is set to true, or if the input is
        detected as being of type "Position_within_Simplex_Product", the input is assumed to be already mapped into
        simplex coordinates for a space which is an n-fold product of k-simplices.  This feature is to allow
        simple unit tests.
    
        Args:
            data in  (R^k)^n (or in  (Delta^k)^n in some cases)
            use_n2k2_optimisation = a flag to request interpretation as (Delta^k)^n
    
        Returns:
            embedding in R^(nk+1)
        """
    
        n, k = data.shape # Valid whether input_is_in_DeltakToN is True or False
    
        if n*k > np.iinfo(Eji_LinComb.INT_TYPE).max:
            raise ValueError("Part of alg stores vertices in Eji_LinComb.INT_TYPE for has so can't work for nk above a threshold.")
    
        if type(data) == Position_within_Simplex_Product:
            delta = data # Don't reformat the data as the data is already in (Delta^k)^n.
        elif input_is_in_DeltakToN:
            # Reformat the data as supposedly the coefficients are correct but the type is not:
            delta = Position_within_Simplex_Product(data)
        else:
            # Embed each vector in a simplex:
            delta = vectors_to_delta(data) # convert (R^k)^n to (Delta^k)^n.
    
        if use_n2k2_optimisation and n == 2 and k == 2:
            c_dc_pairs = make_c_dc_pairs_n2k2(delta)
        else:
            c_dc_pairs = make_c_dc_pairs(delta)
    
        c_dc_pairs_sorted_by_dc = sorted(c_dc_pairs, key=lambda x: x[1], reverse=True) # largest dc first.
        #print("c_dc_pairs .... in order of decreasing dc.") # Only needed for debug.
        #[print(_) for _ in c_dc_pairs_sorted_by_dc] # Only needed for debug.
    
        #simplex = Maximal_Simplex([c for c, _ in c_dc_pairs])  # Only needed for debug.
        # print("simplex (before modding by S(n)) =")  # Only needed for debug.
        # print(simplex) # Only needed for debug.
    
        # Barycentric-subdivide the mother simplex "simplex" into smaller ones, of which the "daugter_simplex"
        # is the one containing our important point.
        len_c_dc_pairs_sorted_by_dc = len(c_dc_pairs_sorted_by_dc)
    
        daughter_simplex_vertices_with_dc = [
            (Eji_LinComb(n,k,[c for c,_ in c_dc_pairs_sorted_by_dc[:i+1]]),
             (c_dc_pairs_sorted_by_dc[i][1] - (c_dc_pairs_sorted_by_dc[i+1][1] if i+1<len_c_dc_pairs_sorted_by_dc else 0)) # The dc values for the daughter simplex are DIFFERENCES of dc values in the mother simplex.
         * (i+1)    # They are also weighted by their index (1,2,3,4,etc) to make them identically distributed.
             ) for i in range(len_c_dc_pairs_sorted_by_dc)
        ]
        #print("\nDaughter simplex vertices before S(n)")
        #[print(ejilc, dc) for ejilc,dc in daughter_simplex_vertices_with_dc]
    
        daughter_simplex_vertices_with_dc_after_Sn = [
            (ejilc.get_canonical_form(), dc) for ejilc,dc in daughter_simplex_vertices_with_dc
        ]
    
        #print("\nDaughter simplex vertices after S(n)")
        #[print(ejilc, dc) for ejilc, dc in daughter_simplex_vertices_with_dc_after_Sn]
    
        big_n_complex = 2 * n * k + 1
        point_in_bigNComplex = sum([dc * pr(ejilc.hash_to_unit_complex_number(), big_n_complex) for ejilc, dc in daughter_simplex_vertices_with_dc_after_Sn]) + np.zeros(big_n_complex)  # Addition of zero
    
        real_array = point_in_bigNComplex.real
        imag_array = point_in_bigNComplex.imag
        assert imag_array[0] == 0 # Since pr's first elt is 1 which has 0 imaginary part
        big_n_real = 4 * n * k + 1 # Note: not 2*(2*n*k+1) because we don't record the always-0 above.
    
        point_in_R_bigNReal = np.concatenate((real_array, imag_array[1:])) # Omitting the first elt of imag_array as it is always zero
    
        assert len(point_in_R_bigNReal) == self.size_from_n_k(n,k)
        assert big_n_real == self.size_from_n_k(n,k)
    
        return point_in_R_bigNReal
    
        # #[print(vertex) for vertex in simplex.get_vertex_list()] # Only needed for debug.
        # #print("simplex_eji_ordering (before mod S(n)) = ") # Only needed for debug.
        # #[print(eji) for eji in simplex.get_Eji_ordering()] # Only needed for debug.
        #
        # #print("simplex_eji_ordering (after mod S(n)) = ") # Only needed for debug.
        # #[print(eji) for eji in simplex.get_Eji_ordering().get_canonical_form()] # Only needed for debug.
        #
        # #print("simplex (after modding by S(n)) =") # Only needed for debug.
        # #print(simplex.get_canonical_form()) # Only needed for debug.
        # #[print(vertex) for vertex in simplex.get_canonical_form().get_vertex_list()] # Only needed for debug.
        #
        # ####TEST_REMOVE#### # Now 'canonicalise' the vertices in c_dc_pairs using that perm:
        # ####TEST_REMOVE#### c_dc_pairs_after_mod_Sn = [ ({ (inverse_perm[j], i) for (j,i) in c }, dc)  for (c,dc) in c_dc_pairs   ]
        #
        # # TEST DON'T APPLY PERM!!!!! It was a mistake!!!!!
        # c_dc_pairs_after_mod_Sn = c_dc_pairs
        #
        # # print("c_dc_pairs (after modding by S(n)) =")
        # # [ print(c) for c in c_dc_pairs_after_mod_Sn ]
        #
        # #shrink = True
        # #c_l_dc_triples = [(c, ell(c, k, shrink=shrink), dc) for (c, dc) in c_dc_pairs_after_mod_Sn]
        # # print("c_l_dc_triples (after modding by S(n)) =")
        # # [ print(c) for c in c_l_dc_triples ]
        #
        # # Please someone re-implement this dot product without so many comprehensions ... or at any rate BETTER:
        # # Want output here to be sum_i pr(r_i, big_n) x_i)
        # # where, in effect, r_i and x_i would be defined by
        # # [ blah for _, r_i, x_i in c_l_dc_triples ]
        #
        # # print("pr(20)=",pr(20,big_n))
        # point_in_R_bigN = sum([d * pr(ell, big_n) for _, ell, d in c_l_dc_triples]) + np.zeros(big_n)  # Addition of zero
        # # term at end ensures that we still get a zero vec (not 0) in the event that c_l_dc_triples is empty!
        #
        # return point_in_R_bigN
    
    def size_from_n_k_generic(self, n: int, k: int) -> int:
        return 4*n*k + 1
    

@dataclass
class Position_within_Simplex:
    """
    A point inside a k-simplex is defined by k reals values delta_i with the property that:

    * 0 <= delta_i for all i in 0,1,2,...,k-1, and
    * sum_{i=0}^{k-1} delta_i <= 1.

    The purpose of this class is to hold and manipulate such delta_i, and to manage conversions
    between such things and other coordinate systems.
    """
    _np_ar : np.array

    def __getitem__(self, item):
        return self._np_ar[item]

    def __init__(self, x):
        self._np_ar = np.array(x)

    def __eq__(self, other):
        return (self._np_ar == other._np_ar).all()

    def check_valid(self):
        #print("  XOO", self._np_ar)
        assert sum(self._np_ar) <= 1, "Total of simplex coordinates should not exceed 1."
        for coeff in self._np_ar:
            #print("    YOO", coeff)
            assert coeff >= 0, "Every simplex coordinate should be non-negative."

@dataclass
class Position_within_Simplex_Product:
    """
    Conceptually this holds "n" copies of Point_in_simplex.
    However, for performance reasons we hold it as a numpy array, with
    first index j in [0,n-1] identifies the simplex, and the
    second index i in [0,k-1] identifies the position within the j-th simplex.
    """
    _np_ar : np.array
    shape : tuple

    def __init__(self, x):
        if isinstance(x, np.ndarray):
            self._np_ar = x
        elif isinstance(x, list) and x and isinstance(x[0], Position_within_Simplex):
            self._np_ar = np.array([pos._np_ar for pos in x])
        else:
            # Hope for best:
            self._np_ar = np.array(x)
        # In all cases:
        self._set_shape()

    def _set_shape(self):
        self.shape = np.shape(self._np_ar)

    def __getitem__(self, item):
        """
        If called with one index, get Position_within_Simplex.
        If called with two indices, get coeff directly from within Position_within_Simplex.
        """
        if isinstance(item, tuple):
            # Get coeff directly from within Position_within_Simplex.
            return self._np_ar[item]

        # Get Position_within_Simplex.
        return Position_within_Simplex(self._np_ar[item])

    def check_valid(self):
        for pos in self._np_ar:
            Position_within_Simplex(pos).check_valid()

from collections import namedtuple
Eji = namedtuple("Eji", ["j", "i"])

@dataclass
class Eji_Ordering:
    """Class to hold eij orderings (biggest first)."""
    _eji_list: list[Eji]

    def __len__(self):
        return len(self._eji_list)

    def __iter__(self):
        return iter(self._eji_list)

    def get_j_order(self, from_right=False):
        """
        Returns a list of j-values in the order in which they first appear.
        By default, we read from beginning of list. Otherwise (if from_right=True) we read from the other end.
        E.g.
                    Eji(2,_), Eji(5,_), Eji(2,_), Eji(9,_), ...
        would generate
                    [2, 5, 9, ...
        since 2 is the first j-value to appear, 5 the next, and 9 is the third (appearing after an ignored j=2 repeat).
        """
        # Do we want 0,1,2,3,4,5 (starting from left) or 5,4,3,2,1,0 (starting from right) ? Select with from_right.

        # Note that setting from_right does not (in general) reverse the answer even though it reverses the input.
        # I.e. perm_from_right(ordering)[::-1] is not in general the same as perm_from_left(ordering[::-1]).

        if from_right:
            return list({j[0]: None for j in self._eji_list[::-1]}) # Uses insertion order preservation! DIC_ORDER_USED
        else:
            # Default!
            return list({j[0]: None for j in self._eji_list      }) # Uses insertion order preservation! DIC_ORDER_USED

    def get_perm(self):
        """
        This returns the perm of [0,1,2,..,nk-1] which would take the observed j-order to
        the canonical order [nk-1,nk-2,...,2,1,0].
        This order is a bit silly (it's backwards). But this is a throwback to early hand calculations and
        a convention that they established (that the largest Eji comes first).  One could just as easily have used
        the opposite order, although a bunch of other things would have to change to compensate.
        TODO: Consider reimplementing with the opposite convention!
        """
        return tools.invert_perm(self.get_j_order(from_right=True))

    def get_canonical_form(self):
        """Returns the ordering in canonical form."""
        # This implementation gets a j-ordering, treats it as a perm, then inverts that perm, then applies it.
        # That seems rather laborious. There is probably a better way.
        # TODO: Find a better way!
        perm = self.get_perm()
        return Eji_Ordering([Eji(perm[j], i) for (j, i) in self._eji_list])

    def check_valid(self):
        j_vals = [e.j for e in self._eji_list]
        i_vals = [e.i for e in self._eji_list]

        presumed_n = max(j_vals, default=-1)+1
        presumed_k = max(i_vals, default=-1)+1
        assert presumed_n*presumed_k == len(self._eji_list), "Expect n*k entries in the ordering!"
        assert set(self._eji_list) == {Eji(j, i)
                                       for j in range(presumed_n)
                                       for i in range(presumed_k)}, "Each Eji value shjould appear once!"

@dataclass
class Maximal_Simplex_Vertex:
    _vertex_set: set[Eji] = field(default_factory=set)

    def __len__(self) -> int:
        return len(self._vertex_set)

    def __iter__(self):
        return iter(self._vertex_set)

    def get_canonical_form(self):
        """Mod out by Sn for this single vertex, ignoring any others."""
        # Method: sort the Eji by the i index, then populate the j's in order.
        sorted_eji_list = sorted(list(self._vertex_set), key=lambda eji: eji.i)
        renumbered_eji_list = [ Eji(j=j, i=eji.i) for j,eji in enumerate(sorted_eji_list)]
        return Maximal_Simplex_Vertex(set(renumbered_eji_list))

    def check_valid(self):
        # every j index in the set must appear at most once
        j_vals = { eji.j for eji in self._vertex_set }
        assert len(j_vals) == len(self._vertex_set)

    def get_permuted_by(self, perm):
        return Maximal_Simplex_Vertex({Eji(perm[eji.j], eji.i) for eji in self._vertex_set})

@dataclass
class Eji_LinComb:
    INT_TYPE = np.uint16 # uint16 should be enough as the eij_counts will not exceed n*k which can therefore reach 65535

    _index : INT_TYPE
    _eji_counts : np.ndarray

    def index(self) -> INT_TYPE:
        """How many things were added together to make this Linear Combination."""
        return self._index

    def hash_to_128_bit_md5_int(self):
        """
        N.B. This hash is of self._eji_counts only -- i.e. it ignores self._index.
        For the purposes to which this hash will be used, that is believed to be apporopriate.
        """
        m = hashlib.md5(self._eji_counts)
        return int.from_bytes(m.digest(), 'big') # 128 bits worth.

    def hash_to_64_bit_reals_in_unit_interval(self):
        """
        An md5 sum is 64 bits long so we get two such reals.
        N.B. This hash is of self._eji_counts only -- i.e. it ignores self._index.
        For the purposes to which this hash will be used, that is believed to be apporopriate.
        """

        x = self.hash_to_128_bit_md5_int()
        bot_64_bits = x & 0xffFFffFFffFFffFF
        top_64_bits = x >> 64
        return np.float64(top_64_bits)/(1 << 64), np.float64(bot_64_bits)/(1 << 64)

    def hash_to_unit_complex_number(self):
        mu,_ = self.hash_to_64_bit_reals_in_unit_interval()
        theta = mu*2*numpy.pi
        return complex(numpy.cos(theta), numpy.sin(theta))

    def __init__(self, n: int, k: int, list_of_Maximal_Simplex_Vertices: list[Maximal_Simplex_Vertex] | None = None):
        self._index = Eji_LinComb.INT_TYPE(0)
        self._eji_counts = np.zeros((n, k), dtype=Eji_LinComb.INT_TYPE, order='C')
        if list_of_Maximal_Simplex_Vertices:
            for msv in list_of_Maximal_Simplex_Vertices: self.add(msv)

    def _setup_debug(self, index: int, eji_counts: np.ndarray): # Really just for unit tests. Don't use in main alg code.
        self._index = Eji_LinComb.INT_TYPE(index)
        self._eji_counts = np.asarray(eji_counts, dtype=Eji_LinComb.INT_TYPE, order='C')

    def add(self, msv: Maximal_Simplex_Vertex):
        self._index += 1
        for j, i in msv: self._eji_counts[j, i] += 1


    def __eq__(self, other: Self):
        return self._index == other._index and np.array_equal(self._eji_counts, other._eji_counts)

    def __ne__(self, other: Self):
        return not self.__eq__(other)

    def get_canonical_form(self) -> Self:
        ans = Eji_LinComb.__new__(Eji_LinComb)
        ans._index = self._index
        ans._eji_counts = tools.sort_np_array_rows_lexicographically(self._eji_counts)
        return ans

@dataclass
class Maximal_Simplex:
    """These are stored in a list which is required to be ordered from big to small
    under the same ordering used to order eji's."""
    _vertex_list: list[Maximal_Simplex_Vertex]
    _eji_ordering_cache: Union["Eji_Ordering", "None"] = field(compare=False) # don't compare caches on equality test
    _canonical_simplex: Union["Maximal_Simplex", "None"] = field(compare=False) # don't compare caches on equality test

    def __init__(self, list_of_vertices):
        self._vertex_list = list_of_vertices
        self._eji_ordering_cache = None
        self._canonical_simplex = None

    def get_canonical_form(self):
        # If already calculated, return the canonical simplex from cache:
        if self._canonical_simplex is not None:
            return self._canonical_simplex

        # Need to initialise canonical simplex cache:
        perm = self.get_Eji_ordering().get_perm()
        self._canonical_simplex = Maximal_Simplex([vertex.get_permuted_by(perm) for vertex in self._vertex_list])
        return self._canonical_simplex

    def get_vertex_list(self) -> list[Maximal_Simplex_Vertex]:
        return self._vertex_list

    def get_Eji_ordering(self) -> Eji_Ordering:
        # If already calculated, return the eji ordering from cache:
        if self._eji_ordering_cache is not None:
            return self._eji_ordering_cache

        # Need to initialise eji ordering cache:
        vertices_and_null = self._vertex_list + [Maximal_Simplex_Vertex()]
        self._eji_ordering_cache = Eji_Ordering([(v1._vertex_set - v2._vertex_set).pop() for v1, v2 in pairwise(
            vertices_and_null)]) # The set c1-c2 should contain only one element, so pop() should return it.
        return self._eji_ordering_cache

    def check_valid(self):
        """This function asserts all sorts of things that should be true for a valid object.
        Use it in unit-tests to raise assertion errors."""
        for v_larger, v_smaller in pairwise(self._vertex_list):
            v_larger.check_valid()
            v_smaller.check_valid()
            #print("LARGER SMALLER ",v_larger,v_smaller)
            gain = v_larger._vertex_set - v_smaller._vertex_set
            loss = v_smaller._vertex_set - v_larger._vertex_set
            assert len(gain) == 1, "Every vertex should have an Eji not contained in the next vertex."
            if loss:  # Something got deleted, so check that what was deleted grew in i
                assert len(loss) == 1
                for eji_loss in loss:
                    break
                for eji_gain in gain:
                    break
                assert eji_loss.j == eji_gain.j
                assert eji_loss.i + 1 == eji_gain.i
            if gain and not loss: # new j got added, so check that what was deleted grew started at base
                for eji_gain in gain:
                    break
                assert eji_gain.i == 0

#@dataclass
#class Maximal_Simplex:
#    """Class to hold any of the big simplices which (before barycentric subdivision)
#    form the beginnings of our simplicial complex.  It holds an eij_ordering."""
#    eji_ordering: Eji_Ordering
#    vertices: Maximal_Simplex_Vertices
#
#    def __init__(self, ):
#        pass

def ell(c, k, shrink=False):
    """
    The function ell(c, k) can (in principle) be anything that satisfies the properties listed below.
    In practice, presumably some choices are better than others for efficiency reasons or reasons of practicality.
    TODO: Therefore one should experiment with alternatives that go beyond the choices made here!

    INPUT REQUIREMENTS:
    * The input k shall be a non-negative integer
    * The input c shall be a (possibly empty) set of (j,i) pairs, with j in [0,n-1] and i in [0,k-1].
    * The set of j values in c shall have no repeates.  I.e. if there are three pairs in c then one pair will have j=0, one pair will have j=1 and one pair will have j=2.  [Aside: this restricts the length of c to be at most n.]

    What properties should ell satisfy? 
   
    * The input c={(0,4),(2,1)} "represents" the sum of basis vectors e^0_4 + e^2_1 which is clearly the same as e^2_1 + e^0_4 as addition is commutative. This is why c is stored as a set.
    * ell should map every possible c to a natural number.
    * at fixed k, it must be the case that ell(c1,k)==ell(c2,k) <==> c1==c2
    
    E.g. we must have all these equal:
    
        ell( {(1,5),(2,42),(0,100)}, 101)
        ell( {(2,42),(1,5),(0,100)}, 101)

    yet these must all differ:

        ell( {(1,5),(2,42),(0,100)}, 101)     (start)
        ell( {(2,5),(1,42),(0,100)}, 101)     (sets differ in content)
        ell( {(1,5),(2,43),(0,100)}, 101)     (43!=42)
        ell( {(1,5),       (0,100)}, 101)     (sets differ in length)
    """
    
    assert isinstance(k,int), "k should be an integer."
    assert k >= 0, "k should be a non-negative integer."
    len_c = len(c)
    assert len({ j for j,_ in c }) == len_c, "No j should repeat."

    #print("c was = ",c)

    i_vals_sorted_by_j_vals = [i for _,i in sorted(list(c)) ] # Fortunately j comes first in (j,i) so this will work.
    # Note that above we do not bother storing the j values as we know they will be 0,1,..,len(c)-1 given
    # the input preconditions.  Thus we will be done if we can hash i_vals_sorted_by_j_vals to a non-colliding value.
    # print("i_vals_sorted_by_j_vals = ",i_vals_sorted_by_j_vals)

    # Either of the two lines below should work identically. Just use one!
    # ell_value = sum(( i*(k**(len_c-1-pos))  for pos,i in enumerate(i_vals_sorted_by_j_vals) ))  +  sum((k**i for i in range(len_c)))
    ell_value_integer = sum(( i*(k**(len_c-1-pos)) + k**pos for pos,i in enumerate(i_vals_sorted_by_j_vals) ))

    # print("ell_value before nonlinear tfm= ",ell_value_integer)
    
    # At this point ell is an integer in 0,1,2,3,4,5, ...
    # We may wish to non-linearly transform it into (-pi/2,pi/2)
   
    if shrink:
        alternate_sign = 2*(ell_value_integer % 2)-1 # Will be +1 or -1 depending on even/oddness of ell_value_ingeger
        ell_value = ell_value_integer * alternate_sign
        ell_value = np.arctan(np.float64(ell_value)/50) # puts ell in (-pi/2,pi/2)
    else:
        ell_value = ell_value_integer

    return ell_value

def make_flat_sums(delta, sort=False, prepend_zero=False):
    #print("MFS",delta)
    n, k = np.shape(delta)
    """
    Conceptually, the dict flat_sums (for k=3) should "represent" the following map:
    
    flat_sums = {
      {(0,0)+(0,1)+(0,2)} : delta[0,0]+delta[0,1]+delta[0,2],
      {      (0,1)+(0,2)} :            delta[0,1]+delta[0,2],
      {            (0,2)} :                       delta[0,2],
      {(1,0)+(1,1)+(1,2)} : delta[1,0]+delta[1,1]+delta[1,2],
      {      (1,1)+(1,2)} :            delta[1,1]+delta[1,2],
      {            (1,2)} :                       delta[1,2],
      ...
      ...
      ...
      {(n-1,0)+(n-1,1)+(n-1,2)} : delta[n-1,0]+delta[n-1,1]+delta[n-1,2],
      {        (n-1,1)+(n-1,2)} :              delta[n-1,1]+delta[n-1,2],
      {                (n-1,2)} :                           delta[n-1,2],
    }.

    However, for practical implementation purposes we "represent" the above by the following more abbreviated structure
    (a list of tuples)
    in which for each key the first index is j in [0,n-1], the second index is
    a tuple containing the i values in the basis sum,
    I.e. key (j,r) below embeds key sum([ (j,i) for i in r]) above. The example below uses k=3 again
    
    flat_sums = [
      (0, (0,1,2,), delta[0,0]+delta[0,1]+delta[0,2]),
      (0,   (1,2,),            delta[0,1]+delta[0,2]),
      (0,     (2,),                       delta[0,2]),
      (1, (0,1,2,), delta[1,0]+delta[1,1]+delta[1,2]),
      (1,   (1,2,),            delta[1,1]+delta[1,2]),
      (1,     (2,),                       delta[1,2]),
      ...
      ...
      ...
      (n-1, (0,1,2,), delta[n-1,0]+delta[n-1,1]+delta[n-1,2]),
      (n-1,   (1,2,),              delta[n-1,1]+delta[n-1,2]),
      (n-1,     (2,),                           delta[n-1,2]),
    ]
    """

    flat_sums=[ (j, tuple(range(i_min, k)), sum([delta[j,i] for i in range(i_min,k) ])) for j in range(n) for i_min in range(k) ]
    # print("flat_sums unsorted",flat_sums)
    if sort:
        flat_sums = sorted(flat_sums, key=lambda x : (x[2], len(x[1])) ) # Sort by delta sum, but break ties in favour of longer sums
        # print("flat_sums sorted",flat_sums)
    if prepend_zero:
        flat_sums = [ ( None, tuple(), 0) ] + flat_sums
    return flat_sums

def make_c_dc_pairs_n2k2(delta: Position_within_Simplex_Product):
    a=delta[0,0]
    b=delta[0,1]
    c=delta[1,0]
    d=delta[1,1]

    if c+d >=d >= a+b >= b: # C'
        c_dc_pairs = [
            (Maximal_Simplex_Vertex({(0,1), (1,1)}),       b),
            (Maximal_Simplex_Vertex({(0,0), (1,1)}),       a),
            (Maximal_Simplex_Vertex({       (1,1)}),   d-a-b),
            (Maximal_Simplex_Vertex({       (1,0)}),       c),
            ]
    elif c+d >= a+b >= b >= d: # A'
        c_dc_pairs = [
            (Maximal_Simplex_Vertex({(0,1), (1,1)}),       d),
            (Maximal_Simplex_Vertex({(0,1), (1,0)}),     b-d),
            (Maximal_Simplex_Vertex({(0,0), (1,0)}),       a),
            (Maximal_Simplex_Vertex({       (1,0)}), c+d-a-b),
            ]
    elif c+d >= a+b >= d >= b: # B'
        c_dc_pairs = [
            (Maximal_Simplex_Vertex({(0,1), (1,1)}),       b),
            (Maximal_Simplex_Vertex({(0,0), (1,1)}),     d-b),
            (Maximal_Simplex_Vertex({(0,0), (1,0)}),   a+b-d),
            (Maximal_Simplex_Vertex({       (1,0)}), c+d-a-b),
            ]
    elif a+b >= c+d >= d >= b: # A
        c_dc_pairs = [
            (Maximal_Simplex_Vertex({(0,1), (1,1)}),       b),
            (Maximal_Simplex_Vertex({(0,0), (1,1)}),     d-b),
            (Maximal_Simplex_Vertex({(0,0), (1,0)}),       c),
            (Maximal_Simplex_Vertex({(0,0)}),        a-c+b-d),
            ]
    elif a+b >= c+d >= b >= d: # B
        c_dc_pairs = [
            (Maximal_Simplex_Vertex({(0,1), (1,1)}),       d),
            (Maximal_Simplex_Vertex({(0,1), (1,0)}),     b-d),
            (Maximal_Simplex_Vertex({(0,0), (1,0)}),   c+d-b),
            (Maximal_Simplex_Vertex({(0,0)}),        a-c+b-d),
            ]
    elif a+b >= b >= c+d >= d: # C
        c_dc_pairs = [
            (Maximal_Simplex_Vertex({(0,1), (1,1)}),       d),
            (Maximal_Simplex_Vertex({(0,1), (1,0)}),       c),
            (Maximal_Simplex_Vertex({(0,1)}),          b-c-d),
            (Maximal_Simplex_Vertex({(0,0)}),              a),
            ]
    else:
        assert False

    return c_dc_pairs

def make_c_dc_pairs(delta : Position_within_Simplex_Product):
    """
    We output a list of c_dc_pairs.  Each pair (c,dc) contains a Maximal_Simplex_Vertex, c, together with a
    coefficient, dc. The simplex vertex, c, is coded as a list of (j,i) values representing e^j_i, i.e. the ith basis
    vector of the j-th simplex. E.g. if the simplex vertex is c=[(0,1),(2,2)] then c represents e^0_1+e^2_2.
    The coefficient dc, attached to c, says how much of c is needed to represent the component of delta in that
    direction.

    Pruning entries with dc=0 might seem like a good idea.  However:

       (1) pruning zeros may also destroy regularity/predictability;  e.g. people might prefer to see c_dc_pairs always have the same length as the number of non-origin simplex points.
       (2) pruning zeros requires a test for equality on a floating point number, which is a bit silly.  Default should therefore be NOT pruning zeros,
       (3) worst of all, pruning zeors makes it impossible to identify and canonicalise the simplex vertices.

    Hence, we do not prune zeros!
    """

    n,k = np.shape(delta)

    flat_sums = make_flat_sums(delta, sort=True)

    dc_vals = [ sum2-sum1 for (j1,i1_vals, sum1),(j2,i2_vals, sum2) in pairwise([(None, tuple(), 0), ]+flat_sums) ]

    c_dc_pairs = [(Maximal_Simplex_Vertex({Eji(j, max(moo)) for j in range(n) if
                    (moo := [min(iis) for (jj, iis, _) in flat_sums[index:] if jj == j])}), dc_vals[index]) for index in
                    range(len(flat_sums))  #if not prune_zeros or dc_vals[index] != 0
                    ]  # See set note below

    """A set rather than a list is used to hold the coordinate vectors because later we want to find out 
    "elements in one set not in another" ... and so if we had used lists we would have to construct a set from a 
    list later anyway.  Fortunately the objects represented are sets anyway (they represent sums of dissimilar
     basis elements which are vertices, and sums are order independent).  The set creation comprehension does not 
     produce duplicate elements squashed by the set, though, so if it's later needed they could be changed back to a 
     list here (instead of set) so long as the later set-difference calculation is done some other way."""

    return c_dc_pairs

def pr(r, big_n):
    return np.power(r, np.arange(big_n)) # Starting at zeroeth power so that r can be both zero and non-zero without constraint.

def vector_to_simplex_point(vec):
    k = len(vec)
    vec = np.array(vec)
    return 1.0/(k*(1.0+np.power(2.0,vec))) # TODO: This somewhat crude parameterisation does not use the WHOLE
    # of the simplex -- so it's a bit wasteful. It also has terrible dynamic range problems and even unit issues.
    # Might want to address all of these points with a better mapping.

def vectors_to_delta(vecs):
    n=len(vecs)
    if len(vecs)==0:
        return Position_within_Simplex_Product([list() for i in range(n)])
    # vecs is not empty
    k = len(vecs[0])
    delta = np.zeros((n,k),)
    for j in range(n):
        vec = vecs[j]
        simplex_point = vector_to_simplex_point(vec)
        k_this = len(vec)
        if k!=k_this:
            raise Exception("Vectors supplied to vectors_to_delta are not all the same dimension!")
        for i in range(k):
            delta[j,i]=simplex_point[i]
    return Position_within_Simplex_Product(delta)

def test_simplex_embedding():

    import test_PKH_alg
    suite = unittest.TestLoader().loadTestsFromModule(test_PKH_alg)
    unittest.TextTestRunner(verbosity=2).run(suite)

    """
    short = embed

    ans1 = short( delta = {  (0,1) : 0.5, (1,2) : 0.25 }, input_is_in_DeltakToN)

    # Next three all similar to each other.
    ans2a = short(n=3, k=3, 
                  delta = {  (0,2) : 0.5, (1,2) : 0.25, (2,2):0.1 }, )
    ans2b = short(n=3, k=3, 
                  delta = {  (0,2) : 0.5, (1,2) : 0.25, (2,2):0.1,  (1,1):0.0001}, )
    ans2c = short(n=3, k=3, 
                  delta = {  (0,2) : 0.5, (1,2) : 0.25, (2,2):0.1001,  }, )

    # Perm invariance
    ans3c1 = short(n=3, k=3, 
                   delta = {  (0,2) : 0.5, (1,0):0.001, (1,2) : 0.25, (2,2):0.1001,  }, )
    ans3c2 = short(n=3, k=3, 
                   delta = {  (1,2) : 0.5, (0,0):0.001, (0,2) : 0.25, (2,2):0.1001,  }, )

    # Trick case:
    ans4 = short(n=3, k=3, 
             delta = {
             (0,0) : 0.5, (0,1) : 0.2, (0,2) : 0.1,    #a point in the 1st simplex
             (1,2) : 0.25,                             #a point in the 2nd simplex 
             (2,0) : 0.1,                              #a point in the 3rd simplex
             })

    # Trick case:
    ans5 = short(n=7, k=3, 
             delta = {
             (0,0) : 0.5, (0,1) : 0.2, (0,2) : 0.1,    #a point in the 1st simplex
             (1,2) : 0.25,                             #a point in the 2nd simplex 
             (2,0) : 0.1,                              #a point in the 3rd simplex
             })

    # Zero cases:
    ans6a = short(n=7, k=3, 
                  delta = { (0,1) : 0.0, (1,2) : 0.0, })
    ans6b = short(n=7, k=3,
                  delta = dict()                       )
    ans6c = short(n=7, k=3,
                  delta = vectors_to_delta( [] )       )

    ans7 = short(n=3, k=2,
                  delta = vectors_to_delta( [
                     np.array([1,2]), # k-vector 1 of n
                     np.array([1,0]), # k-vector 2 of n
                     np.array([5,2]), # k-vector 3 of n
                  ] )
                )

    enc1 = embed(np.array([[1,2],[1,0],[5,2]]))
    """

    # n=4
    # k=3
    # delta=dict()
    # delta[(0,0)]=0.00
    # delta[(0,1)]=0.03
    # delta[(0,2)]=0.22
    # delta[(1,0)]=0.00
    # delta[(1,1)]=0.11
    # delta[(1,2)]=0.10
    # delta[(2,0)]=0.10
    # delta[(2,1)]=0.00
    # delta[(2,2)]=0.50
    # delta[(3,0)]=0.01
    # delta[(3,1)]=0.03
    # delta[(3,2)]=0.20

    delta = Position_within_Simplex_Product([[0.00, 0.03, 0.22],
                                             [0.00, 0.11, 0.10],
                                             [0.10, 0.00, 0.50],
                                             [0.01, 0.03, 0.20]])

    ans8 = embed(delta)

    """
    print("Ans1 was ",ans1)
    print("Ans2a was ",ans2a)
    print("Ans2b was ",ans2b)
    print("Ans2c was ",ans2c)
    print("Ans3c1 was ",ans3c1)
    print("Ans3c2 was ",ans3c2)
    print("Ans5 was ",ans5)
    print("Ans6a was ",ans6a)
    print("Ans6b was ",ans6b)
    print("Ans6c was ",ans6c)
    print("Ans7 was ",ans7)
    print("Ans8 was ",ans8)

    print("enc1 was ",enc1)
    """

    print()
    print()
    print("Ans8 was ",ans8)

if __name__ == "__main__":
    test_simplex_embedding()
       #unittest.main(exit=False)
