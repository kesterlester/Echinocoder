from vertex_matches import (
    generate_canonical_vertex_matches,
    generate_all_vertex_matches,
    generate_all_vertex_match_signatures,
    _smallest_odd_number_greater_than_or_equal_to,
    generate_all_vertex_matches_given_equivalent_places,
    generate_all_vertex_matches_given_equivalent_places_IMPLEMENTATION_A,
    generate_all_vertex_matches_given_equivalent_places_IMPLEMENTATION_B,
    generate_viable_vertex_match_matrices,
)

from _vertex_matches import (
    generate_all_useful_vertex_matches,
    generate_all_vertex_matches_given_perming_places,
    generate_all_useful_vertex_matches_given_perming_places,
    _old_generate_all_vertex_match_signatures,
)

from equivalent_places import Equivalent_Places
from itertools import zip_longest

def test_versions():

    import numpy as np
    np_version_string = np.__version__
    np_ver_tup = tuple(int(x) for x in np_version_string.split("."))
    print(f"Numpy version {np_ver_tup}")

    assert np_ver_tup[0:3] >= (2,3,1), "Needed to fix https://github.com/numpy/numpy/pull/29223 and https://github.com/numpy/numpy/issues/28687 ."



# Vertex matches have an even number of +1 and and odd number of -1 entries, and others zero. Their total number of entries is M, the numnber of bad bats.
# The "signature" of a vertex match is how many ones, minus ones and zeros it contains.
# "Useful" vertex matches have at least k+1 non-zero entries (because all sums of <=k linearly dependent non-zero things in k-dimes are non-zero).
# A "canonical" vertex match is one where all the ones come before all the minus ones which come before all the zeros WITHIN any positions which are otherwise equivalent. . E.g., of all position are equivalent, then (1,1,-1,-1,-1,0) is a canonical match.

M0_all_vertex_matches_expected = [
    (),
]

M1_all_vertex_matches_expected = [
    (0,),
    (1,),
]

M2_all_vertex_matches_expected = [
    (0,0,),
    (0,1,),
    (1,0,),
    (1,1,),
]

M3_all_vertex_matches_expected = [
    (0, 0, 0),
    (0, 0, 1),
    (0, 1, 0),
    (1, 0, 0),
    (0, 1, 1),
    (1, 0, 1),
    (1, 1, 0),
    (1, 1, 1),
]
k5M8_all_vertex_match_signatures_expected = [
    # (0, 0, 8),
    # (1, 0, 7),
    # (2, 0, 6),
    # (3, 0, 5),
    # (4, 0, 4),
    # (5, 0, 3),
    (6, 0, 2),
    (7, 0, 1),
    (8, 0, 0),
]
k4M8_all_vertex_match_signatures_expected = [
    # (0, 0, 8),
    # (1, 0, 7),
    # (2, 0, 6),
    # (3, 0, 5),
    # (4, 0, 4),
    (5, 0, 3),
    (6, 0, 2),
    (7, 0, 1),
    (8, 0, 0),
]
k3M8_all_vertex_match_signatures_expected = [
    # (0, 0, 8),
    # (1, 0, 7),
    # (2, 0, 6),
    # (3, 0, 5),
    (4, 0, 4),
    (5, 0, 3),
    (6, 0, 2),
    (7, 0, 1),
    (8, 0, 0),
]
k2M8_all_vertex_match_signatures_expected = [
    # (0, 0, 8),
    # (1, 0, 7),
    # (2, 0, 6),
    (3, 0, 5),
    (4, 0, 4),
    (5, 0, 3),
    (6, 0, 2),
    (7, 0, 1),
    (8, 0, 0),
]
k5M7_all_vertex_match_signatures_expected = [
    # (0, 0, 7),
    # (1, 0, 6),
    # (2, 0, 5),
    # (3, 0, 4),
    # (4, 0, 3),
    # (5, 0, 2),
    (6, 0, 1),
    (7, 0, 0),
]
k4M7_all_vertex_match_signatures_expected = [
    # (0, 0, 7),
    # (1, 0, 6),
    # (2, 0, 5),
    # (3, 0, 4),
    # (4, 0, 3),
    (5, 0, 2),
    (6, 0, 1),
    (7, 0, 0),
]
k3M7_all_vertex_match_signatures_expected = [
    # (0, 0, 7),
    # (1, 0, 6),
    # (2, 0, 5),
    # (3, 0, 4),
    (4, 0, 3),
    (5, 0, 2),
    (6, 0, 1),
    (7, 0, 0),
]
k2M7_all_vertex_match_signatures_expected = [
    # (0, 0, 7),
    # (1, 0, 6),
    # (2, 0, 5),
    (3, 0, 4),
    (4, 0, 3),
    (5, 0, 2),
    (6, 0, 1),
    (7, 0, 0),
]
M3_all_vertex_match_signatures_expected = [
    (0, 0, 3),
    (1, 0, 2),
    (2, 0, 1),
    (3, 0, 0)
]

M4_all_vertex_matches_expected = [
     ( 0, 0, 0, 0,),

     (0, 0, 0, 1,),
     (0, 0, 1, 0,),
     (0, 1, 0, 0,),
     (1, 0, 0, 0,),

    (0, 0, 1, 1,),
    (0, 1, 0, 1,),
    (1, 0, 0, 1,),
    (0, 1, 1, 0,),
    (1, 0, 1, 0,),
    (1, 1, 0, 0,),

    (1, 1, 1, 0,),
    (1, 1, 0, 1,),
    (1, 0, 1, 1,),
    (0, 1, 1, 1,),

    (1, 1, 1, 1,),
    ]

M4_left_pair_right_pair = [
     #( 0, 0, 0,-1,),
     ( 0, 0,-1, 0,),#
     #( 0,-1, 0, 0,),
     (-1, 0, 0, 0,),#

     #( 0, 1, 1,-1,),
     ( 0, 1,-1 ,1,),#
     #( 0,-1, 1, 1,),
     #( 1, 0, 1,-1,),
     #( 1, 0,-1 ,1,),
     (-1, 0, 1, 1,),
     #( 1, 1, 0,-1,),
     #( 1,-1, 0, 1,),
     (-1, 1, 0, 1,),#
     ( 1, 1,-1, 0,),#
     #( 1,-1 ,1, 0,),
     #(-1, 1, 1, 0,),

     #( 0,-1,-1,-1,),
     (-1, 0,-1,-1,),#
     #(-1,-1, 0,-1,),
     (-1,-1,-1, 0,),#
    ]

k2M4_all_useful_vertex_match_signatures_expected = [
    #(0, 0, 4,),
    #(1, 0, 3,),
    #(2, 0, 2,),
    (3, 0, 1,),
    (4, 0, 0,),
    ]
M4_all_vertex_match_signatures_expected = [
    (0, 0, 4,),
    (1, 0, 3,),
    (2, 0, 2,),
    (3, 0, 1,),
    (4, 0, 0,),
]

k2M3_all_useful_vertex_matches_expected = [
    #(0, 0, 0),
    #(0, 0, 1),
    #(0, 1, 0),
    #(1, 0, 0),
    #(0, 1, 1),
    #(1, 0, 1),
    #(1, 1, 0),
    (1, 1, 1),
]
k2M4_all_useful_vertex_matches_expected = [
    (1,1,1,0,),
    (1,1,0,1,),
    (1,0,1,1,),
    (0,1,1,1,),

    (1,1,1,1,),
    ]
k2M4_all_useful_vertex_matches_with_1_perming_places = [
     #( 0, 0, 0,-1,),
     #( 0, 0,-1, 0,),
     #( 0,-1, 0, 0,),
     #(-1, 0, 0, 0,),

     #( 0, 1, 1,-1,),
     #( 0, 1,-1 ,1,),
     ( 0,-1, 1, 1,),
     #( 1, 0, 1,-1,),
     #( 1, 0,-1 ,1,),
     (-1, 0, 1, 1,),
     #( 1, 1, 0,-1,),
     ( 1,-1, 0, 1,),
     #(-1, 1, 0, 1,),
     #( 1, 1,-1, 0,),
     #( 1,-1 ,1, 0,),
     #(-1, 1, 1, 0,),

     ( 0,-1,-1,-1,),
     #(-1, 0,-1,-1,),
     #(-1,-1, 0,-1,),
     (-1,-1,-1, 0,),
]
k2M4_all_useful_vertex_matches_with_2_perming_places = [
     #( 0, 0, 0,-1,),
     #( 0, 0,-1, 0,),
     #( 0,-1, 0, 0,),
     #(-1, 0, 0, 0,),

     #( 0, 1, 1,-1,),
     ( 0, 1,-1 ,1,),
     ( 0,-1, 1, 1,),
     #( 1, 0, 1,-1,),
     ( 1, 0,-1 ,1,),
     (-1, 0, 1, 1,),
     #( 1, 1, 0,-1,),
     ( 1,-1, 0, 1,),
     (-1, 1, 0, 1,),
     ( 1, 1,-1, 0,),
     #( 1,-1 ,1, 0,),
     #(-1, 1, 1, 0,),

     ( 0,-1,-1,-1,),
     (-1, 0,-1,-1,),
     #(-1,-1, 0,-1,),
     (-1,-1,-1, 0,),
]
k3M4_all_useful_vertex_matches_expected = [
     (1, 1, 1, 1,),
    ]

"""
All USEFUL matches in k=2 dimensions, given M=4 bad bats, for perming_places=1
   1:    (0, -1, -1, -1)
   2:    (-1, -1, -1, 0)
   3:    (0, 1, 1, -1)
   4:    (-1, 1, 1, 0)
   5:    (1, 1, -1, 0)

All USEFUL matches in k=2 dimensions, given M=4 bad bats, for perming_places=2
   1:    (-1, 0, -1, -1)
   2:    (0, -1, -1, -1)
   3:    (-1, -1, -1, 0)
   4:    (-1, 0, 1, 1)
   5:    (0, -1, 1, 1)
   6:    (0, 1, 1, -1)
   7:    (1, 0, 1, -1)
   8:    (-1, 1, 1, 0)
   9:    (1, -1, 1, 0)
   10:    (1, 1, -1, 0)

"""

def test_helper_functions():

    for x,y in [
            (-0.5, 1),
            (0, 1),
            (1, 1),
            (1.1, 3),
            (2, 3),
            (3, 3),
            (3.0001, 5),
            (3.2, 5),
            (4.2, 5),
            (4.9999, 5),
            (5, 5),
            (5.0001, 7),
            ]:
        z = _smallest_odd_number_greater_than_or_equal_to(x)
        print(f"We hope that the smallest odd number greater than or equal to {x} is {y} and is also {z}")
        assert z==y
 
def test_signatures():
  for method in (
                generate_all_vertex_match_signatures,
                #_old_generate_all_vertex_match_signatures,
                ):
    test_programme = [
        (3, None, M3_all_vertex_match_signatures_expected, "M3 signatures"),
        (4, None, M4_all_vertex_match_signatures_expected, "M4 signatures"),
        (4, 2, k2M4_all_useful_vertex_match_signatures_expected, "k2M4 useful signatures"),

        (7, 2, k2M7_all_vertex_match_signatures_expected, "k2M7 useful signatures"),
        (7, 3, k3M7_all_vertex_match_signatures_expected, "k3M7 useful signatures"),
        (7, 4, k4M7_all_vertex_match_signatures_expected, "k4M7 useful signatures"),
        (7, 5, k5M7_all_vertex_match_signatures_expected, "k5M7 useful signatures"),

        (8, 2, k2M8_all_vertex_match_signatures_expected, "k2M8 useful signatures"),
        (8, 3, k3M8_all_vertex_match_signatures_expected, "k3M8 useful signatures"),
        (8, 4, k4M8_all_vertex_match_signatures_expected, "k4M8 useful signatures"),
        (8, 5, k5M8_all_vertex_match_signatures_expected, "k5M8 useful signatures"),
        ]
    for M, k, expected_signature, name in test_programme:
        LHS = sorted(list(method(M=M, k=k)))
        RHS = sorted(expected_signature)

        print(f"Test '{name}' has LHS and RHS:")
        print("[")
        for i,j in zip_longest(LHS,RHS):
            print(f"({i}, {j}), ")
        print("]")

        print()

        assert LHS==RHS

def test_main_generators():

    test_programme = [
        (None, 0, generate_all_vertex_matches(M=0), M0_all_vertex_matches_expected, "M0 all",),
        (None, 1, generate_all_vertex_matches(M=1), M1_all_vertex_matches_expected, "M1 all",),
        (None, 2, generate_all_vertex_matches(M=2), M2_all_vertex_matches_expected, "M2 all",),
        (None, 3, generate_all_vertex_matches(M=3), M3_all_vertex_matches_expected, "M3 all",),
        (None, 4, generate_all_vertex_matches(M=4), M4_all_vertex_matches_expected, "M4 all",),
        (2, 3, generate_all_useful_vertex_matches(k=2, M=3), k2M3_all_useful_vertex_matches_expected, "k2M3 useful",),
        (2, 4, generate_all_useful_vertex_matches(k=2, M=4), k2M4_all_useful_vertex_matches_expected, "k2M4 useful",),
        (3, 4, generate_all_useful_vertex_matches(k=3, M=4), k3M4_all_useful_vertex_matches_expected, "k3M4 useful",),

        (2, 4, generate_all_useful_vertex_matches_given_perming_places(k=2, M=4,perming_places=4), k2M4_all_useful_vertex_matches_expected, "k2M4 useful but testing perming_places=4}",),
        (2, 4, generate_all_useful_vertex_matches_given_perming_places(k=2, M=4, perming_places=4), list(generate_all_useful_vertex_matches(M=4, k=2, permute=True)), "k2M4 useful but testing perming_places=4}",), # permute=True is like perming_places=M
        (2, 4, generate_all_useful_vertex_matches_given_perming_places(k=2, M=4, perming_places=3), list(generate_all_useful_vertex_matches(k=2, M=4, permute=True)), "k2M4 useful but testing perming_places=3}",), # permute=True is also like perming_places=M-1
        (2, 4, generate_all_useful_vertex_matches_given_perming_places(k=2, M=4, perming_places=2), k2M4_all_useful_vertex_matches_with_2_perming_places, "k2M4 useful but testing perming_places=2}",),
        (2, 4, generate_all_useful_vertex_matches_given_perming_places(k=2, M=4, perming_places=1), k2M4_all_useful_vertex_matches_with_1_perming_places, "k2M4 useful but testing perming_places=1}",),
        (2, 4, generate_all_useful_vertex_matches_given_perming_places(k=2, M=4, perming_places=0), list(generate_all_useful_vertex_matches(M=4, k=2, permute=False)), "k2M4 useful but testing perming_places=0}",),

        (None, 4, generate_all_vertex_matches_given_perming_places(M=4, perming_places=4), M4_all_vertex_matches_expected, "M4 but testing perming_places=4}",),
        (None, 4, generate_all_vertex_matches_given_perming_places(M=4, perming_places=4), list(generate_all_vertex_matches(M=4, permute=True)), "M4 but testing perming_places=4}",), # permute=True is like perming_places=M
        (None, 4, generate_all_vertex_matches_given_perming_places(M=4, perming_places=3), list(generate_all_vertex_matches(M=4, permute=True)), "M4 but testing perming_places=3}",), # permute=True is also like perming_places=M-1
        (None, 4, generate_all_vertex_matches_given_perming_places(M=4, perming_places=0), list(generate_all_vertex_matches(M=4, permute=False)), "M4 but testing perming_places=0}",),


        (None, 4, generate_all_vertex_matches_given_equivalent_places(equivalent_places = Equivalent_Places(size=4, all_equivalent=True) ), list(generate_all_vertex_matches(M=4, permute=False)), "M4 but testing equivalent_places ALL",),
        # which should be the same as
        (None, 4, generate_all_vertex_matches_given_equivalent_places(equivalent_places = Equivalent_Places(equivalents_with_singletons=( (0,1,2,3,),    )) ), list(generate_all_vertex_matches_given_perming_places(M=4, perming_places=0)), "M4 but using equivalent_places to regenerate perming places 0",),
        # which leads to this sequence:
        (None, 4, generate_all_vertex_matches_given_equivalent_places(equivalent_places = Equivalent_Places(equivalents_with_singletons=( (0,),(1,2,3,),    )) ), list(generate_all_vertex_matches_given_perming_places(M=4, perming_places=1)), "M4 but using equivalent_places to regenerate perming places 1",),
        (None, 4, generate_all_vertex_matches_given_equivalent_places(equivalent_places = Equivalent_Places(equivalents_with_singletons=( (0,),(1,),(2,3,),    )) ), list(generate_all_vertex_matches_given_perming_places(M=4, perming_places=2)), "M4 but using equivalent_places to regenerate perming places 2",),
        (None, 4, generate_all_vertex_matches_given_equivalent_places(equivalent_places = Equivalent_Places(equivalents_with_singletons=( (0,),(1,),(2,),(3,),    )) ), list(generate_all_vertex_matches_given_perming_places(M=4, perming_places=3)), "M4 but equivalent_places to regenerate perming places 3",),
        (None, 4, generate_all_vertex_matches_given_equivalent_places(equivalent_places = Equivalent_Places(equivalents_with_singletons=( (0,),(1,),(2,),(3,),    )) ), list(generate_all_vertex_matches_given_perming_places(M=4, perming_places=4)), "M4 but equivalent_places to regenerate perming places 4",),
        # which should be the same as:
        (None, 4, generate_all_vertex_matches_given_equivalent_places(equivalent_places = Equivalent_Places(size=4, none_equivalent=True) ), list(generate_all_vertex_matches(M=4, permute=True)), "M4 but testing equivalent_places NONE",),

        (None, 4, generate_all_vertex_matches_given_equivalent_places(equivalent_places = Equivalent_Places(equivalents_with_singletons=( (0,1,),(2,3,),    )) ), M4_left_pair_right_pair, "M4 left pair right pair",),

        (None, 10,
     generate_all_vertex_matches_given_equivalent_places_IMPLEMENTATION_A(equivalent_places = Equivalent_Places(equivalents_with_singletons=( (0,7,),(1,4,5,6),(2,3,9,), (8,),  ))), 
     generate_all_vertex_matches_given_equivalent_places_IMPLEMENTATION_B(equivalent_places = Equivalent_Places(equivalents_with_singletons=( (0,7,),(1,4,5,6),(2,3,9,), (8,),  ))), 
      "Big 10 test.",),

        ]

    for k, M, gen, expected, name in test_programme:
 

        gen = tuple(gen)
        print(f"Test '{name}' has LHS and RHS:")
        print(f"gen {gen}")
        print(f"expected {expected}")

        LHS = sorted(gen)
        RHS = sorted(expected)

        print("[")
        for idx, (i,j) in enumerate(zip_longest(LHS,RHS)):
            print(f"{idx+1}:  ({i} {'==' if i==j else '!='} {j}), ")
        print("]")

        print()

        assert LHS==RHS

def test_canonical_order():
    # Output shoud be in order when permute=False
    for gen in (
          generate_canonical_vertex_matches(M=4, ),
          generate_canonical_vertex_matches(M=6, ),
          generate_canonical_vertex_matches(M=10, ),

          generate_canonical_vertex_matches(M=4, k=2),
          generate_canonical_vertex_matches(M=6, k=2),

          generate_canonical_vertex_matches(M=4, k=3),
          generate_canonical_vertex_matches(M=6, k=3),
          generate_canonical_vertex_matches(M=7, k=3),

          generate_canonical_vertex_matches(M=4, k=3, start=(0,0,1,1)),
          generate_canonical_vertex_matches(M=6, k=3, start=(0,0,0,0,1,1)),
          generate_canonical_vertex_matches(M=7, k=3, start=(0,1,1,1,1,1,1)),
          ):
        print("==================")
        print("New test")
        for vm in gen:
            print(f"We hope that this tuple is in non-decreasing order: {vm}")
            assert tuple(sorted(vm)) == vm
    

def test_canonical_order_stream():
   
    for M in range(1,15):
        for k in range(0,5):
          print(f"\nNew CANONICAL vertex match ORDER Test for M={M} and k={k}.")
          last_match = None
          for match in generate_canonical_vertex_matches(M=M, k=k):
            print(f"  order test saw match={match}")
            if last_match != None:
                # Check_order!
                assert last_match < match
            last_match = match

def test_order_stream():
    
    for gen in (
          generate_all_vertex_matches(M=4),
          generate_all_vertex_matches(M=6),
          # WILL NOTE WORK: generate_all_vertex_matches_given_equivalent_places_IMPLEMENTATION_A(equivalent_places = Equivalent_Places(equivalents_with_singletons=( (0,7,),(1,4,5,6),(2,3,9,), (8,),  ))),
          # WILL NOTE WORK: generate_all_vertex_matches_given_equivalent_places_IMPLEMENTATION_B(equivalent_places = Equivalent_Places(equivalents_with_singletons=( (0,7,),(1,4,5,6),(2,3,9,), (8,),  ))),
          ):
        print(f"\n New Order Test")
        last_signature, last_vertex_match = None, None
        for vertex_match in gen:
            signature = vertex_match.count(-1), vertex_match.count(0), vertex_match.count(+1)
            print(f"order test saw signature={signature} and vm={vertex_match}")
            # Only check_order progression when within a given signature due to perms:
            if last_signature != None and signature == last_signature:
                # Check_order!
                assert last_vertex_match < vertex_match
            last_signature, last_vertex_match = signature, vertex_match


def test_start_vertex_match_signatures():
  for method_with_start, method_without_start in (
       (generate_all_vertex_match_signatures, generate_all_vertex_match_signatures), 
       (_old_generate_all_vertex_match_signatures, _old_generate_all_vertex_match_signatures), 
    ):
    
     test_programme = [
       (10,2,(4,3,3)),
       (20,2,(2,13,5)),
       (20,2,(18,1,1)),
     ]
    
     from itertools import chain
     for M, k, start in test_programme:
         print(f"================== Testing startup of vertex SIGNATURE generators with M={M}, k={k}, start={start} =======")
         start_pos = None
         for i, (lesters, itertoolss) in enumerate(zip_longest(method_with_start(M=M, k=k), method_without_start(M=M, k=k))):
             print(f"{i}:         {lesters} {'==' if lesters==itertoolss else '!='} {itertoolss}")
             if lesters == start and start_pos == None:
                start_pos = i
             assert lesters==itertoolss
         print("Match confirmed!")
         assert start_pos is not None
         print(f"Start pos determined to be {start_pos}.")
         for i, (lesters, itertoolss) in enumerate(zip_longest(chain(iter([None,]*start_pos),method_with_start(M=M, k=k, start=start)), method_without_start(M=M, k=k))):
             if i < start_pos:
                 print(f"{i}:         {lesters} ...  {itertoolss}")
             else:
                 print(f"{i}:         {lesters} {'==' if lesters==itertoolss else '!='} {itertoolss}")
                 assert lesters==itertoolss

                
def test_start_vertex_matches():
    
     test_programme = [
       (3,1,(-1,-1,-1)),
       (6,3,(-1,0,1,1,1,1)),
       (7,2,(-1,-1,-1,0,0,0,0)),
       (7,2,(-1,0,0,0,0,1,1)),
       (7,2,(-1,1,1,1,1,1,1)),
     ]
    
     from itertools import chain
     for M, k, start in test_programme:
         for method_with_start, method_without_start in (
                 (generate_canonical_vertex_matches, generate_canonical_vertex_matches),
                 # TODO - UNCOMMENT WHEN START IMPLEMENTED FOR THESE! (generate_all_vertex_matches, generate_all_vertex_matches),
             ):
             print(f"================== Testing startup of vertex MATCH generators with M={M}, k={k}, start={start} =======")
             start_pos = None
             for i, (lesters, itertoolss) in enumerate(zip_longest(method_with_start(M=M, k=k), method_without_start(M=M, k=k))):
                 print(f"{i}:         {lesters} {'==' if lesters==itertoolss else '!='} {itertoolss}")
                 if lesters == start and start_pos == None:
                    start_pos = i
                 assert lesters==itertoolss
             print("Match confirmed!")
             assert start_pos is not None
             print(f"Start pos determined to be {start_pos}.")
             for i, (lesters, itertoolss) in enumerate(zip_longest(chain(iter([None,]*start_pos),method_with_start(M=M, k=k, start=start)), method_without_start(M=M, k=k))):
                 if i < start_pos:
                     print(f"{i}:         {lesters} ...  {itertoolss}")
                 else:
                     print(f"{i}:         {lesters} {'==' if lesters==itertoolss else '!='} {itertoolss}")
                     assert lesters==itertoolss

def test_start_vertex_matches_given_equivalent_places():
    
     test_programme = [
       (Equivalent_Places(size=3, all_equivalent=True), None, (-1,-1,-1)),
       (Equivalent_Places(size=3, all_equivalent=True), None, (-1,0,0)),
       (Equivalent_Places(size=3, all_equivalent=True), None, (-1,1,1)),

       (Equivalent_Places(size=3, all_equivalent=True), 1, (-1,-1,-1)),
       (Equivalent_Places(size=6, all_equivalent=True), 3, (-1,0,1,1,1,1)),
       (Equivalent_Places(size=7, all_equivalent=True), 2, (-1,-1,-1,0,0,0,0)),
       (Equivalent_Places(size=7, all_equivalent=True), 2, (-1,0,0,0,0,1,1)),
       (Equivalent_Places(size=7, all_equivalent=True), 2, (-1,1,1,1,1,1,1)),

       (Equivalent_Places(size=3, equivalents_without_singletons=((0,2),)),         1, (-1,-1,-1)),
       (Equivalent_Places(size=6, equivalents_without_singletons=((1,2), (3,4),)),  3, (-1,0,1,1,1,1)),
       (Equivalent_Places(size=7, equivalents_without_singletons=((2,4,6),(1,3),)), 2, (-1,-1,-1,0,0,0,0)),
       (Equivalent_Places(size=7, equivalents_without_singletons=((2,4,6),(1,3),)), 2, (-1,0,0,0,0,1,1)),
       (Equivalent_Places(size=7, equivalents_without_singletons=((2,4,6),(1,3),)), 2, (-1,1,1,1,1,1,1)),

       (Equivalent_Places(size=3, equivalents_without_singletons=((0,2),)),         None, (-1,-1,-1)),
       (Equivalent_Places(size=6, equivalents_without_singletons=((1,2), (3,4),)),  None, (-1,0,1,1,1,1)),
       (Equivalent_Places(size=7, equivalents_without_singletons=((2,4,6),(1,3),)), None, (-1,-1,-1,0,0,0,0)),
       (Equivalent_Places(size=7, equivalents_without_singletons=((2,4,6),(1,3),)), None, (-1,0,0,0,0,1,1)),
       (Equivalent_Places(size=7, equivalents_without_singletons=((2,4,6),(1,3),)), None, (-1,1,1,1,1,1,1)),

       (Equivalent_Places(size=3, none_equivalent=True), 1, (-1,-1,-1)),
       (Equivalent_Places(size=6, none_equivalent=True), 3, (-1,0,1,1,1,1)),
       (Equivalent_Places(size=7, none_equivalent=True), 2, (-1,-1,-1,0,0,0,0)),
       (Equivalent_Places(size=7, none_equivalent=True), 2, (-1,0,0,0,0,1,1)),
       (Equivalent_Places(size=7, none_equivalent=True), 2, (-1,1,1,1,1,1,1)),

       (Equivalent_Places(size=3, none_equivalent=True), None, (-1,-1,-1)),
       (Equivalent_Places(size=6, none_equivalent=True), None, (-1,0,1,1,1,1)),
       (Equivalent_Places(size=7, none_equivalent=True), None, (-1,-1,-1,0,0,0,0)),
       (Equivalent_Places(size=7, none_equivalent=True), None, (-1,0,0,0,0,1,1)),
       (Equivalent_Places(size=7, none_equivalent=True), None, (-1,1,1,1,1,1,1)),
     ]
    
     from itertools import chain
     for e_places, k, start in test_programme:
         for method_with_start, method_without_start in (
                 (generate_all_vertex_matches_given_equivalent_places, generate_all_vertex_matches_given_equivalent_places),
                 # TODO UNCOMMENT WHEN start IMPLEMENTED FOR _A: (generate_all_vertex_matches_given_equivalent_places_IMPLEMENTATION_A, generate_all_vertex_matches_given_equivalent_places),
                 (generate_all_vertex_matches_given_equivalent_places_IMPLEMENTATION_B, generate_all_vertex_matches_given_equivalent_places),
             ):
             print(f"================== Testing startup of vertex MATCH generators with e_places={e_places}, k={k}, start={start} =======")
             start_pos = None
             for i, (lesters, itertoolss) in enumerate(zip_longest(method_with_start(equivalent_places=e_places, k=k), method_without_start(equivalent_places=e_places, k=k))):
                 print(f"{i}:         {lesters} {'==' if lesters==itertoolss else '!='} {itertoolss}")
                 if lesters == start and start_pos == None:
                    start_pos = i
                 assert lesters==itertoolss
             print("Match confirmed!")
             assert start_pos is not None
             print(f"Start pos determined to be {start_pos}.")
             for i, (lesters, itertoolss) in enumerate(zip_longest(chain(iter([None,]*start_pos),method_with_start(equivalent_places=e_places, k=k, start=start)), method_without_start(equivalent_places=e_places, k=k))):
                 if i < start_pos:
                     print(f"{i}:         {lesters} ...  {itertoolss}")
                 else:
                     print(f"{i}:         {lesters} {'==' if lesters==itertoolss else '!='} {itertoolss}")
                     assert lesters==itertoolss

                    

def test_matrice_some_other_way():
    print("== Test of Matrix Generation =========")
    import sympy as sp
    def max_row_requirement(mat, max_rows):
        return sp.shape(mat)[0] <= max_rows

    def max_row_requirement(mat, max_rows):
        return sp.shape(mat)[0] <= max_rows

    for M,k in (
        (2,2),
        (3,2),
        (5,2),
        ):
        from functools import partial

        mat_gen_slow = generate_viable_vertex_match_matrices(
            M=M,
            k=k,
            old_yield_matrix = partial(max_row_requirement, max_rows=4),
            )

        mat_gen_fast = generate_viable_vertex_match_matrices(
            M=M,
            k=k,
            go_deeper = partial(max_row_requirement, max_rows=3),
            )

        print("=================================================")
        print(f"Will check if two methods agree for M={M}, k={k}:")
        print("Doing fast calc ...")
        fast = tuple(mat_gen_fast)
        print(f" ... len(fast)={len(fast)}.")
        print("Doing slow calc ...")
        slow = tuple(mat_gen_slow)
        print(f" ... len(slow)={len(slow)}.")
        assert fast == slow
        print("Fast agreed with slow")

        once = True
        for i, mat in enumerate(fast):
            if i<10 or i>len(fast)-10:
                print(i, mat)
            else:
                if once:
                    print(".....")
                    once=False
                continue
        print("===========================================")



