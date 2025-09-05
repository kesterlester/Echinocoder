import sympy as sp
import math
from more_itertools import distinct_permutations
from distinct_partitions_with_start import distinct_partitions_with_start as distinct_partitions
from bi_range import bi_range_with_maxes
from equivalent_places import Equivalent_Places
from functools import partial
import sympy_tools

"""
Vertex matches have an even number of +1 and and odd number of -1 entries, and others zero. Their total number of entries is M, the numnber of bad bats.

"Useful" vertex matches have at least k+1 non-zero entries (because all sums of <=k linearly dependent non-zero things in k-dimes are non-zero).

A "canonical" vertex match is one where the elements in the tuple never decrease reading left to right. E.g., if all position are equivalent, then (1,1,-1,-1,-1,0) is a canonical match.

Sometimes it is not worth permuting vertex matches over every bad bat because other matches in the existing context may not yet have broken any symmetries between the bats.

E.g. if the first two places are different to each other (and to any other place) but the last four places are all equivalent to each other, then one need only consider these orderings of the letters "Speedo" in that context:

    Sp.eedo
    Se.pedo
    Sd.peeo
    So.peed

    pS.eedo
    pe.Sedo
    pd.Seeo
    po.Seed

    eS.pedo
    ep.pedo
    ee.pedo  # Yes, this group has five on account of the double "e" !
    ed.pedo
    eo.Sped

    dS.peeo
    dp.Seep
    de.Speo
    do.Spee

    oS.peed
    op.Seed
    oe.Sped
    od.Spee
"""

def _smallest_odd_number_greater_than_or_equal_to(x):
    return 2*math.ceil((x+1)/2)-1 

def generate_all_vertex_match_signatures(
    M, #number of bad bats
    k = None, # k=dimension of space (supply k if you want to calculate only useful matches, otherwise omit)
    start = None,
    ):
    """
    The signature of a vertex match is how many ones, minus ones and zeros it has.
    We yield triplets of numbers having the order (number_ones, number_minus_ones, number_zeros).

    If suplied and not None, start must be a tuple containing a signature which the method would ordinarily generate, and as a consequence the generator will start here.

    Vertex matches have M entries in total, comprising an even number of +1 and and odd number of -1 entries, and others zero.

    "Useful" vertex matches have at least k+1 non-zero entries (because all sums of <=k linearly dependent non-zero things in k-dimes are non-zero).

    We choose to generate the signatures in the order such that a canonical vertex match (i.e. a tuple like (-1,0,0,0,1,1,1,1)
    in which the elements are non-decreasing) which represents each signature would come out in ascending tuple order.

    IMPLEMENTATION:

    All signatures have at least one "minus 1". We will break it out and add it back in only at the very end (see (*)).
    In addition, there is a 'free for all' region of the signature where minus ones come in paris, ones come in pairs and zeros come in pairs.
    This region has a total of 2*((M-1)//2) digits, with a totally unconstrained composition (other than that all are in pairs).
    Finally, there is also an extra packing zero when M is even ... however we don't really need to think much about this as 
    we know that e+o+z = M so we can just calculate z from z = M - (e+o) at the end.

    The next few lines of code implements the above.

    NOTATION:

    We abbreviate the number_of_ones as "e" as there are an Even number of them.
    We abbreviate the number_of_minus_ones as "o" as there are an Odd number of them.
    We abbreviate the number_of_zeros as "z".
    """

    # if False: # Uncomment to use old method for testing purposes
    #    import _vertex_matches
    #    yield from _vertex_matches._old_generate_all_vertex_match_signatures(M=M, k=k, start=start)
    #    return

    total_among_pairs = 2*((M-1)//2)

    if start is not None:
        if len(start) != 3:
            raise ValueError(f"len(start) should equal 3 but is {len(start)}. start={start}.")
        if sum(start) != M:
            raise ValueError(f"sum(start) should equal {M} but is {sum(start)}. start={start}.")
        if True in ((int(c) != c or c<0) for c in start):
            raise ValueError(f"start should be a tuple of non-negative integers but start={start}.")

        start_e, start_o, start_z = start

        # Work backwards to start of e_among_pairs and start of z_among_pairs by reverse engineering yield line:
        start_e_among_pairs = start_e
        start_o_among_pairs = start_o - 1 # See (*) above and (*) below.
        start_z_among_pairs = total_among_pairs - start_e_among_pairs - start_o_among_pairs
        starting = True
    else:
        starting = False

    if k is None:
        # This is the k->minus_infinity limit of the "k is not None" case below.
        for o_among_pairs in range(start_o_among_pairs if starting else total_among_pairs, -1, -2):
            for z_among_pairs in range(start_z_among_pairs if starting else total_among_pairs - o_among_pairs, -1, -2):
                starting = False
                e_among_pairs = total_among_pairs - o_among_pairs - z_among_pairs
                e, o = e_among_pairs, o_among_pairs + 1 # (*) There is always one extra o.
                z = M - (e+o)
                yield e,o,z
    else:
        if k % 2:
            k = k+1 # This correction ensures z_among_pairs starts at an even number (and in fact the correct number!).
        for o_among_pairs in range(start_o_among_pairs if starting else total_among_pairs, -1, -2):
            for z_among_pairs in range(start_z_among_pairs if starting else min(total_among_pairs - o_among_pairs, total_among_pairs - k), -1, -2):
                starting = False
                e_among_pairs = total_among_pairs - o_among_pairs - z_among_pairs
                e, o = e_among_pairs, o_among_pairs + 1 # (*) There is always one extra o.
                z = M - (e+o)
                yield e,o,z

def generate_canonical_vertex_matches(
        M, # M=number of bad bats
        k=None, # k=dimension of space (supply k if you want to calculate only useful matches, otherwise omit)
        start = None,
        ):
        """
        M should be a non-negative integer. It specifies how long each generated tuple will
        be. (The number of bad bats!)

        k, if not None, should be the dimension of space.  Supply k if you want to generate
        only the useful matches for that spatial dimension.

        This method generates all tuples with the following properties:

           * the tuple has length M,
           * the tuple has an even number of ones,
           * the tuple has an odd number of minus ones,
           * the tuple is composed only of ones, minus ones and zeros,
           * if k is not None, then the tuple shall have AT LEAST k+1 non-zero entries, and
           * the tuple is in canonical form (i.e. sorted into non-decreasing order).
        
        If start is supplied (which must be something which the stream would ordinarily output) 
        the generator should start from there instead of starting at the beginnning.
        """
        return generate_all_vertex_matches(M=M, k=k, permute=False, start=start)

def generate_all_vertex_matches(
        M, # M=number of bad bats
        k=None, # k=dimension of space (supply k if you want to calculate only useful matches, otherwise omit)
        permute = True,
        start = None,
        ):
        """
        M should be a non-negative integer. It specifies how long each generated tuple wil
        be. (The number of bad bats!)

        k, if not None, should be the dimension of space.  Supply k if you want to generate only
        the useful matches for that spatial dimension.

        This method generates all tuples with the following properties:

           * the tuple has length M,
           * the tuple has an even number of ones,
           * the tuple has an odd number of minus ones,
           * the tuple is composed only of ones, minus ones and zeros, and
           * if k is not None, then the tuple shall have AT LEAST k+1 non-zero entries.
        
        By default, each tuple is yielded in every possible distinct ordering of its elements.
        However, this perming can be disabled by setting permute=False. This will result in each
        tuple being yielded once only in a canonical form (ie. all minus ones followed by all zeros
        followed by all ones). Note canonical form is also sorted into non-decreasing order!

        If start is supplied (which must be something which the stream would ordinarily output) 
        the generator should start from there instead of starting at the beginnning.
        At present non-none start is only implemented for permute=False, so you will get an exception
        if you try to use it with permute=True.
        """

        if permute and start is not None:
            raise NotImplementedError("Sorry, we don't yet implement start when permute=True")

        if start is not None:
            start_signature = start.count(1), start.count(-1), start.count(0)

        for number_of_ones, number_of_minus_ones, number_of_zeros in generate_all_vertex_match_signatures(M, k=k, start=start_signature if (start is not None) else None):

            ones = (1,)*number_of_ones
            minus_ones = (-1,)*number_of_minus_ones
            zeros = (0,)*number_of_zeros

            tup = minus_ones + zeros + ones # Note numerical order!

            if permute:
                if start is not None:
                    raise NotImplementedError("Sorry, we don't yet implement start when permute=True")
                for match in distinct_permutations(tup):
                    yield match
            else:
                yield tup

def generate_all_useful_vertex_matches(
    M, # M=number of bad bats
    k, # k=dimension of space 
    permute = True,
    start = None,
    ):
    return generate_all_vertex_matches(M=M, k=k, permute=permute, start=start)

def generate_all_vertex_matches_given_equivalent_places_IMPLEMENTATION_A(
    equivalent_places : Equivalent_Places,
    # M=None, # M = number of bad bats. (Can be derived from equivalent_places, so no longer supplied)
    k=None, # k=dimension of space (supply k if you want to calculate only useful matches, otherwise omit)
    start=None,
    ):

    if start is not None:
        raise NotImplementedError

    M = equivalent_places.size
    if int(M) != M or M<0:
        raise RuntimeError("Equivalent_Places is not behaving!")

    if M==0:
        return

    assert M>0

    ##############################################
    def _generate_dicts_for(e_places, signature):
        # Caller must guarantee e_places is non-empty, M>0 and signature consistent with e_places
        tot = sum(signature)
        assert e_places
        number_of_ones, number_of_minus_ones, number_of_zeros = signature
        assert tot == sum(len(e_place) for e_place in e_places)
    
        perming_places = len(e_places[0])
        non_perming_places = tot - perming_places
    
        # Recursion will stop when non_perming_places reaches 0 .... which should be the same as wheb e_places is length 1.
        # Let's check the above statement.
        assert (non_perming_places > 0 and len(e_places) > 1) or (non_perming_places == 0 and len(e_places)==1)
    
        for perming_ones, non_perming_ones in bi_range_with_maxes(number_of_ones, max_first=perming_places, max_second=non_perming_places):
            for perming_minus_ones, non_perming_minus_ones in bi_range_with_maxes(number_of_minus_ones, max_first = perming_places-perming_ones, max_second=non_perming_places - non_perming_ones):
                perming_zeros = perming_places - (perming_ones + perming_minus_ones)
                non_perming_zeros = number_of_zeros - perming_zeros
    
    
                assert perming_zeros >=0
                assert non_perming_zeros >=0
    
                perming_part = (-1,)*perming_minus_ones + (0,)*perming_zeros + (1,)*perming_ones # Note numerical order!
                assert len(perming_part) == len(e_places[0])
    
                perming_part_dict = dict(zip(e_places[0], perming_part))
    
                if non_perming_places == 0:
                    yield perming_part_dict
                else:
                    new_signature = (non_perming_ones, non_perming_minus_ones, non_perming_zeros)
                    for non_perming_part_dict in _generate_dicts_for(e_places[1:], new_signature):
                       yield perming_part_dict | non_perming_part_dict
    ##############################################

    e_places = equivalent_places._equivalent_places_with_singletons
    for signature in generate_all_vertex_match_signatures(M,k=k):

        for d in _generate_dicts_for(e_places, signature):
            yield tuple(d[i] for i in range(M))

def generate_all_vertex_matches_given_equivalent_places_IMPLEMENTATION_B(
        equivalent_places : Equivalent_Places,
        # M=None, # M = number of bad bats. (Can be derived from equivalent_places, so no longer supplied)
        k=None, # k=dimension of space (supply k if you want to calculate only useful matches, otherwise omit)
        start=None,
        ):

    M = equivalent_places.size
    if int(M) != M or M<0:
        raise RuntimeError("Equivalent_Places is not behaving!")

    if M==0:
        return

    assert M>0

    e_places = equivalent_places._equivalent_places_with_singletons
    splitting = tuple( len(group) for group in e_places )

    starting = start is not None

    # As vertex_matches (below) are canonical, they are sorted, so:
    start_vertex_match = tuple(sorted(start)) if starting else None

    start_partition = None
    if starting:
        #raise NotImplementedError()

        """
        The intial canonical vertex_match might be, say

              start_vertex_match=(-1,0,0,1,1)

        and if the equivalent places were 

              e_places=((1,3,4),(0,),(2,))

        then we would have

              splitting = (3,1,1)

        and so could concevably have encountered the partition

               partition = ( (-1,0,1), (1,), (0,) )

        which would have pushed

            (-1,0,1) into pos(1,3,4)
            (1,) into pos(0,) and
            (0,) into pos(2,)

        resulting in the following yield:

            y = (1,-1,0,0,1).

        When constructing the starting_partition our job would be:

            GIVEN y, e_places, splitting and start_vertex_match, FIND partition.

        How can we do this?
        """
        start_partition = tuple( tuple( start[pos] for pos in pos_group  )  for pos_group in e_places ) 

    # Now start the actual iteration:

    workspace = [None]*M

    for vertex_match in generate_canonical_vertex_matches(M=M, k=k, start=start_vertex_match):
        for partition in distinct_partitions(vertex_match, splitting, start=start_partition if starting else None):
            starting = False # This is important! We only start once.

            # That's all the looking done. We now just need to fill in the workspace ....
            assert len(partition) == len(e_places)
            for payload_group, pos_group in zip(partition, e_places):
                assert len(payload_group) == len(pos_group)
                for payload, pos in zip(payload_group, pos_group):
                    assert 0 <= pos < M
                    assert payload in (+1, -1, 0)
                    workspace[pos] = payload
            # .... OK the workspace is now filled so:

            yield tuple(workspace)

def generate_all_vertex_matches_given_equivalent_places(
        equivalent_places : Equivalent_Places,
        # M=None, # M = number of bad bats. (Can be derived from equivalent_places, so no longer supplied)
        k=None, # k=dimension of space (supply k if you want to calculate only useful matches, otherwise omit)
        start=None
        ):
    #return generate_all_vertex_matches_given_equivalent_places_IMPLEMENTATION_A(equivalent_places, k=k, start=start)
    return generate_all_vertex_matches_given_equivalent_places_IMPLEMENTATION_B(equivalent_places, k=k, start=start)

def generate_viable_vertex_match_matrices(
    M, # M = number of bad bats. 
    k, # k=dimension of space.
    remove_obvious_collapses = True, # Discards matrices whose RRE form have a row with betwen 1 and k non-zero elements. (Nb: ihis setting forces rre to be calcualted.)
    debug_test_max_rows = True,
    return_mat = False,
    return_rre = False,
    return_rre_pivots = False,
    return_hashable_rre = False,
    return_confusable_sets = False,
    remove_duplicates_via_hash = False, # Making this true could crash your program via memory usage. Beware!  This setting forces rre to be calculated -- so no harm in also choosing to return it.
    go_deeper    = None, # If present, then the branch topped by matrix "mat" is only explored more deeply if go_deeper(mat) is True. Does not affect whether mat itself is yielded.
    old_yield_matrix = None, # If present, then the matrix "mat" is only yielded if old_yield_matrix(mat) is True.  If not yielded, further branch exploration is suppressed. Note that, other things being equal, and if it is physically possible, it is better to use "go_deeper" (with or without old_yield_matrix) than "old_yield_matrix" alone.
    confusable_sets_or_None_function = None, # Like old_yield_matrix, except that confusable_sets_or_None(mat) is assumed to return either None or a  tuple containing a two confusable multisets. And if there ARE confusable multisets (rather than None) it is assumed that the user would like these to be yielded. So none None her acts like old_yield_matrix == True.
    debug = False,
    ):
    """
    Generate sympy.Matrix objects which represent constraints on lattice alignments of red/blue vertices. 
    Uses depth-first traversal of rows.
    
    To abort the current branch after yielding a given matrix, the user may send True to the generator.

    Alteratively, the user may specify which matrices and branches should be explored by supplying one or both of the yield_matrix and go_deeper arguments.

    It is far better to kill a branch before generating its daughter matrixes than to kill a branch by killing/vetoing each daughter matrix.  This, if it is possible to do so, it is far better to use "go_deeper" (with or without  "yield_matrix") to kill a whole branch in one test, than to use only "yield_matrix".
    """

    if return_confusable_sets and confusable_sets_or_None_function == None:
        raise ValueError("You cannot ask to return confusable sets unless you also supply a confusable_sets_or_None_function.")

    max_rows = sympy_tools.max_rows_for_viable_stripped_RRE_matrix(M=M, k=k)

    calculate_hashable_rre_early = remove_duplicates_via_hash
    calculate_hashable_rre_late = return_hashable_rre and not calculate_hashable_rre_early

    calculate_rre_early = calculate_hashable_rre_early or remove_obvious_collapses
    calculate_rre_late = (return_rre or return_rre_pivots) and not calculate_rre_early

    hashable_rre_seen = set()

    def calc_rre(mat):
        rre, pivots = mat.rref()
        stripped = sympy_tools.strip_zero_rows(rre)
        #print(f"{pivots},{repr(rre)}")
        return stripped, pivots

    # TODO: consider making prefix a SymPy matrix natively, so that we are not always converting, and can more easily get different views. Maybe this would speed somet hings up??
    def dfs(prefix, start_row):
        # "prefix" is a list or rows, each of which is a tuple. 
        # "mat" is a Sympy representation of prefix.

        # Yield the current matrix (if prefix is non-empty)
        if prefix:
            mat = sp.Matrix(prefix)

            if calculate_rre_early:
                rre, rre_pivots = calc_rre(mat)
 
            if remove_obvious_collapses:
                assert calculate_rre_early
                if sympy_tools.some_row_causes_collapse(rre, k):
                    if debug: print(f"VETO as row collapse in {rre} with k={k}.")
                    # Some row causes collapse!
                    # Skip deeper evaluation or return of it!
                    return
                if not sympy_tools.pivot_positions_are_all_viable_for_stripped_RRE_matrix(rre.shape, rre_pivots, k):
                    if debug: print(f"VETO as bad pivot positions in {rre} for k={k}.")
                    # Some row causes collapse!
                    # Skip deeper evaluation or return of it!
                    return

            if calculate_hashable_rre_early:
                hashable_rre = sp.ImmutableMatrix(rre)

            if remove_duplicates_via_hash:
                assert calculate_hashable_rre_early
                if hashable_rre in hashable_rre_seen:
                    if debug: print(f"VETO as already seen {rre}")
                    # We already saw this one, so don't need to produce it again!
                    # Skip deeper evaluation or return of it!
                    return
                else:
                    # record that we have seen this item:
                    hashable_rre_seen.add(hashable_rre) # Note this is a sort of voluntary memory leak. Users use this at their own risk!

            # Our own standard checks are complete! Now allow external user checks on mat. (TODO -- allow user to check RRE too?)

            if old_yield_matrix is not None and not old_yield_matrix(mat):
                return # Skip deeper exploration without yielding mat as the user's pre-yield test is not passed.


            EE, OO = None, None
            if confusable_sets_or_None_function is not None:
                confusable_sets_or_None = confusable_sets_or_None_function(mat)
                if confusable_sets_or_None == None:
                    return # Skip deeper exploration without yielding mat as the user's pre-yield test is not passed.

                assert len(confusable_sets_or_None) == 2
                # Let's extract the confusable sets:
                EE, OO = confusable_sets_or_None
                assert len(EE)==len(OO)
                assert len(EE)>0
                print(f"CONFUSABLE SET SIZE WAS {len(EE)}")

            # At this point we know we have to return things, so finish any late computations, if required:
            if calculate_rre_late:
                rre, rre_pivots = calc_rre(mat)
            if calculate_hashable_rre_late:
                hashable_rre = sp.ImmutableMatrix(rre)

            ans = []
            if return_mat:
                ans.append(mat)
            if return_rre:
                ans.append(rre)
            if return_rre_pivots:
                ans.append(rre_pivots)
            if return_hashable_rre:
                ans.append(hashable_rre)
            if return_confusable_sets:
                assert EE != None
                assert OO != None
                assert len(EE)==len(OO)
                assert len(EE)>0
                ans.append( (EE, OO) )
            
            # Now we pass all our outputs to the caller:
            user_aborted_this_branch = (yield ans)

            # Don't go deeper if caller says not to, or for any other reason:
            if (
                user_aborted_this_branch
                or
                (remove_obvious_collapses and debug_test_max_rows and rre.shape[0]>=max_rows) # We are already as tall as we ever can be!
                or
                (go_deeper is not None and not go_deeper(mat))
               ):
                return  # Skip deeper exploration

            columns_of_mat_as_tuples = [tuple(mat.col(i)) for i in range(mat.cols)]  # as tuples so that they will be hashable and thus usable as dictionary keys
            e_places = Equivalent_Places(exemplar = columns_of_mat_as_tuples)
        else:
            assert not prefix
            e_places = Equivalent_Places(size=M, all_equivalent=True)

        # Start the rows at the given start_row:
        row_gen = generate_all_vertex_matches_given_equivalent_places(equivalent_places = e_places, k=k, start=start_row)

        for row in row_gen:
            # Avoid repeating the start_row itself at the top of recursion
            if start_row is not None and row == start_row:
                continue
            # Recurse with the new row appended to prefix
            yield from dfs(prefix + [row], row)

    # Start with no prefix and no lower bound
    yield from dfs([], None)

def alpha_attacking_matrix(
                L_matrix : sp.Matrix,  # Any number of rows, but M columns.
                bat_matrix : sp.Matrix, # M rows, and k columns, so that each row is a bat vector.
                ) -> sp.Matrix :
    """
    This method generates the matrix A for which the solns of A.(vec of alphas) are the same as the solutions to L.(alpha1 w1, alpha2 w2, ... , alphaM, wM) where w1 is the first bat (i.e. first row of bat matrix) and w2 the second, and so on.
    """

    M, k = bat_matrix.shape 
    R, M_L = L_matrix.shape

    assert M == M_L, f"L and B must work with the same no. of vectors. Wanted M({M}) == M_L({M_L})"

    effective_row = 0
    ans = sp.zeros(R*k, M)
    for i in range(R):
        for kk in range(k):
            for j in range(M):
                ans[effective_row, j] = L_matrix[i,j]*bat_matrix[j,kk]
            effective_row += 1

    return ans

def demo():
    M=10
    print(f"All ***SIGNATURES*** given M={M} bad bats are:")
    for i, match in enumerate(generate_all_vertex_match_signatures(k=None, M=M)):
       print(f"   {i+1}:    {match}")
    print()

    M=10
    start=(0,3,7)
    print(f"The ***SIGNATURES*** starting at {start} given M={M} bad bats are:")
    for i, match in enumerate(generate_all_vertex_match_signatures(k=None, M=M, start=start)):
       print(f"   {i+1}:    {match}")
    print()

    M=10
    start=(4,3,3)
    print(f"The ***SIGNATURES*** starting at {start} given M={M} bad bats are:")
    for i, match in enumerate(generate_all_vertex_match_signatures(k=None, M=M, start=start)):
       print(f"   {i+1}:    {match}")
    print()

    M=4
    print(f"All matches given M={M} bad bats are:")
    for i, match in enumerate(generate_all_vertex_matches(k=None, M=M)):
       print(f"   {i+1}:    {match}")
    print()

    k=2

    print(f"All useful matches in k={k} dimensions, given M={M} bad bats, but ignoring permutations are:")
    for i,match in enumerate(generate_all_vertex_matches(k=k, M=M, permute=False)):
       print(f"   {i+1}:    {match}")
    print()



    for equivalent_places in ( Equivalent_Places(size=M, none_equivalent=True), ):
        print(f"All USEFUL matches in k={k} dimensions, given M={M} bad bats, for equivalent_places={equivalent_places}")
        for i,match in enumerate(generate_all_vertex_matches_given_equivalent_places(k=k, equivalent_places=equivalent_places)):
           print(f"   {i+1}:    {match}")
        print()

    for equivalent_places in ( Equivalent_Places(size=M, none_equivalent=True), ):
        print(f"All matches given M={M} bad bats, for equivalent_places={equivalent_places}")
        for i,match in enumerate(generate_all_vertex_matches_given_equivalent_places(k=None, equivalent_places=equivalent_places)):
           print(f"   {i+1}:    {match}")
        print()

    print("timing tests are now in timing_tests.py")





    print("== Test of Matrix Generation =========")
    def max_row_requirement(mat, max_rows):
        return sp.shape(mat)[0] <= max_rows

    mat_gen_slow = generate_viable_vertex_match_matrices(
        M=5,
        k=2,
        yield_matrix = partial(max_row_requirement, max_rows=4),
        ) 

    mat_gen_fast = generate_viable_vertex_match_matrices(
        M=5,
        k=2,
        go_deeper = partial(max_row_requirement, max_rows=3),
        ) 

    print("Will check if two methods agree:")
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


if __name__ == "__main__":
    M=7 #16 #7
    k=2 # 5 #2
    k=3 # 5 #2
    print(f"All useful matches in k={k} dimensions, given M={M} bad bats, but ignoring permutations are:")
    for i,match in enumerate(generate_all_vertex_matches(k=k, M=M, permute=False)):
       print(f"   {i+1}:    {match}")
    for i,match in enumerate(generate_all_vertex_match_signatures(k=k, M=M)):
       print(f"   {i+1}:    {match}")
    print()
    #demo()
