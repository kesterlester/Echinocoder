# This is taken from https://github.com/more-itertools/more-itertools/blob/edcafcfa58b1f2edc4b8588e98fba3f41784b746/more_itertools/more.py#L675 ... but with
# some modifications by Christopher Lester to add leftovers.

from functools import partial, total_ordering
from collections import defaultdict
from itertools import cycle

def distinct_permutations(iterable, r=None, output_leftovers=False):
    """Yield successive distinct permutations of the elements in *iterable*.

        >>> sorted(distinct_permutations([1, 0, 1]))
        [(0, 1, 1), (1, 0, 1), (1, 1, 0)]

    Equivalent to yielding from ``set(permutations(iterable))``, except
    duplicates are not generated and thrown away. For larger input sequences
    this is much more efficient.

    Duplicate permutations arise when there are duplicated elements in the
    input iterable. The number of items returned is
    `n! / (x_1! * x_2! * ... * x_n!)`, where `n` is the total number of
    items input, and each `x_i` is the count of a distinct item in the input
    sequence. The function :func:`multinomial` computes this directly.

    If *r* is given, only the *r*-length permutations are yielded.

        >>> sorted(distinct_permutations([1, 0, 1], r=2))
        [(0, 1), (1, 0), (1, 1)]
        >>> sorted(distinct_permutations(range(3), r=2))
        [(0, 1), (0, 2), (1, 0), (1, 2), (2, 0), (2, 1)]

    *iterable* need not be sortable, but note that using equal (``x == y``)
    but non-identical (``id(x) != id(y)``) elements may produce surprising
    behavior. For example, ``1`` and ``True`` are equal but non-identical:

        >>> list(distinct_permutations([1, True, '3']))  # doctest: +SKIP
        [
            (1, True, '3'),
            (1, '3', True),
            ('3', 1, True)
        ]
        >>> list(distinct_permutations([1, 2, '3']))  # doctest: +SKIP
        [
            (1, 2, '3'),
            (1, '3', 2),
            (2, 1, '3'),
            (2, '3', 1),
            ('3', 1, 2),
            ('3', 2, 1)
        ]
    """

    # Algorithm: https://w.wiki/Qai
    def _full(A):
        while True:
            # Yield the permutation we have
            yield tuple(A) if not output_leftovers else (tuple(A), ()) # Was just "yield tuple(A)" before CGL MOD

            # Find the largest index i such that A[i] < A[i + 1]
            for i in range(size - 2, -1, -1):
                if A[i] < A[i + 1]:
                    break
            #  If no such index exists, this permutation is the last one
            else:
                return

            # Find the largest index j greater than j such that A[i] < A[j]
            for j in range(size - 1, i, -1):
                if A[i] < A[j]:
                    break

            # Swap the value of A[i] with that of A[j], then reverse the
            # sequence from A[i + 1] to form the new permutation
            A[i], A[j] = A[j], A[i]
            A[i + 1 :] = A[: i - size : -1]  # A[i + 1:][::-1]

    # Algorithm: modified from the above
    def _partial(A, r):
        # Split A into the first r items and the last r items
        head, tail = A[:r], A[r:]
        right_head_indexes = range(r - 1, -1, -1)
        left_tail_indexes = range(len(tail))

        while True:
            # Yield the permutation we have
            if output_leftovers: # CGL ADDED THIS LINE
                yield tuple(head), tuple(tail) # CGL ADDED THIS LINE
            else: # CGL ADDED THIS LINE
                yield tuple(head) # CGL INDENTED THIS LINE

            # Starting from the right, find the first index of the head with
            # value smaller than the maximum value of the tail - call it i.
            pivot = tail[-1]
            for i in right_head_indexes:
                if head[i] < pivot:
                    break
                pivot = head[i]
            else:
                return

            # Starting from the left, find the first value of the tail
            # with a value greater than head[i] and swap.
            for j in left_tail_indexes:
                if tail[j] > head[i]:
                    head[i], tail[j] = tail[j], head[i]
                    break
            # If we didn't find one, start from the right and find the first
            # index of the head with a value greater than head[i] and swap.
            else:
                for j in right_head_indexes:
                    if head[j] > head[i]:
                        head[i], head[j] = head[j], head[i]
                        break

            # Reverse head[i + 1:] and swap it with tail[:r - (i + 1)]
            tail += head[: i - r : -1]  # head[i + 1:][::-1]
            i += 1
            head[i:], tail[:] = tail[: r - i], tail[r - i :]

    items = list(iterable)

    try:
        items.sort()
        sortable = True
    except TypeError:
        sortable = False

        # New implementation by CGL.
        @total_ordering
        class Wrapper(object):
            def __init__(self, sorting_cue, payload):
                # sorting_cue should be something for which comparisons and sorting work, as surrogate for the payload
                self.sorting_cue = sorting_cue
                self.payload = payload

            def __eq__(self, other):
                return self.sorting_cue == other.sorting_cue

            def __lt__(self, other):
                return self.sorting_cue < other.sorting_cue

        items = [Wrapper(items.index(item), item) for item in items ]

    size = len(items)
    if r is None:
        r = size

    # functools.partial(_partial, ... )
    algorithm = _full if (r == size) else partial(_partial, r=r)

    if 0 < r <= size:
        if sortable:
            yield from algorithm(items)
        else:
            if output_leftovers:
                for wrapped_items, wrapped_leftovers in algorithm(items):
                    yield [wrapped_item.payload for wrapped_item in wrapped_items], [wrapped_leftover.payload for wrapped_leftover in wrapped_leftovers]     
            else:
                for wrapped_items in algorithm(items):
                    yield [wrapped_item.payload for wrapped_item in wrapped_items],     

    return iter(() if r else ((),))


def tost():
    thing=[3,0,3]
    print("distinct perms of",thing, "are")
    expected_perms = [ (0,3,3), (3,0,3), (3,3,0) ]
    for n, perm in enumerate(distinct_permutations(thing)):
        print(f"Expected perm {expected_perms[n]} got perm {perm}.")
        assert expected_perms[n] == perm

if __name__ == "__main__":
    tost()

