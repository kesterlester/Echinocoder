from more_itertools import unique_everseen
from operator import itemgetter
from more_itertools import distinct_combinations 
from bisect import bisect_left

#### Here is the implementation from itertools before it got modified to add "start"
####
#### def distinct_combinations_as_more_itertools(iterable, r):
####     """Yield the distinct combinations of *r* items taken from *iterable*.
#### 
####         >>> list(distinct_combinations([0, 0, 1], 2))
####         [(0, 0), (0, 1)]
#### 
####     Equivalent to ``set(combinations(iterable))``, except duplicates are not
####     generated and thrown away. For larger input sequences this is much more
####     efficient.
#### 
####     """
####     if r < 0:
####         raise ValueError('r must be non-negative')
####     elif r == 0:
####         yield ()
####         return
####     pool = tuple(iterable)
####     generators = [unique_everseen(enumerate(pool), key=itemgetter(1))]
####     current_combo = [None] * r
####     level = 0
####     while generators:
####         try:
####             cur_idx, p = next(generators[-1])
####         except StopIteration:
####             generators.pop()
####             level -= 1
####             continue
####         current_combo[level] = p
####         if level + 1 == r:
####             yield tuple(current_combo)
####         else:
####             generators.append(
####                 unique_everseen(
####                     enumerate(pool[cur_idx + 1 :], cur_idx + 1),
####                     key=itemgetter(1),
####                 )
####             )
####             level += 1

def distinct_combinations_with_start(iterable, r, start=None):
    """Yield the distinct combinations of r items taken from iterable.

    If `start` is provided: (1) it must bea tuple that this generator would normally produce,
    and (2) the generator will then begin yielding from that combination efficently -- 
    i.e. without calculating and discarding earlier results.

    Examples:
        >>> list(distinct_combinations_with_start("hello", 2))
        [('h','e'), ('h','l'), ('h','o'), ('e','l'), ('e','o'), ('l','l'), ('l','o')]

        >>> list(distinct_combinations_with_start("hello", 2, start=('l','l')))
        [('l','l'), ('l','o')]
    """
    if r < 0:
        raise ValueError("r must be non-negative")
    if r == 0:
        yield ()
        return

    pool = tuple(iterable)
    n = len(pool)
    if r > n:
        return

    def make_gen(from_index):
        # Generator over first occurrences (by value) in the suffix [from_index:]
        return unique_everseen(enumerate(pool[from_index:], from_index), key=itemgetter(1))

    # No start: behave like the normal distinct-combinations generator
    if start is None:
        generators = [make_gen(0)]
        current_combo = [None] * r
        level = 0
        while generators:
            try:
                cur_idx, p = next(generators[-1])
            except StopIteration:
                generators.pop()
                level -= 1
                continue
            current_combo[level] = p
            if level + 1 == r:
                yield tuple(current_combo)
            else:
                generators.append(make_gen(cur_idx + 1))
                level += 1
        return

    # Validate start length
    if len(start) != r:
        raise ValueError("start must be a tuple of length r")

    # Precompute occurrence lists for each value: value -> [indices...]
    occ = {}
    for i, v in enumerate(pool):
        occ.setdefault(v, []).append(i)

    # Map `start` values to the earliest valid indices i0 < i1 < ... using bisect
    index_path = []
    last_index = 0
    for v in start:
        if v not in occ:
            raise ValueError(f"value {v!r} from start not found in iterable")
        idx_list = occ[v]
        pos = bisect_left(idx_list, last_index)  # first occurrence >= last_index
        if pos == len(idx_list):
            # No occurrence of v at or after last_index -> invalid start
            raise ValueError(f"start {start!r} is not a valid distinct combination for this iterable")
        idx = idx_list[pos]
        index_path.append(idx)
        last_index = idx + 1

    # Recreate the stack of generators positioned at the start path.
    # For levels 0..r-2: advance their generators *through* index_path[level].
    # For the final level r-1: also advance to consume its element so we can yield `start`
    # once, then continue naturally afterward.
    generators = []
    current_combo = [None] * r

    for level, idx in enumerate(index_path):
        base_from = 0 if level == 0 else index_path[level - 1] + 1
        g = make_gen(base_from)

        # Advance this generator until it yields the pair with cur_idx == idx
        for cur_idx, p in g:
            if cur_idx == idx:
                current_combo[level] = p  # record chosen value at this level
                break
        else:
            # Should never happen if index_path is consistent
            raise RuntimeError("Failed to align generator with start index.")

        generators.append(g)

    # Yield the start tuple exactly once
    yield tuple(current_combo)

    # Continue the normal generation loop from the reconstructed state
    level = r - 1
    while generators:
        try:
            cur_idx, p = next(generators[-1])
        except StopIteration:
            generators.pop()
            level -= 1
            continue
        current_combo[level] = p
        if level + 1 == r:
            yield tuple(current_combo)
        else:
            generators.append(make_gen(cur_idx + 1))
            level += 1

def demo():

    print("=================== hello r=2 more_itertools ====================")
    for i, x in enumerate(distinct_combinations("hello", 2)):
        print(f"{i}:     {x}")

#    print("=================== hello r=2 as_more_itertools ====================")
#    for i, x in enumerate(distinct_combinations_as_more_itertools("hello", 2)):
#        print(f"{i}:     {x}")

    print("=================== hello r=2 with_start ====================")
    for i, x in enumerate(distinct_combinations_with_start("hello", 2)):
        print(f"{i}:     {x}")

    print("=================== hello r=2 start=('e', 'l') ====================")
    for i, x in enumerate(distinct_combinations_with_start("hello", 2, start=('e', 'l'))):
        print(f"{i}:     {x}")

    print("=================== hello r=2 more_itertools ====================")
    for i, x in enumerate(distinct_combinations("Christopher", 4)):
        print(f"{i}:     {x}")

    print("=================== hello r=2 start=('h', 'i', 'p', 'h') ====================")
    for i, x in enumerate(distinct_combinations_with_start("Christopher", 4, start=('h', 'i', 'p', 'h'))):
        print(f"{i}:     {x}")

    print("=================== hello r=2 start=('s', 't', 'h', 'r') ====================")
    for i, x in enumerate(distinct_combinations_with_start("Christopher", 4, start=('s', 't', 'h', 'r'))):
        print(f"{i}:     {x}")

if __name__ == "__main__":
    demo()
