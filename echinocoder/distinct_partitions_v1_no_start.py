from more_itertools import distinct_combinations

def distinct_partitions_v1_no_start(iterable, splitting, validate_splitting=True):
    
    """
    This returns a generator which will yield all the partitions of "iterable" which split it according to
    the sizes given in the (ordered!) tuple "splitting".
    Within each group, the order is preserved according to first-come first-served in the original iterable.

    splitting shall be a tuple of non-negative integers whose sum shall be the length of the iterable.

    The word "distinct" is part of the name of the function to indicate that identical items in the input
    are not distinguished.  See how the example below is affected by the repeated "l" in "hello" to explain this further.

    By default, the function validdates the splitting input against the iterable. However, if you are confident
    that your splitting and iterable inputs are compatible, you may turn off that validation for some extra
    speed. Note, however, that undefined behaviour may ensue if your confidence was misplaced!

    The result for

        distinct_partitions( "hello", (2,3) )

    should be an iterator which would successively yield the 7 elements of this tuple:

        (
            (('h', 'a'), ('p', 'p', 'y')),
            (('h', 'p'), ('a', 'p', 'y')), # Appears only once despite p in first group and p in second group!
            (('h', 'y'), ('a', 'p', 'p')),
            (('a', 'p'), ('h', 'p', 'y')), # Appears only once despite p in first group and p in second group!
            (('a', 'y'), ('h', 'p', 'p')),
            (('p', 'p'), ('h', 'a', 'y')), 
            (('p', 'y'), ('h', 'a', 'p')), # Appears only once despite p in first group and p in second group!
        )

    Note that had we used a five-letter word without repeated letters, such as "Spain", then
    we would have found 5-choose-2 = 10 elements in the list not 7 as with "hello":

        distinct_partitions( "Spain", (2,3) )

        iterates over:

        (
            (('S', 'p'), ('a', 'i', 'n')),
            (('S', 'a'), ('p', 'i', 'n')),
            (('S', 'i'), ('p', 'a', 'n')),
            (('S', 'n'), ('p', 'a', 'i')),
            (('p', 'a'), ('S', 'i', 'n')),
            (('p', 'i'), ('S', 'a', 'n')),
            (('p', 'n'), ('S', 'a', 'i')),
            (('a', 'i'), ('S', 'p', 'n')),
            (('a', 'n'), ('S', 'p', 'i')),
            (('i', 'n'), ('S', 'p', 'a')),
        )

    This non-repetition of repeated elements is what distinguishes distinct_partitions
    from other partitioners which consider all elements of the input to be distinguishable
    (e.g. more-itertools.partitions).

    """

    items = tuple(iterable)  # fixed-size storage

    if validate_splitting:
        if not isinstance(splitting, tuple):
            raise ValueError("distinct_partitions requires that splitting be a tuple.")
            
        total = 0
        for count in splitting:
            if int(count) != count or count < 0:
                raise ValueError("distinct_partitions requires that splitting is a tuple of non-negative integers.")
            total += count

        if len(items) != total:
            raise ValueError(f"distinct_partitions requires that len(items)==sum(splitting). You have {len(items)} != {total}.")

    if len(splitting) == 1:
        ans = (tuple(iterable),)
        yield ans
        return 

    assert len(splitting)>1
   
    for comb in distinct_combinations(items, splitting[0]):
        first_group = tuple(comb)

        # Build complement without storing all seen
        remaining_items = list(items)
        for c in first_group:
            remaining_items.remove(c)  # safe: remove one occurrence per element in comb

        for other_groups in distinct_partitions_v1_no_start(remaining_items, splitting[1:], validate_splitting=False): # Don't need to recheck splitting.
            yield (first_group,) + other_groups


def demo():
    print("distinct_partitions.py demo:")
    print()
    for splitting in ((2,3), (2,1,2)):
        for data in ("happy", "Spain"):
            print(f"distinct_partitions({data}, {splitting}) gives:")
            for split in distinct_partitions_v1_no_start(data, splitting):
                print(f"{split},")


if __name__ == "__main__":
    demo()
