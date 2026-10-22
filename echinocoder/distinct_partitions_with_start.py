from distinct_combinations_with_start import distinct_combinations_with_start

def distinct_partitions_with_start(iterable, splitting, validate_splitting=True, start=None):
    """
    Yield all the partitions of `iterable` split according to `splitting` (tuple of group sizes).
    Order is determined by first-come-first-served in the original iterable.
    Identical elements are not distinguished.

    If `start` is provided, it must be a tuple this generator would normally yield,
    and generation will begin there directly without generating/discarding earlier results.

    splitting shall be a tuple of non-negative integers whose sum shall be the length of the iterable.

    The word "distinct" is part of the name of the function to indicate that identical items in the input
    are not distinguished.  See how the example below is affected by the repeated "l" in "hello" to explain this further.

    By default, the function validdates the splitting input against the iterable. However, if you are confident
    that your splitting and iterable inputs are compatible, you may turn off that validation for some extra
    speed. Note, however, that undefined behaviour may ensue if your confidence was misplaced!

    The result for

        distinct_partition_with_start( "hello", (2,3) )

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

    items = tuple(iterable)

    if validate_splitting:
        if not isinstance(splitting, tuple):
            raise ValueError("splitting must be a tuple.")
        if any((not isinstance(c, int)) or c < 0 for c in splitting):
            raise ValueError("splitting must be a tuple of non-negative integers.")
        total = sum(splitting)
        if total != len(items):
            raise ValueError(
                f"len(items) ({len(items)}) != sum(splitting) ({total}) when items={items} and splitting={splitting}."
            )

    if len(splitting) == 1:
        ans = (tuple(iterable),)
        if start is None or ans >= start:
            yield ans
        return

    first_group_target = start[0] if start else None
    skipping = start is not None

    comb_gen = distinct_combinations_with_start(
        items, splitting[0],
        start=first_group_target if skipping else None
    )

    for first_group in comb_gen:
        remaining_items = list(items)
        for c in first_group:
            remaining_items.remove(c)

        if skipping and first_group == first_group_target:
            new_start = start[1:]
        else:
            new_start = None
            skipping = False

        for other_groups in distinct_partitions_with_start(
            remaining_items, splitting[1:],
            validate_splitting=False,
            start=new_start
        ):
            yield (first_group,) + other_groups


def demo():
    print("distinct_partitions.py demo:")
    print()
    for splitting in ((2, 3), (2, 1, 2)):
        for data in ("happy", "Spain"):
            print(f"distinct_partitions_with_start({data}, {splitting}) gives:")
            for split in distinct_partitions_with_start(data, splitting):
                print(f"{split},")
    print("\nDemo starting from a given partition:")
    start_point = (('p', 'i'), ('S', 'a', 'n'))
    for split in distinct_partitions_with_start("Spain", (2, 3), start=start_point):
        print(split)


if __name__ == "__main__":
    demo()
