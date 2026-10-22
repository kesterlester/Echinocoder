from distinct_partitions_v1_no_start import distinct_partitions_v1_no_start
from distinct_partitions_with_start import distinct_partitions_with_start
from itertools import zip_longest, chain

def test_distinct_partitions():

    def tuple_of_strings_converter(tup_of_strings):
        return tuple( ( tuple(ch for ch in s) for s in tup_of_strings)  )

    test_programme = [
        ( "make", (2,2), (
                      ("ma", "ke"), # Note that the other way around ("ke", "ma") appears later too. That is intentional.
                      ("mk", "ae"),
                      ("me", "ak"),
                      ("ak", "me"),
                      ("ae", "mk"),
                      ("ke", "ma"),
                     ), tuple_of_strings_converter, "set",
                     ("me","ak"),
        ),
        ( "make", (4,0), (
                      ("make", ""), 
                     ), tuple_of_strings_converter, "set",
                     ("make", ""),
        ),
        ( "make", (0,4), (
                      ("", "make"), 
                     ), tuple_of_strings_converter, "set",
                      ("", "make"), 
        ),
        ( "happy", (2,3), (
                      ("ha", "ppy"),
                      ("hp", "apy"),
                      ("hy", "app"),
                      ("ap", "hpy"),
                      ("ay", "hpp"),
                      ("pp", "hay"),
                      ("py", "hap"),
                     ), tuple_of_strings_converter, "set",
                      ("hy", "app"),
        ),
        ( "happy", (2,3), (
                      (('h', 'a'), ('p', 'p', 'y')),
                      (('h', 'p'), ('a', 'p', 'y')),
                      (('h', 'y'), ('a', 'p', 'p')),
                      (('a', 'p'), ('h', 'p', 'y')),
                      (('a', 'y'), ('h', 'p', 'p')),
                      (('p', 'p'), ('h', 'a', 'y')),
                      (('p', 'y'), ('h', 'a', 'p')),
                          ), None, "exact",
                      (('a', 'y'), ('h', 'p', 'p')),
        ),
        ( "Spain", (2,3), (
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
                          ), None, "exact",
            (('p', 'i'), ('S', 'a', 'n')),
        ),
        ( "happy", (2,1,2), (
(('h', 'a'), ('p',), ('p', 'y')),
(('h', 'a'), ('y',), ('p', 'p')),
(('h', 'p'), ('a',), ('p', 'y')),
(('h', 'p'), ('p',), ('a', 'y')),
(('h', 'p'), ('y',), ('a', 'p')),
(('h', 'y'), ('a',), ('p', 'p')),
(('h', 'y'), ('p',), ('a', 'p')),
(('a', 'p'), ('h',), ('p', 'y')),
(('a', 'p'), ('p',), ('h', 'y')),
(('a', 'p'), ('y',), ('h', 'p')),
(('a', 'y'), ('h',), ('p', 'p')),
(('a', 'y'), ('p',), ('h', 'p')),
(('p', 'p'), ('h',), ('a', 'y')),
(('p', 'p'), ('a',), ('h', 'y')),
(('p', 'p'), ('y',), ('h', 'a')),
(('p', 'y'), ('h',), ('a', 'p')),
(('p', 'y'), ('a',), ('h', 'p')),
(('p', 'y'), ('p',), ('h', 'a')),
                          ), None, "exact",
(('p', 'p'), ('y',), ('h', 'a')),
        ),
        ( "Spain", (2,1,2), (
(('S', 'p'), ('a',), ('i', 'n')),
(('S', 'p'), ('i',), ('a', 'n')),
(('S', 'p'), ('n',), ('a', 'i')),
(('S', 'a'), ('p',), ('i', 'n')),
(('S', 'a'), ('i',), ('p', 'n')),
(('S', 'a'), ('n',), ('p', 'i')),
(('S', 'i'), ('p',), ('a', 'n')),
(('S', 'i'), ('a',), ('p', 'n')),
(('S', 'i'), ('n',), ('p', 'a')),
(('S', 'n'), ('p',), ('a', 'i')),
(('S', 'n'), ('a',), ('p', 'i')),
(('S', 'n'), ('i',), ('p', 'a')),
(('p', 'a'), ('S',), ('i', 'n')),
(('p', 'a'), ('i',), ('S', 'n')),
(('p', 'a'), ('n',), ('S', 'i')),
(('p', 'i'), ('S',), ('a', 'n')),
(('p', 'i'), ('a',), ('S', 'n')),
(('p', 'i'), ('n',), ('S', 'a')),
(('p', 'n'), ('S',), ('a', 'i')),
(('p', 'n'), ('a',), ('S', 'i')),
(('p', 'n'), ('i',), ('S', 'a')),
(('a', 'i'), ('S',), ('p', 'n')),
(('a', 'i'), ('p',), ('S', 'n')),
(('a', 'i'), ('n',), ('S', 'p')),
(('a', 'n'), ('S',), ('p', 'i')),
(('a', 'n'), ('p',), ('S', 'i')),
(('a', 'n'), ('i',), ('S', 'p')),
(('i', 'n'), ('S',), ('p', 'a')),
(('i', 'n'), ('p',), ('S', 'a')),
(('i', 'n'), ('a',), ('S', 'p')),
                          ), None, "exact",
(('p', 'a'), ('n',), ('S', 'i')),
        ),
        ( "Spain", (2,3), (
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
                          ), None, "exact",
            (('S', 'p'), ('a', 'i', 'n')),
        ),
        ( "make", (5,), (
                      ("", "make"), 
                     ), tuple_of_strings_converter, "throws",
                     None,
        ),
        ( "make", (5,-1), (
                      ("", "make"), 
                     ), tuple_of_strings_converter, "throws",
                     None,
        ),
        ( "make", (-1), (
                      ("", "make"), 
                     ), tuple_of_strings_converter, "throws",
                     None,
        ),
        ( "make", "a", (
                      ("", "make"), 
                     ), tuple_of_strings_converter, "throws",
                     None,
        ),
    ]

  
    
    for data, splitting, expected, converter, typ, start in test_programme:
        print("==============================")
        print(f"Starting test with data={data}, splitting={splitting}, type={typ}  and start (where used) of {start}.")
        print("==============================")
        if typ=="throws":
            try:
               distinct_partitions_v1_no_start(data, splitting)
               assert True
            except:
               assert False
            continue

        result_exact = distinct_partitions_v1_no_start(data, splitting)
        result_list = tuple(distinct_partitions_v1_no_start(data, splitting))
        if converter is not None:
            expected = tuple(map( converter, expected ))

        if typ=="exact":
            start_pos = None
            for i, (lesters, itertoolss) in enumerate(zip_longest(distinct_partitions_with_start(data, splitting), distinct_partitions_v1_no_start(data, splitting))):
                print(f"{i}:         {lesters} {'==' if lesters==itertoolss else '!='} {itertoolss}")
                if lesters == start and start_pos == None:
                   start_pos = i
                assert lesters==itertoolss
            print("Match confirmed!")
            assert start_pos is not None
            print(f"Start pos determined to be {start_pos}.")
            for i, (lesters, itertoolss) in enumerate(zip_longest(chain(iter([None,]*start_pos),distinct_partitions_with_start(data, splitting, start=start)), distinct_partitions_v1_no_start(data, splitting))):
                if i < start_pos:
                    print(f"{i}:         {lesters} ...  {itertoolss}")
                else:
                    print(f"{i}:         {lesters} {'==' if lesters==itertoolss else '!='} {itertoolss}")
                    assert lesters==itertoolss


            continue

        if typ=="set":
            print(f"got {result_list}")
            print(f"expected {expected}")
            assert set(result_list) == set(expected)
            assert len(result_list) == len(expected)
            continue

        assert False



