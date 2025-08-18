from distinct_combinations_with_start import distinct_combinations_with_start as LESTER # Supports "start".
from more_itertools import distinct_combinations as ITERTOOLS # Does not support "start".
from itertools import zip_longest
from itertools import chain

class X:
    def __init__(self, name): self.name = name
    def __repr__(self): return f"X({self.name!r})"

a = X('a')
b = X('b')
c = X('c')

def test_things():
    for data, start in (
         ((1,2,2,3), (2,)),
         (("S","p","e","e","d","o"),  ("p", "d", "o")),
         (("Christopher"),  tuple("hrhr")),
         ((a,b,c,b,c),    (c,)),
         ((a,b,b,c),    (b,b)),
         ((a,b,b,c),    (b,b)),
         ):
         r = len(start)
         print("==============================")
         print(f"Starting test with data={data}, r={r} and start (where used) of {start}.")
         print("==============================")
         start_pos = None
         for i, (lesters, itertoolss) in enumerate(zip_longest(LESTER(data, r), ITERTOOLS(data, r))):
             print(f"{i}:         {lesters} {'==' if lesters==itertoolss else '!='} {itertoolss}")
             if lesters == start and start_pos == None:
                start_pos = i
             assert lesters==itertoolss
         print("Match confirmed!")
         assert start_pos is not None
         print(f"Start pos determined to be {start_pos}.")
         for i, (lesters, itertoolss) in enumerate(zip_longest(chain(iter([None,]*start_pos),LESTER(data, r, start=start)), ITERTOOLS(data, r))):
             if i < start_pos:
                 print(f"{i}:         {lesters} ...  {itertoolss}")
             else:
                 print(f"{i}:         {lesters} {'==' if lesters==itertoolss else '!='} {itertoolss}")
                 assert lesters==itertoolss

if __name__ == "__main__":
    test_things()
