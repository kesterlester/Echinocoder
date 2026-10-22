# Embed list of real numbers treated as multiset.
# E.g. this method can embed things like:
#
#        [3,4,-2,]
#
# representing the multiset
#
#        {{ 3,4,-2 }}
#
# Although this implementation claims to only do reals, it does in fact (privately) embed complex lists, albeit to complex outputs.
#
# This implementation as it has TERRIBLE run-time scaling with the length of the imput list, so don't 
# use it unless your your input list has at most ~15 elements.  Use the Cinf_numpy_... version for long lists.
# This implementation can work with arbitrary precision integers, however, which the Cinf_numpy_.... implementation cannot.

from math import prod
from itertools import combinations
#import tools


def embed(data, debug=False):
    if debug:
        print(f"Thinking about {data}")

    ans = [ sum( [ prod(c)  for c in combinations(data, r+1) ] ) for r in range(len(data)) ], len(data), None
    if debug:
        print(f"Expect to produce {ans}.")

    return ans
    
    # Alternative
    # 
    # ans = np.array([ sum( [ prod(c)  for c in combinations(data, r+1) ] ) for r in range(len(data)) ])
    #
    # All embedders have to output lists of real numbers (at least for now) so:
    #if np.iscomplexobj(ans):
    #  ans=tools.expand_complex_to_real_pairs(ans)
    #
    # return ans



