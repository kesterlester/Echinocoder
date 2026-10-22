# Embed list of real numbers treated as multiset.
# E.g. this method can embed things like:
#
#        [3,4,-2,]
#
# representing the multiset
#
#       {{ 3, 4, -2 }}
#
# Although this implementation claims to only do reals, it does in fact (privately) embed complex lists, albeit to complex outputs.

import numpy as np

def embed(data, debug=False):
    if debug:
        print("Asked to encode ",data)
    return np.polynomial.polynomial.polyfromroots(-data)[-2::-1], len(data), None  # The -1 in [-2::-1] reverses the order of the list so that the terms linear in the roots come first. The -2 at the front forced the list to start from the coeffient for x^(n-1). (The coefficient of x^n is always 1 and we dont need it!)


    ## Alternative:
    #
    # import tools
    # ans = np.polynomial.polynomial.polyfromroots(-data)[-2::-1] # The -1 in [-2::-1] reverses the order of the list so that the terms linear in the roots come first. The -2 at the front forced the list to start from the coeffient for x^(n-1). (The coefficient of x^n is always 1 and we dont need it!)

    ## All embedders have to output lists of real numbers (at least for now) so:
    #if np.iscomplexobj(ans):
    #  ans=tools.expand_complex_to_real_pairs(ans)
    #
    #return ans
