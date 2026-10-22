# Embed list of real  numbers treated as multiset.
# E.g. this method can embed things like:
#
#        [3,4,-2,]
#
# representing the multiset
#
#        {{ 3,4,-2, }}
#
# into
#
#        [-2, 3, 4].

def embed(data, debug=False):
    return sorted(data), len(data), None
