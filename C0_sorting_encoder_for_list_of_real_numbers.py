# Encode list of real  numbers treated as multiset.
# E.g. this method can encode things like:
#
#        [3,4,-2,]
#

name="C0_sorting_encoder"

def encode(data):
    return sorted(data)
