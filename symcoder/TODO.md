Move def _force_positive_side_key(op, flavour): out of symcoder/encode.py and into symatom somwhere?   Yes it's valid where it already is, but conceptually if there is a need for sign-blind atoms or operations (this is not the same thing as taking their modulus when evaluated as +eps(a,b,c) can still evaluate to a negative number!) then it is porbably better to have then in symatom as the current implementation is using implementation details of the symatom FO representation to achieve its end.  It is potentially a worry that if the symatom FO or atom back-end changedm it would break this piece of code in symcoder.

Is there a better way of doing the 
print("BEWARE: a deduplication event fired in encode.py despite it being thought that this could/should never happen with the implementation at time of coding.") in [encode.py](encode.py) ?
            

Get rid of the 8-fold deduplication in [encode.py](encode.py) if it is really not needed.
