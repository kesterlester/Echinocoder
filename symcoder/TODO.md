Move def _force_positive_side_key(op, flavour): out of symcoder/encode.py and into symatom somwhere?   Yes it's valid where it already is, but conceptually if there is a need for sign-blind atoms or operations (this is not the same thing as taking their modulus when evaluated as +eps(a,b,c) can still evaluate to a negative number!) then it is porbably better to have then in symatom as the current implementation is using implementation details of the symatom FO representation to achieve its end.  It is potentially a worry that if the symatom FO or atom back-end changedm it would break this piece of code in symcoder.

Is there a better way of doing the 
print("BEWARE: a deduplication event fired in encode.py despite it being thought that this could/should never happen with the implementation at time of coding.") in [encode.py](encode.py) ?
            

Get rid of the 8-fold deduplication in [encode.py](encode.py) if it is really not needed.


(1) I like packages to have simple demo routines (like symatom/demo.p) with sample output also stored in the repsoitory, o show what's going on.  We are on the point of being able to have such a thing for this encoding. It would start with some numerical vectors a=[2,4,-2], b=[2,6,4],p=[-1,-1,1] or whatever and an E={{"a","b"},{"p"}} plan of some kind (maybe more complex than this one) and would then encode this to  perm and rot invariant numbers, but saying what it is doing on the way. E,g, "FO(blah) for rows 3+5 with canonical orbit representative are BLH BLAH end encode to BLAH" then nex tline .etc.  If we went for such a thing, it would need some discussion on what to use to keep things complex enough that a human can folow through the logic, but not so simple as to make the encoding not show critical features of the algorithm.  Probably a way of achieving that would be to work in a 2D system where eps2() is a thing.

(2) Is symatom/SPEC.md  up to date?

(3) At present symcoder doesn't hava a SPEC.md at all --- but maybe that's because it's evolving and we should write the spec for it after it has stoped evolving.

(4) There is no top-level documentation explaining how symatom and symcode responsibilites are split. Some parts of this is indirectly explained in lower level docs, but there is no top-down attempt.

(5) There are two potentially major disruptive changes I have mostly kept quiet about so far (for fear of creating confusion). They have the possibility of cutting the embedding size down a lot.  One relates to signs, one relates to pre-row knowledge. Either of them could be discussed.

(6) If there was a demo for symcoder, it is likely that the end-to-end-output could potentially become a unit test BUT this is also a bad idea as there are lots of allowed choices (types of canoncialisation, and optimisations hinted at in (5) above, that would change the output of an end to end test by removing redundancy, so maybe we don't really want to over-fixate on having tests that will break when output becomes less redund
