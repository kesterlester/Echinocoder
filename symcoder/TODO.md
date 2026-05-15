## Make a fully fledged encoder object to pass through in metadata, rather than creating capability and then passing that back, as a user could tweak the capability (maybe) breaking the encoding. or he wrong capability could be sent back ....

## Move 11/12/21/22 stuff out of group.py

## remove-pair-with-self or do something better

This next OVERLAP_BLOCK shows that pair-with-self is still sneaking in, not deleted by 5c prining as here there are four jets {u,v,w,x}.

[120:132]  dot×dot  SS  u=(0,0,2)  v=(0,0,2)  shared=(0,0,0)  len=12  | dot(u,v)×dot(w,x)
    [4.2750, 4.2750, 0.0000, 14.3759, -6.4346, 6.4346, 0.5108, 0.0000, 15.5026, 15.5026, 0.0000, -55.0714]
[132:132]  dot×dot  SS  u=(0,0,2)  v=(0,0,2)  shared=(0,0,1)  len=0  NULL_ENCODING(deducible_from_uv_overlap_block)  | dot(u,v)×dot(u,w)
    [(empty)]
[132:144]  dot×dot  SS  u=(0,0,2)  v=(0,0,2)  shared=(0,0,2)  len=12  | dot(u,v)×dot(u,v)
    [4.2750, 4.2750, 0.0000, 1.9035, 17.6316, -17.6316, 6.5976, 0.0000, -14.7793, -14.7793, 0.0000, 3.1666]

Notes to self:

I was surprised that the code was ever iterating over row-paired-with-self because I had never done that "by hand" since indeed it adds no information. Now, in our ~few hours ago discussion on that point I mentioned that these pair-with-self rows are thrown out by 5c since "pair-with-self rows will form monolithic OVERLAP_BLOCKS.

We can see now that I was only half right.I was right that self-pair rows add no information and should be thrown out, but I was wrong to say that they alwaus form monolithic OVERLAP_BLOCKS and so get thrown out by 5c.  I was fooled into thinking that at the time as I was looking at a dot(p,q) with dot(p,q) pairing where p and q where the only vectors in the muons. WHereas the example above (with four vectors in the jets {u,v,w,x}) shows that pair-with-self of dot(u,v) is NOT thrown out by 5c as the extra jets provide extra overlap structure.

So, this tells us that we could modify 5c so that it is not just:

(i) "throw out the largest ASSOCIATION within each OVERLAP BLOCK as it is fixed by complementarity, but keep the N-1 others."

but can instead be:

(i) "within each OVERLAP BLOCK throw out every pair-with-self ASSOCIATION as it contains no information byeond that already in the orbit storage for that row,
(ii) then keep all but the largest of the ASSOCIATIONS that remain, if any remain.

That certainly could be done (though it would need edits to the latex encode.tex pedagogical document in parallel to explain why and when this happens and why it is correct) .... BUT ... it's just a particular special case. pair-with-self is uniquely bad and so should be thrown out.  Your correct point above is that the pair-not-quite-with-self things still have a lot of structure, and so can likely be not thrown out but at least reduced themselves by some as-yet not defined process.


## Idea

Move def _force_positive_side_key(op, flavour): out of symcoder/encode.py and into symatom somwhere?   Yes it's valid where it already is, but conceptually if there is a need for sign-blind atoms or operations (this is not the same thing as taking their modulus when evaluated as +eps(a,b,c) can still evaluate to a negative number!) then it is porbably better to have then in symatom as the current implementation is using implementation details of the symatom FO representation to achieve its end.  It is potentially a worry that if the symatom FO or atom back-end changedm it would break this piece of code in symcoder.

Is there a better way of doing the 
print("BEWARE: a deduplication event fired in encode.py despite it being thought that this could/should never happen with the implementation at time of coding.") in [encode.py](encode.py) ?
            

Get rid of the 8-fold deduplication in [encode.py](encode.py) if it is really not needed.


(1) I like packages to have simple demo routines (like symatom/demo.p) with sample output also stored in the repsoitory, o show what's going on.  We are on the point of being able to have such a thing for this encoding. It would start with some numerical vectors a=[2,4,-2], b=[2,6,4],p=[-1,-1,1] or whatever and an E={{"a","b"},{"p"}} plan of some kind (maybe more complex than this one) and would then encode this to  perm and rot invariant numbers, but saying what it is doing on the way. E,g, "FO(blah) for rows 3+5 with canonical orbit representative are BLH BLAH end encode to BLAH" then nex tline .etc.  If we went for such a thing, it would need some discussion on what to use to keep things complex enough that a human can folow through the logic, but not so simple as to make the encoding not show critical features of the algorithm.  Probably a way of achieving that would be to work in a 2D system where eps2() is a thing.

(2) Is symatom/SPEC.md  up to date?

(3) At present symcoder doesn't hava a SPEC.md at all --- but maybe that's because it's evolving and we should write the spec for it after it has stoped evolving.

(4) There is no top-level documentation explaining how symatom and symcode responsibilites are split. Some parts of this is indirectly explained in lower level docs, but there is no top-down attempt.

(5) There are three potentially major disruptive changes I have mostly kept quiet about so far (for fear of creating confusion). They have the possibility of cutting the embedding size down a lot.
	(5a) One relates to signs[implemented 0b9f6a2ec69ae35fa95e],
	(5b) I keep forgetting. Ah yes, it is that there are BLAH^2 zips associations, and one you know all but one by class, and provided you know the column and row headings (via single row encoding) you are done. [This also means that all-against-all comparisons don't need recording propvided you have the full row.] See longer description belopw headed "Descrption of 5b in detail".
	(5c) Removes 1/2 of the elements in ORBIT rows as they are +- sign repeats.
	(5d) one relates to pre-row knowledge being subtracted off
Any of them could be discussed.

(6) If there was a demo for symcoder, it is likely that the end-to-end-output could potentially become a unit test BUT this is also a bad idea as there are lots of allowed choices (types of canoncialisation, and optimisations hinted at in (5) above, that would change the output of an end to end test by removing redundancy, so maybe we don't really want to over-fixate on having tests that will break when output becomes less redund



* For any plan (context, groups, etc) the encoding symbol calculations in symatom will only be done once regardless of the numerical values that the encoder will eval. I.e. if we have a thing that can encode 2electrons, 1muon, 4jets then the output eval numbers will vary from event to event, but the structure of the encoding is always the same ... and so one CHUNK of the encoded numbers will be known to have the meaning of, say "the zip of a row with the following FO pair and canonical orbit representative (dot(a,b), dot(a,q)) and so any later DEcoder would be able to unwidn all this knowing just the plan.  But it would be good nonetheless if when a nuimerical encod coding was built a human/machine readable set of metadata was built along side it to record what the meaning of the numbers just output were. As I say, this wordy metadata is technically independent of the  numerical evals themselves, and is just a function of the plan. Thus it could either be generated on-the-fly by the first eval pass (and then skipped  if subsequent events are encoded with the same plan) .... or perhaps the skipping could be omitted until/unless it were fonud to be too costsly. [No point in premature optimisation.] If it were costly, it could be pre-generated by a single pass though with (say) all zeros as input data. Hey, that's cool! :) Yes that would be how to do it.

* The output order of the row-pair-zips is a bit odd (see below, describe_demo.sh) with dot(p,q)xdot(p,q) coming out first (so first operation with last group rather than first operation with first group).  It's not an important fix, but it would be nicer if they came out first-first. 

* Also I am surprised that dot(p,q)xdot(p,q) seems to be a thing, though am pleased that it compressed to just 1 element. How did it get there? It's not a row pair? [But it is worth recording given 5c)


* At some point we should modify describe.py so that describe_demo.sh would produce something like this:
    (venv) CGL-M4-2222:symcoder lester$ ./describe_demo.sh 
    [0:1]  dot(p,q)×dot(p,q)  SS  u=(0,2)  v=(0,2)  shared=(0,2)  len=1
    [1:7]  dot(p,q)×dot(a,p)  SS  u=(0,2)  v=(1,1)  shared=(0,1)  len=6
    [7:10]  dot(p,q)×dot(a,b)  SS  u=(0,2)  v=(2,0)  shared=(0,0)  len=3
    [10:13]  dot(p,q)×eps3(a,p,q)  SA  u=(0,2)  v=(1,2)  shared=(0,2)  len=3
    [13:19]  dot(p,q)×eps3(a,b,p)  SA  u=(0,2)  v=(2,1)  shared=(0,1)  len=6
    [19:20]  dot(p,q)×eps3(a,b,c)  SA  u=(0,2)  v=(3,0)  shared=(0,0)  len=1
    [20:32]  dot(a,p)×dot(b,q)  SS  u=(1,1)  v=(1,1)  shared=(0,0)  len=12
    [32:44]  dot(a,p)×dot(b,p)  SS  u=(1,1)  v=(1,1)  shared=(0,1)  len=12
	...

rather than something like this which it does at present:
    (venv) CGL-M4-2222:symcoder lester$ ./describe_demo.sh 
    [0:1]  dot×dot  SS  u=(0,2)  v=(0,2)  shared=(0,2)  len=1
    [1:7]  dot×dot  SS  u=(0,2)  v=(1,1)  shared=(0,1)  len=6
    [7:10]  dot×dot  SS  u=(0,2)  v=(2,0)  shared=(0,0)  len=3
    [10:13]  dot×eps3  SA  u=(0,2)  v=(1,2)  shared=(0,2)  len=3
    [13:19]  dot×eps3  SA  u=(0,2)  v=(2,1)  shared=(0,1)  len=6
    [19:20]  dot×eps3  SA  u=(0,2)  v=(3,0)  shared=(0,0)  len=1
    [20:32]  dot×dot  SS  u=(1,1)  v=(1,1)  shared=(0,0)  len=12
    [32:44]  dot×dot  SS  u=(1,1)  v=(1,1)  shared=(0,1)  len=12
	...
Probably it would be best to put the decorated operations as an extra column (rather than modifying the first column as I have done above) as the new info is redundant output given to everything earlier on in the output, and nice to keep that earlier output all independent. Also it does not matter what the vectors are called to decode the output, so the extra letters are decoration not fundamental.

## Description of (5b) in detail:
ALTHOUGH ... I have just been thinking about how to explain 5b to you, and I realise that it is going to be easier if the row pairs are explored in a  very parcicular order (if that's easy to arrange). The display output shows that at present the order is this:

[0:1]  dot×dot  SS  u=(0,2)  v=(0,2)  shared=(0,2)  len=1
[1:7]  dot×dot  SS  u=(0,2)  v=(1,1)  shared=(0,1)  len=6
[7:10]  dot×dot  SS  u=(0,2)  v=(2,0)  shared=(0,0)  len=3
[10:13]  dot×eps3  SA  u=(0,2)  v=(1,2)  shared=(0,2)  len=3
[13:19]  dot×eps3  SA  u=(0,2)  v=(2,1)  shared=(0,1)  len=6
[19:20]  dot×eps3  SA  u=(0,2)  v=(3,0)  shared=(0,0)  len=1
[20:32]  dot×dot  SS  u=(1,1)  v=(1,1)  shared=(0,0)  len=12
[32:44]  dot×dot  SS  u=(1,1)  v=(1,1)  shared=(0,1)  len=12
[44:50]  dot×dot  SS  u=(1,1)  v=(1,1)  shared=(1,0)  len=6
[50:56]  dot×dot  SS  u=(1,1)  v=(1,1)  shared=(1,1)  len=6
[56:62]  dot×dot  SS  u=(1,1)  v=(2,0)  shared=(0,0)  len=6
[62:74]  dot×dot  SS  u=(1,1)  v=(2,0)  shared=(1,0)  len=12
[74:86]  dot×eps3  SA  u=(1,1)  v=(1,2)  shared=(0,1)  len=12
[86:92]  dot×eps3  SA  u=(1,1)  v=(1,2)  shared=(1,1)  len=6
[92:98]  dot×eps3  SA  u=(1,1)  v=(2,1)  shared=(0,0)  len=6
[98:104]  dot×eps3  SA  u=(1,1)  v=(2,1)  shared=(0,1)  len=6
[104:116]  dot×eps3  SA  u=(1,1)  v=(2,1)  shared=(1,0)  len=12
[116:128]  dot×eps3  SA  u=(1,1)  v=(2,1)  shared=(1,1)  len=12
[128:134]  dot×eps3  SA  u=(1,1)  v=(3,0)  shared=(1,0)  len=6
[134:140]  dot×dot  SS  u=(2,0)  v=(2,0)  shared=(1,0)  len=6
[140:143]  dot×dot  SS  u=(2,0)  v=(2,0)  shared=(2,0)  len=3
[143:146]  dot×eps3  SA  u=(2,0)  v=(1,2)  shared=(0,0)  len=3
[146:152]  dot×eps3  SA  u=(2,0)  v=(1,2)  shared=(1,0)  len=6
[152:164]  dot×eps3  SA  u=(2,0)  v=(2,1)  shared=(1,0)  len=12
[164:170]  dot×eps3  SA  u=(2,0)  v=(2,1)  shared=(2,0)  len=6
[170:173]  dot×eps3  SA  u=(2,0)  v=(3,0)  shared=(2,0)  len=3
[173:179]  eps3×eps3  AA  u=(1,2)  v=(1,2)  shared=(0,2)  len=6
[179:182]  eps3×eps3  AA  u=(1,2)  v=(1,2)  shared=(1,2)  len=3
[182:188]  eps3×eps3  AA  u=(1,2)  v=(2,1)  shared=(0,1)  len=6
[188:200]  eps3×eps3  AA  u=(1,2)  v=(2,1)  shared=(1,1)  len=12
[200:203]  eps3×eps3  AA  u=(1,2)  v=(3,0)  shared=(1,0)  len=3
[203:215]  eps3×eps3  AA  u=(2,1)  v=(2,1)  shared=(1,0)  len=12
[215:227]  eps3×eps3  AA  u=(2,1)  v=(2,1)  shared=(1,1)  len=12
[227:233]  eps3×eps3  AA  u=(2,1)  v=(2,1)  shared=(2,0)  len=6
[233:239]  eps3×eps3  AA  u=(2,1)  v=(2,1)  shared=(2,1)  len=6
[239:245]  eps3×eps3  AA  u=(2,1)  v=(3,0)  shared=(2,0)  len=6
[245:246]  eps3×eps3  AA  u=(3,0)  v=(3,0)  shared=(3,0)  len=1


which appears to go lex sort order by u, then v then shared with some other priority that puts dot before eps3 and presumably if we had multiple ops with identical rank (which we don't have in this example) then they would also be ordered somehow.

Is that really the order.

To explain 5b  efficiently it would be convenient if it is guaranteed :

* for a fixed opxop in column 2,
* and for a fixed u and v in columns 4 and 5 of the display output

That all shared values for those overlaps came consecutively in an unbroken block. I.e. when all shared combinations of u=(3,4) and v=(4,1) for a given operator pairing WILL COME TOGETHER (not interrupted by, say, changes of operation or changes of u or v).

Is that an order we are currently guaranteed?

-- YES
Yes, this is definitively guaranteed. The sort key is:

key=lambda pf: (
    _op_fl_key(pf.op_u, pf.flavour_u),   # (name, rank, parity, sym, counts)
    _op_fl_key(pf.op_v, pf.flavour_v),
    pf.overlap,
)
overlap is the third component of the key, so for any fixed (op_u, flavour_u, op_v, flavour_v) the first two components are identical and Python's sorted is stable — all overlap tuples for that combination land in one unbroken block, sorted lexicographically within it.

The only subtlety: PairFlavour.__init__ may swap the sides so the "smaller" side (by _op_fl_key) is always op_u. That's why you see dot×eps3 everywhere rather than eps3×dot. But after that normalization the three-level sort is strict: op_u key, then op_v key, then overlap.

So the structure you want for 5b — enumerate over (op×op, flavour_u, flavour_v) blocks, and within each block see all overlaps consecutively — is an exact description of the current output order. Ready to hear 5b.




OK - let's call one of those consecutive blocks an "OVERLAP BLOCK" as it covers all forms of overlap between those two operators each having a fixed flavour. And call each member of an OVLEAP BLOCK an "ASSOCIATION" -- this name because it tells us which atoms from the top row area ASSOCIATED with which other atoms from the bottom row, in the sense that they "move together" under the group action.

Note also that if the pair is (u,v) then G.u is a set (often a small one) but the orbit of the pair, G.(u,v), is often longer than G.u because the top elements of G.(u,v) [I'll call this. the TOP ROW] could contain (say) 4 copies of every u=a.b if each a.b is associated with four things v. So theTOP ROW  is a mulitiset (with a pre-computable duplucation factor) wqhile G.u is just a set and lacks that factor.

Preamble out of the way, we could imagine unpacking all the ASSOCIATIONS in one OVERLAP BLOCK into an ASSOCTIATION TABLE, where we label the columms with the elements of G.x (a set, no dupliucation) and we label the rows with the elements of the set G.y (no duplications in there either).  Initially the table is empty, but when we read our first ASSOCIATION from the OVERLAP BLOCK we could insert a number "1" into each cell were the appropruate elements are assocaited. Then we could move to the 2nd ASSCITATION in the ASSOCIATION block, and put a 2 in every cell of the TABLE where those ele,memts are assopciated. None of these will overlap with an exising 1. {check out can see why]

Keep going.

By the end of the OVERLAP BLOCK every cell of the table will be filled by the integer telling us whether that association was made by the 1st, or the 2nd or the 3rd (etc) ASSOCIATION within the block.

Hopefully you can see now, that

IF:

we already  know the atoms on the column lables (G.x) (on their own, not in any association), and have a way of storing their evals in a permutation invariant way,   and  if we already know the atoms which are the row labels (on their own, not in any association) and can store their evals in a perm invariant way,

THEN:

It is not necessary to do encode ONE of the associtions in the OVERLAP BLOCK .... since once you know where ALL but one of the associations fill the table, you will know by process of elimination where the last one goes!

The last set of associations could be reconstructed, information theoretically, form the column headings, row headings, and the first N-1 associations in the block.

And you don't have to choose the association wihtin the block to drop to be the last one. You can choose it to be THE WOST ONE -- i.e, the one that yields the largest number of coefficients.

One corollary of this is if there is only one ASSOCIATION in the BLOCK, then the BLOCK does not need storing at all (provided that the row headings and column headings have been recorded).

Another amazing thing about this is that very often in large OVERLAP_BLOCKS one of the ASSOCIATIONS is much larger than the others, so there is often a huge saving here. Indeed The cost of saving the heading (which are one dimensional, growing like the larger of the table width and height) is always better than saving the worst table ASSOCIATION (which grows as the table's are).

Finally, how would one store the row (or column) headings in a perm invariant way?

Well, one COULD use a 1-D (k=1) polynomial encoder over the reals which would put the values into an identical number of reals. Or you could (my preference) just sort the values into numerical order.  BUT this choice should be selectable for reasons I have not told you about yet.

BUT BUT BUT (and this is the final cherry)  the storage of that G.u or G.v row/col header should NOT be done when processing this OVERLAP BLOCK because almost certainly there is another OVERLAP BLOCK in the encoding that will re-use the same G.u (possibluy this time as a row or a column header) and the same G.v.  Since the G.u and the G.v are just USED BY BUT NOT PART OF the association, the right place to store them all is BEFORE any OVERLAP BLOCKS are processed.

I.e. phase one of storage is

FIRST store (either) the polynomial encoded (or) sort-encoded BUT ALWAYS non-associated orbits of each canonically different element of repL FIRST , and then with this new non-association coding block out of the way [self documenting with its desribed]

THEN: each OVERLAP BLOCK can be investigated for the worst-sized ASSOCIATION, and then all associations EXCEPT that one would be encoded.

Phew! And this may perhaps need a latex level document with nice tables to show what's happening unless a good clear way of explaining it in asciii is possilbe.


* Contemplate removing symcoder/SPEC.md --- maybe no-one needs it.
