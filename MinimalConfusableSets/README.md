Run compare_M7k3.sh to compare two different match matrices.

OK - I have decided to bite the bullet and write some code (I suspect you will contribute to it too) to do what we talked about yesterday. I.e.:

	(1) enumerate all matches (and multi-matches) in a sensible order, then

	(2) generate L matrices from them, then

	(3) [[optionally]] canonicalise or reduce the L's to some common forms, e.g. perhaps by putting into RRE form and sorting columns where possible, maybe even hashing them) before,

	(4) passing these L’s to a decider function “f” (or one of a set of related decider functions f) whose job is to work out whether the null space requires (or does not require) any of the alphas, betas, gammas (etc) to be zero.

The decider function(s) will probably need to use floats and may wish to imagine concrete distributions of vectors on a Fibbonaci sphere. Or they might work other ways. I propose that there be a library of different “f”s, using different assumptions or methods, to allow them to be checked / validated against each other. E.g. some might be reliable but slow and only able to work in low n,k … while others might perhaps purport to work for large n,k but maybe are more complex etc.

The common rule is that all f’s must expect AT A MNIMUM to be given as standard arguments:

	(i) an integer k (the number of dimensions) and
	(ii) a SymPy L matrix of rationals (potentially this could even be integers as we could in principle multiply through by the LCM of the denominators …. but let’s just assume rationals for the moment).

Individual f’s may take extra arguments that allow them to function, adding to those above. E.g. there could be an f which takes as a third argument a set of actual locations for vectors on the fibbonaci sphere — and if this f was used you would always know where it had put the vectors. Whereas there could be a different f which decides for itself where to put vectos on the fibonaci sphere (perhaps it even does so multiple times) — so it would not take such a list of vectors as one of its inputs.

Individual f’s may also declare extra pre-conditions that they anticipate will have been fulfilled in the supplied L’s. E.g. an “f” implementation that REQUIRES its input L’s to be in RRE Form, say, would report this in a function property — whereas other “f”s that don’t have such a requirement would not set that property.

I propose to write the following within my Echinocoder githib repository.
https://github.com/kesterlester/Echinocoder

You can therefore contribute by forking that and then issuing pull requests, etc. However you could also develop code separately and send it to me if you preferred (though this may be less efficient in the long run).

I plan to start on (1) and (2) first and will not start (4) until I am happy with (1)+(2).  So if you were able to make a start at one or more “f”s for (4) then this parallel work could allow us to put things together efficiently.

You may suppose that the L matrices being fed to any “f” you write would look something like this:

```python
from sympy import Matrix, Rational
# Define a matrix with Rational entries
L = Matrix([
   [Rational(1, 1), Rational(2,3), Rational(1, 3), Rational(3,1)],
   [Rational(0, 1), Rational(0,1), Rational(1, 1), Rational(2,1)]
])
```

# Unit Testing:
NEWSFLASH: some global vars temporarily block the running of pytest framework other than by the workaround script run_tests.sh. So for the moment run pytest via that script. Otherwise, in the past (and in the future):

Unit tests are organised for pytest, so run with
```
pytest
```
or if you don't want stdout suppressed use
```
pytest -s
```

# Notation:
```python
We work with (multu)sets of $n$ vectors in $k$-dimensions.

Encoding dots are done with $m$ directions: $D=\{\vec d_1, \vec d_2, ... , \vec d_m\}$

There are $M$ vectors ($M=ceil(m/(k-1)$) which are orthogonal to the vectors in $D$ (these are the $M$ bad-bats): $B=\{\vec q_1, \vec q_2, \cdots, \vec q_M\}$. E.g. $\vec q_1 \cdot m_i = 0$ for all $i\in \{ 1, 2, \cdots, k-1 \}$.

The (multi)sets $E$ and $O$ have even and odd sums of the elements of $B$:

$E=\{ \vec 0, \vec q_1 + \vec q_2, \cdots \}$
$O=\{  \vec q_1, \vec q_2, \cdots, \vec q_2 + \vec q_5 + \vec q_9, \cdots \}$
```
# Integration Lester -- Ruane

```
from ruane import factory_something

checker_func = factory_something(M=M, k=k, ...)

# checker_func will be a function that can be called with a sympy matrix L and other things specified by Tom but which already "knows" it will be used with M and k so that are not supplied.

# Lester will use
for L in some_matrices_L(M=M, k=k):
	tom_answer = checker_func(L, ... other switches) # in tom_answer ideally somewherte is a statement about whether (a) there is collapse to "nothing" (or not) at the alpha,beta,gamma level, (b) [for future] at the L -integer level one can have inconsitent L , so some info on this might come out.
	do_stuf_with(tom_answer) 
```

