# dual_attacks


# [DP23]


**Lemma 8.** Let Λ be a lattice of volume 1, and r > 0 such that r < λ1(Λ)
2 GH(n) .
Then, for a target t uniform in R^n/Λ, it holds with probability rn that t is at
distance at most r GH(n) from the lattice.

Let us use this lemma in the case of a random lattice of volume 1, or more
specifically, one where we expect λ1(Λ) ≈ GH(n). Using Lemma 8 with r = 0.49,
the probability that a uniform target lies in the ball of radius r GH(n) equals
r^n. When taking T = 2.05^n uniform samples, on expectation we have T · r^n >
1.004^n >> 1 of the uniform samples to fall in this ball and therefore with high 
probability there will be one such target at most r GH(n) away from a lattice
point. More concretely, the probability that none of these targets lies in this ball
is (1 - r^n)^T → e^(-1.004^n)
(as n → ∞), so with overwhelming probability there is
a uniform target in the ball of radius r GH(n).
On the other hand, recall that the actual BDD target had an expected length
of σ√n ≈ 0.89 GH(n) > r GH(n). We note that 0.89 > r, so we expect one
uniformly sampled candidate lying closer to the lattice than the solution we are
looking for. However, the score function f_W is precisely meant to associate larger
score to closer targets, so we expect that this uniform sample will get a higher
score than the BDD sample and thus, the algorithm gives with overwheling
probability a wrong result.


**Improvement.**  How to get closer to the 0.89 coefficient ? (instead of 0.49). Use the 
following heuristic claim instead of the Lemma 8.

**Heuristic Claim 4.** Let Λ be a random lattice of volume 1, and r ∈ (0, 1).
For a target t sampled uniformly in the torus R^n/Λ, we have a probability of
r^n · (1 − n^{O(1)} r^n) that t is at distance at most r GH(n) from the lattice.

**The Contradictory Regime.** The idea now is that we end up at a contradiction whenever
we instantiate the claim from above with the smallest possible r such that it is
likely there is such a point among the T uniform samples.
For a given number of random samples T , we will pick σ for the BDD sample
as before, i.e. such that ln T = N ε^2 where ε = exp(−2π2σ2`2). By the above
Heuristic Claim 4, among those T targets, with constant probability there exists
a target at distance at most r · GH(n) from the lattice, where r = 1/T^(1/n) . When
for a given T , the length r GH(n) is smaller than the expected length of the
BDD sample, i.e. √n · σ, this contradicts Heuristic Claim 3 and we say it is
in the contradictory regime. This concrete contradictory regime is depicted in
Figure 2.




# [SP24]

**6. Comparison with [23]’s Contradictory Regime**

Their reasoning is as follows: assume that we have a large number of random
candidates (the t^{(i)}_{unif}) and one point close to the lattice L (the point e), 
then Heuristic 2 says that we can always distinguish e from the candidates (since
it has maximum value of g_W ). The contradiction comes from the fact that in
reality, for T large enough, many of candidates will be closer to L than e and
therefore no algorithm can distinguish them [21]. This gives rise to what [23]
calls the “contradictory regime” where an algorithm would somehow be able to
distinguish indistinguishable distributions.


**6.2 On the distribution of targets**
The authors of [23] decided to modularize the algorithm by separating the lattice
in which dual-distinguishing is done (Lq(A_dual)) from the part of the lattice
that is enumerated over (Lq(A_guess)). In fact, Heuristic 1 only mentions the
dual-distinguishing and not the enumeration. This however, poses a difficulty
because it is clear that the “targets” (b − A_guess \tilde{s}_guess) in our terminology, (t^(i)_
unif in Heuristic 1) are not arbitrary but have some structure.
The authors of [23] decided to model the statistics of the targets in a way that
is independent of the actual choice of A_guess: they chose the uniform distribution
over the fundamental domain of Lq(A_dual). In the case of [37] and our algorithm,
the algorithm exclusively works over integers which is why we propose Heuristic 2
as an integer-version of Heuristic 1. This means that we now have two different
settings:
– In Heuristic 2, t^(i)_unif is sampled uniformly in Z^m/L.
– In reality, t^(i)_unif = e + x^(i) where x^(i) can be any vector in L'\qZ^m where L'
is another random q-ary lattice, chosen independently of L. In our algorithm,
L = Lq(A_dual) and L' = Lq(A_guess).
Indeed, a key point in the proof of Theorem 6 is to show that points of the form
e + x^(i) as described are always far away from L, a fact that does not hold for
completely uniform targets. As a result, with high probability over the choice of
A, the targets (except for the correct guess) are all bounded away from 0 in the
dual lattice. For uniform targets, the argument of [23] is statistical in nature:
while there can be very short vectors, they are unlikely and the contradiction
comes from the fact that if we try too many targets, we will eventually find
a short one and get a false-positive. On the other hand, our algorithm and
analysis is not statistical: for the vast majority of choices of A, all targets satisfy
the bound unconditionally and we can safely look at all targets without the risk
of any false-positive.



et W = (w1, ... , wN ) ∈ L^⊥_q(A_dual)^N be sampled
according to D_{L^⊥_q(A_dual),qs}.






# [DP24]

To be able to work with the output of a lattice sieve in our analysis, we will
make a heuristic assumption about the distribution of the lattice vectors that a
lattice sieve algorithm outputs.

**Heuristic 2.** Given a unit-volume lattice Λ, the output distribution of a lattice
sieve with a saturation radius of rsat ≥ 1 and a saturation ratio of fsat ∈ (0, 1],
is a list of vectors W ⊂ R^n of size N = 1/2 f_sat r^n_sat, where its elements are inde-
pendently sampled uniformly at random from the ball of radius r_sat GH(n).


This heuristic makes two simplifications on the output of a lattice sieve.
The first simplification made in the heuristic is that not much changes when
going from a list of N vectors that are the output of a sieve, to a list of N
vectors that are each sampled uniformly from Λ ∩ r_sat GH(n) B^n. Although the
heuristic allows for duplicates, many distinct values will be sampled, which gives
a motivation for this simplification. As an illustration, consider the following:
sampling n times uniformly from a set of size n yields in expectation a set of
size n(1 − 1/e) ≈ 0.63n as n → ∞.
The second simplification is that the lattice structure W ⊂ Λ is ignored in
the output. The Gaussian Heuristic predicts that the norm of the lattice vectors
follows a similar distribution as the norm distribution of points uniformly from
the ball. When running a sieve on a random lattice, this assumption seems fair,
because we do not expect the lattice to be distorted in any particular direction.
Note that this heuristic does not assume all vectors are of the same length
r_sat √n, as done in e.g. [GJ21]. Instead the heuristic predicts there may be some
shorter vectors w, albeit with smaller probability, but still most of the vectors
are close to the boundary of the ball. The very short vectors are more beneficial
in the dual attack than longer vectors, so the analysis will be more conservative
when taking shorter vectors into account.