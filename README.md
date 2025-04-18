# dual_attacks



**Lemma 8.** Let Λ be a lattice of volume 1, and r > 0 such that r < λ1(Λ)
2 GH(n) .
Then, for a target t uniform in Rn/Λ, it holds with probability rn that t is at
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






