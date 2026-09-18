# Stochastic Simulation

Sampling and Monte Carlo methods implemented for the graduate course in Stochastic
Simulation, UNAM (MSc Mathematical Sciences, 2025). Each part pairs a derivation with an
implementation and a numerical experiment.

## Contents

**1. Sampling from log-concave densities**

- Rejection sampler built on a fixed exponential envelope, with expected cost derived in closed form
- Adaptive rejection sampling: piecewise-linear tangent hull on `-log f`, refined on every rejection
- Runtime comparison over `N(0,1)`, shifted Gamma and shifted Beta, 10^5–10^7 samples

**2. Bias and multilevel Monte Carlo for a random recursion**

- Existence of a stationary distribution for `X_n = X_{n-1} U_n + S_n`, via convergence of the backward series
- Bias of `X_n` as an estimator of the stationary mean, and the cost of reaching MSE `eps^2`
- An MLMC estimator for the same target, with its complexity compared against the single-level cost

**3. Gibbs sampler for a Gaussian mixture model**

- Full conditionals derived for allocations, component parameters and mixture weights
- Gibbs sampler over the posterior, with burn-in and configurable hyperparameters
- Posterior mean density estimated and fitted to the `galaxy` dataset

**4. ULA vs MALA on a bimodal target**

- Both samplers implemented for `V(x) = x^4/4 - 9x^2/2`
- ULA diverges at every step size tried, from `h = 1` down to `h = 1e-5`
- MALA stays bounded and recovers the target; a large step lowers acceptance but never explodes

## Notes

Coursework for Stochastic Simulation, UNAM Faculty of Sciences — instructors María Fernanda
Gil Leyva Villa and Jorge González Cázares.

## References

- Gilks & Wild (1992), *Adaptive Rejection Sampling for Gibbs Sampling*
- Giles (2008), *Multilevel Monte Carlo Path Simulation*
- Roberts & Tweedie (1996), *Exponential convergence of Langevin distributions and their discrete approximations*
