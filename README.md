# Stochastic Simulation

Implementations of sampling and Monte Carlo methods, written for the graduate course in
Stochastic Simulation at UNAM (MSc Mathematical Sciences, 2025). Four self-contained
studies: adaptive rejection sampling, multilevel Monte Carlo for a random recursion, a Gibbs
sampler for Bayesian mixture models, and a comparison of discretized Langevin samplers.

Each part pairs a derivation with an implementation and a numerical experiment that tests
whether the theory holds in practice.

---

## 1. Sampling from log-concave densities

A density `f` is log-concave when `-log f` is convex. Two samplers are implemented and compared.

**Fixed envelope.** For a log-concave `f` with mode at 0, the exponential bound

```
f(x) <= f(0) * min{1, exp(1 - x*f(0)/(1-q))}   for x >= 0,   q = F(0)
f(x) <= f(0) * min{1, exp(1 + x*f(0)/q)}       for x <  0
```

gives a rejection sampler whose expected cost is computed in closed form.

**Adaptive rejection sampling.** Because `-log f` is convex, a piecewise-linear function
tangent to it at a finite set of points lies below it, so `exp(-g_n)` dominates `f`. Every
rejected proposal is added to the tangent set, the hull tightens, and the acceptance
probability tends to 1. The envelope improves at the cost of an increasingly expensive
proposal step.

**Experiment.** Both samplers draw 10^5–10^7 variates from `N(0,1)`, shifted
`Gamma(a,1)` and shifted `Beta(a,b)` across several parameter settings, with wall-clock
times compared.

**Finding.** [Which method wins, and where the crossover sits — the adaptive method only
pays off once the fixed envelope's acceptance rate is low enough to offset the cost of
maintaining the hull.]

## 2. Bias and multilevel Monte Carlo for a random recursion

For the chain `X_n = X_{n-1} U_n + S_n` with iid `(S_n, U_n)`:

- A stationary distribution is shown to exist by rewriting `X_n` in backward form and
  proving the resulting series converges absolutely, almost surely and in `L^1`.
- The bias of using `X_n` as an approximation to the stationary mean is computed, giving
  the cost of reaching mean squared error `eps^2` with a plain MCMC estimator.
- An MLMC estimator is constructed for the same target, with its complexity derived and
  compared against the single-level cost.

This is the one part with a direct line to derivatives pricing: the bias–cost tradeoff here
is the same one that governs discretization error in path-dependent Monte Carlo.

## 3. Gibbs sampler for a Gaussian mixture model

A Bayesian mixture `f(x) = sum_j w_j N(x | mu_j, 1/tau_j)` with Normal-Gamma priors on the
component parameters and a Dirichlet prior on the weights.

- Full conditionals are derived for the allocation variables `d_i` (categorical), the
  component parameters `theta_j` (Normal-Gamma), and the weights `w` (Dirichlet).
- The Gibbs sampler takes the data, burn-in length `B`, post-burn-in length `T` and the
  hyperparameters, and returns the chain.
- The posterior mean density is estimated by averaging the sampled mixtures, and is shown
  to be a density itself.

**Experiment.** The `galaxy` dataset, with `m = 7` components, `T = 10000` iterations after
`B = 2000` burn-in. Output is the estimated density plotted over the data histogram.

## 4. ULA explodes, MALA does not

Target `pi(x) ∝ exp(-V(x))` with `V(x) = x^4/4 - 9x^2/2`: bimodal, modes near `±2.12`,
barrier at 0. The associated Langevin diffusion is ergodic with `pi` as its unique invariant
law — but discretizing it is not automatically safe.

- **ULA** (unadjusted): `X_{n+1} = X_n - h ∇V(X_n) + sqrt(2h) Z`. For *every* step size
  tried — `h = 1` down to `h = 1e-5` — the trajectory eventually diverges to `±∞`. Sometimes
  within a hundred steps, sometimes after thousands. It never converges to `pi`.
- **MALA** (Metropolis-adjusted): the same proposal with an accept/reject step. The chain
  stays bounded even at `h = 0.2`–`0.3`, visits both modes, and its empirical histogram
  matches `pi`. Too large a step only drives the acceptance rate down; it never explodes.

The point is that the cubic drift makes the Euler step unstable at large `|x|` regardless of
`h`, and the Metropolis correction is what rescues it. Trajectories and terminal histograms
are plotted for each `(h, algorithm)` pair.

---

## Running

```
[pip install -r requirements.txt]
[python exercise1.py]
```

[One line per script: what it produces and roughly how long it takes.]

## Notes

Coursework for Stochastic Simulation, UNAM Faculty of Sciences — instructors María Fernanda
Gil Leyva Villa and Jorge González Cázares. The problem set is included as
[`Tarea_Simulacion_Estocastica.pdf`](.).

## References

- Gilks & Wild (1992), *Adaptive Rejection Sampling for Gibbs Sampling*.
- Giles (2008), *Multilevel Monte Carlo Path Simulation*.
- Roberts & Tweedie (1996), *Exponential convergence of Langevin distributions and their
  discrete approximations*.
