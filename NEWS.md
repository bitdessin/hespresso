# hespresso 1.1.0

## NEW FEATURES

* Added the HOBIT MLE implementation as `hobit_mle()`.
* Added HOBIT MCMC (strict), available through `hobit_mcmc(mode = "strict")`.

## SIGNIFICANT USER-VISIBLE CHANGES

* Renamed `ExpMX` to `SeqCountData`.
* Renamed the original HOBIT implementation to HOBIT MCMC (fast),
  which available through `hobit_mcmc(mode = "fast")`.
* Changed `hobit()` to use HOBIT MLE (`hobit_mle()`) by default.
* Moved `cmdstanr` from Imports to Suggests.
* Removed the exported `norm_counts()` function. HOBIT MCMC (fast)
  now performs the same normalization internally.
* Redesigned the artificial RNA-seq count simulator to improve the realism
  of benchmarking datasets and the definition of ground-truth homeolog
  expression ratio shifts.
* Changed argument name `seed_expmx` to `seed_counts` in the `sim_homeolog_counts()` function.


# hespresso 1.0.6

## BUG FIXES

* Fixed bug in group name handling during dispersion estimation.


# hespresso 1.0.5

## BUG FIXES

* Fixed formatting bug in HomeoRoq output.


# hespresso 1.0.0

## NEW FEATURES

* Release of the package.

