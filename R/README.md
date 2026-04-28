# R source files

- `bs2_core.R`: statistical core: BS2 simulation, log-likelihood, profiled ML, modified moments, conditional pseudo-likelihood, Fisher information, LR/score/Wald/MM-Wald tests, and bootstrap helpers.
- `sim_designs.R`: design grids for the estimation, size, power, bootstrap and GOF/influence experiments.
- `sim_parallel_utils.R`: command-line parsing, deterministic task seeds, checkpointed RDS task construction, `future` backend setup and numerical-thread pinning.

The public bundle keeps only the current score-corrected implementation. Obsolete patch files and old script copies were intentionally excluded to avoid ambiguity.
