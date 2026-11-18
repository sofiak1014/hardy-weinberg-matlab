# Hardy-Weinberg Equilibrium Explorer (MATLAB)

A MATLAB tool for exploring the Hardy-Weinberg principle: simulate random mating,
sample genotypes, and test whether data fit equilibrium.

## Files
| File | What it does |
|---|---|
| hw_expected.m | Expected genotype frequencies for allele frequency p (and inbreeding coefficient F) |
| hw_random_mating.m | Simulates one generation of random mating |
| hw_sample.m | Draws a random sample of n individuals and counts genotypes |
| hw_chisq.m | Chi-square test for Hardy-Weinberg equilibrium |
| hw_tests.m | Verification script. Run it, and every line should print PASS |

## How to run
Put all files in one folder, open MATLAB, and run `hw_tests`.
Base MATLAB only (no toolboxes required).

## Status
In progress. Figures and GUI coming.
