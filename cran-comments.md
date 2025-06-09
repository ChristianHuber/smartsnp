## Resubmission

This is an update to the smartsnp package (version 1.2.0).

Changes since the previous CRAN version:

* Replaced all uses of the deprecated `vegan::adonis()` with `vegan::adonis2()`.
* Explicitly set `by = "terms"` in `adonis2()` to preserve prior behavior.
* Updated example and test code to ensure compatibility with `vegan` ≥ 2.6.6.
* Fixed a bug in PCA projections involving single-individual input.
* Cleaned up documentation and Rd files to comply with CRAN policies.

## Test environments

* Local macOS, R 4.4.0
* Win-builder (R-devel and R-release)
* Fedora Linux (devel)

## R CMD check results

I have checked the package using `R CMD check --as-cran` and the current version of R-devel on win-builder (https://win-builder.r-project.org).

There were no ERRORs or WARNINGs.

One NOTE:
* "unable to verify current time" — this appears to be related to local system configuration and does not affect package functionality.