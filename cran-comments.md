## Update

This is smartsnp 1.2.1, an update to CRAN version 1.2.0.

* Preserve matrix dimensions when filtering SNPs in the four analysis functions,
  fixing errors when projecting a single sample (PR #13).
* Add regression coverage for single-sample projections with SNP filtering.
* Restore usage signatures in four Rd files to address the current CRAN notes.
* Replace the broken EIGENSOFT URL with its official GitHub repository.

## Local validation

Ubuntu 24.04, R 4.5.3 (conda-forge), R CMD check --as-cran:
0 ERRORs, 0 WARNINGs, 2 NOTEs.

* unable to verify current time (environment/network restriction).
* non-portable compilation flag -march=nocona (conda toolchain default).

Examples, regression tests, vignettes, and PDF/HTML manuals passed.
Remote incoming checks were disabled for this final local run following network
failures in the preceding full run. The preceding run identified a broken
EIGENSOFT URL, which has been corrected. PDF checks used R_RD4PDF=times,hyper
because the local TeX installation lacks the inconsolata font.

## Pending before submission

Run the prepared GitHub Actions matrix on current R-release (Linux) and R-devel
(Windows), with incoming checks enabled. Build the submitted archive with the
current R release. Record those results here before submitting to CRAN.
