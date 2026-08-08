# cran-comments

## Submission

`itsne` 0.2.0 — new submission to CRAN.

The package provides nonlinear dimension reduction and visualization for
interval-valued observations (symbolic data analysis), including four
extensions of t-SNE, interval PCA baselines, interval-aware neighborhood
evaluation, plotting helpers, and reproducible analysis workflows.

This package has never been on CRAN. Version 0.1.0 existed only as a GitHub
release and was never submitted here.

## Test environments

* local Windows 10 x64 (build 19045), R 4.6.0 — `R CMD check --as-cran` on the
  built tarball
* GitHub Actions:
  * ubuntu-latest, R devel
  * ubuntu-latest, R release
  * ubuntu-latest, R oldrel-1
  * windows-latest, R release
  * macos-latest, R release

## R CMD check results

0 ERRORs | 0 WARNINGs | 1 NOTE

```
* checking CRAN incoming feasibility ... NOTE
Maintainer: 'Han-Ming Wu <hanmingwu1103@gmail.com>'

New submission
```

This is the standard note for a first submission. It is the only note; there are
no other notes, warnings, or errors.

## Additional notes for the reviewer

* **Bundled data.** `inst/extdata/` contains three small CSV files used by the
  examples, tests, and vignette so that all checks run offline and
  deterministically. Source, citation, licence, derivation, column meanings, and
  SHA-256 checksums are documented in `inst/extdata/README.md`. The Digits files
  derive from the UCI *Optical Recognition of Handwritten Digits* data
  (CC BY 4.0, attributed); `data-raw/prepare_digits_interval.R` documents the
  derivation. The Face data is a published symbolic-data benchmark distributed
  via the `dataSDA` package, of which the maintainer of this package is also the
  maintainer.

* **No writes outside `tempdir()`.** `run_simulation_studies()` and
  `run_real_data_analysis()` default `output_dir` to a subdirectory of
  `tempdir()`. Nothing is written to the user's working directory, home
  directory, or the package library. A test asserts this.

* **`\dontrun{}` usage.** The examples for those two workflow functions are
  wrapped in `\dontrun{}` rather than `\donttest{}`. They run the complete
  manuscript-scale analyses (eight simulation scenarios and two real datasets,
  including full vertex expansion) and take far longer than is acceptable in an
  automated check. Every other exported function has a runnable example.

* **`skip_on_cran()` usage.** One test is skipped on CRAN: it executes the full
  simulation workflow to confirm that output is confined to the supplied
  directory. The cheap invariant it protects — that both defaults resolve under
  `tempdir()` — is asserted by a separate test that always runs.

* **`itsne_mm()` warning.** `itsne_mm()` deliberately emits a warning when the
  fitted embedding contains endpoint-order violations. The method controls order
  with a soft penalty, so a finite penalty weight cannot guarantee validity; the
  package reports violations through `$violation_summary` and never silently
  repairs the optimized endpoints. Examples and tests use `suppressWarnings()`
  where a violating fit is intentionally exercised.
