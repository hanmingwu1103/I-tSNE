# itsne 0.2.0

First CRAN-targeted release.

## Breaking changes

* `run_simulation_studies()` and `run_real_data_analysis()` no longer write into
  the working directory. `output_dir` now defaults to a session-specific
  subdirectory of `tempdir()` (`itsne-simulations` and `itsne-real-data`
  respectively). Supply an explicit path to write elsewhere. This brings the
  package into line with CRAN policy; code that relied on the previous
  `"results/simulations"` and `"results/real-data"` defaults must now pass the
  path explicitly.

## New features

* `itsne_mm()` gains a post-fit endpoint-order diagnostic. The returned list now
  has a third component, `violation_summary`, reporting `n_violations`,
  `n_coordinates`, `prop_violations`, `n_objects_violating`, `max_violation`,
  and `repaired`. A warning is emitted when at least one reduced-space
  coordinate violates `lower <= upper`.

  **The optimized endpoints are never modified.** No reordering, swapping,
  sorting, clipping, coordinatewise hull, or `pmin()`/`pmax()` repair is applied;
  `embedding` always contains the raw optimizer output. The order penalty is
  soft, so for any finite `penalty_lambda` a local numerical solution is not
  guaranteed to be order-valid, and this diagnostic makes that visible rather
  than silently hiding it. The no-repair contract is pinned by dedicated tests.

* Added a test suite covering interval-data validation, standardization, vertex
  and quantile expansion, expanded-cloud perplexity admissibility, the MM
  no-repair contract, CR positive radii, seed determinism, distance matrices,
  LCMC, plotting with one, two, and many groups, workflow output confined to
  temporary directories, and example-data integrity.

## Documentation

* `itsne_mm()` now documents all three returned components, and its examples show
  how to inspect `violation_summary`.
* Added `inst/extdata/README.md` recording the source, citation, licence,
  derivation, column meaning, and SHA-256 checksum of each bundled dataset, and
  `data-raw/prepare_digits_interval.R` documenting the Digits derivation.
* Added `inst/CITATION`.
* Expanded the `Description` field to state what the package does, distinguish
  the expansion-based and direct variants, and cite the underlying methods.
* Examples that write files now write only to `tempdir()`.

## Infrastructure

* Added a GitHub Actions `R-CMD-check` matrix over Ubuntu (release, oldrel-1,
  devel), Windows (release), and macOS (release).
* Declared `grDevices`, `stats`, and `utils` in `Imports`, matching actual use.
* Added `Language: en-US`.

# itsne 0.1.0

* Initial public GitHub release (2026-04-08). Not submitted to CRAN.
* Interval-valued t-SNE variants `itsne_vm()`, `itsne_qm()`, `itsne_mm()`, and
  `itsne_cr()`.
* Interval PCA baselines `ipca_vm()`, `ipca_qm()`, and `ipca_cr()`.
* Interval-aware LCMC evaluation utilities.
* Plotting helpers and reproducible analysis workflows.
