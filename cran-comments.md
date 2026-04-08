## Test environments

- Windows 10 x64, R 4.5.0

## R CMD check results

I ran:

- `R CMD build I-tSNE`
- `R CMD check --as-cran itsne_0.1.0.tar.gz`

The final CRAN-style check was rerun from an ASCII-only temporary path on
Windows because the working project path contains non-ASCII characters. Under
that ASCII-only path, `R CMD check --as-cran` completed successfully with two
NOTEs:

1. `New submission`
2. `unable to verify current time`

There were no WARNINGs or ERRORs in the final ASCII-path check.

## Notes

The package provides interval-valued t-SNE methods, interval-aware evaluation
utilities, and reproducible analysis workflows associated with the accompanying
manuscript.
