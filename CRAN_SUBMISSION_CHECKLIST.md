**CRAN Submission Checklist For `itsne` 0.1.0**

1. Run `R CMD build I-tSNE`.
2. Run `R CMD check --as-cran itsne_0.1.0.tar.gz` from an ASCII-only working path.
3. Confirm the check result has no ERRORs or WARNINGs.
4. Submit `itsne_0.1.0.tar.gz` to CRAN.
5. Include `cran-comments.md` with the submission.
6. After acceptance, push the release commit and tag to GitHub.
7. Update the installation instructions to include the CRAN install path as the primary option.

**Current Local Status**

- Package version: `0.1.0`
- Source package built successfully
- `R CMD check --as-cran` passed from an ASCII-only path with only two benign NOTEs:
  - `New submission`
  - `unable to verify current time`
- CRAN-facing package surface no longer depends on `mdsOpt`
