**Exact CRAN Upload Package / Files Checklist For `itsne` 0.1.0**

**Primary file to upload to CRAN**

Upload this source package file:

- [itsne_0.1.0.tar.gz](/d:/08-MyProjects/00-MasterStudent/2024-吳玟樺/LaTeX/Article/Interval-tSNE_Wu&Wu_20260402/itsne_0.1.0.tar.gz)

This is the file that should be submitted through the CRAN web form.

**Do not upload these to CRAN**

Do not upload the following as the package submission file:

- `itsne_0.1.0.zip` (Windows binary; not the CRAN source submission)
- `cran-comments.md` (supporting text, not the package tarball)
- `CRAN_SUBMISSION_CHECKLIST.md` (local checklist only)
- `CRAN_UPLOAD_FILES_CHECKLIST.md` (local checklist only)

**Supporting submission text**

Keep this file ready when filling the CRAN submission form:

- [cran-comments.md](/d:/08-MyProjects/00-MasterStudent/2024-吳玟樺/LaTeX/Article/Interval-tSNE_Wu&Wu_20260402/I-tSNE/cran-comments.md)

Use its contents in the CRAN comments box or as the submission notes text.

**Recommended final pre-upload check**

Before uploading, rerun:

- `R CMD check --as-cran itsne_0.1.0.tar.gz`

from an ASCII-only working path.

**Exact items to have ready at submission time**

1. The source tarball `itsne_0.1.0.tar.gz`
2. The text from `I-tSNE/cran-comments.md`
3. Maintainer email: `hanmingwu1103@gmail.com`
4. Repository URL: `https://github.com/hanmingwu1103/I-tSNE`
5. Release tag: `v0.1.0`

**Current release state**

- GitHub repo updated: `origin/master`
- Release commit pushed: `793e7cb`
- Release tag pushed: `v0.1.0`
- CRAN-style check from an ASCII-only path: passed with only two benign NOTEs
