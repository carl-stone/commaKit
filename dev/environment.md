# Development environment

commaKit uses a reproducible local container based on the same Bioconductor
3.23 image as CI. Package dependencies come from `DESCRIPTION`; exact versions
for the development image are recorded in `renv.lock`.

Use `dev/run` rather than invoking Docker directly:

```bash
./dev/run build
./dev/run shell
./dev/run R
./dev/run test
./dev/run test tests/testthat/test-parsers.R
./dev/run document
./dev/run validate fast
./dev/run validate pr
./dev/run validate release
```

`fast` runs Air, lintr, and the full test suite. `pr` also checks generated
roxygen files, builds the source tarball, and runs `R CMD check`. `release`
adds `--as-cran`, the PDF manual, and `BiocCheck`.

The image is rebuilt automatically when `DESCRIPTION`, `renv.lock`, or the
container configuration changes. Fast and PR validation run without network
access. All profiles operate on a temporary source copy and write logs under
`.commakit/artifacts/`.

Docker containers run as the host UID/GID. If the current user cannot access
the Docker socket directly, `dev/run` uses passwordless `sudo docker` when it
is available.
