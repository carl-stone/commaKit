# Development environment

The fast development loop uses native R 4.6.1 and Bioconductor 3.23. Install
dependencies from `DESCRIPTION` when setting up or changing dependencies:

```r
if (!requireNamespace("pak", quietly = TRUE)) install.packages("pak")
pak::local_install_dev_deps(upgrade = FALSE)
if (!requireNamespace("devtools", quietly = TRUE)) install.packages("devtools")
```

Use `devtools::load_all()` for interactive work, `devtools::test()` for tests,
and `devtools::document()` when documentation changes. From a terminal:

```bash
./dev/run test-native
./dev/run test-native parsers  # test files matching "parsers"
```

`test-native` uses `devtools::test()` and does not invoke Docker. It does not
install the package; this keeps the edit/test loop fast. For a built-package
check before a PR, use the local container based on the same Bioconductor 3.23
image as CI. `DESCRIPTION` remains the package dependency contract;
`renv.lock` records one known-good set of package versions for this local
validation environment. CI instead installs dependencies from `DESCRIPTION`.

Use `dev/run` rather than invoking Docker directly:

```bash
./dev/run build
./dev/run shell
./dev/run R
./dev/run test                  # container tests, including package installation
./dev/run test tests/testthat/test-parsers.R
./dev/run document
./dev/run validate fast
./dev/run validate pr
./dev/run validate release
```

`fast` runs Air, lintr, and the full test suite. `pr` also checks generated
roxygen files, builds the source tarball, and runs `R CMD check`. `release`
adds `--as-cran`, the PDF manual, and `BiocCheck`.

The image is rebuilt only when `renv.lock` or the container configuration
changes. Editing `DESCRIPTION` (including its version) does not rebuild the
image. Before running a container command, the current `DESCRIPTION`
dependencies are checked against the installed packages; if a new dependency
is missing, update `renv.lock` and rebuild the image. Container tests and
validation run without network access. Validation profiles operate on a
temporary source copy and write logs under `.commakit/artifacts/`.
An interactive `./dev/run shell` skips this check so a missing dependency can
be diagnosed; scripted shell commands do not.

Docker containers run as the host UID/GID. If the current user cannot access
the Docker socket directly, `dev/run` uses passwordless `sudo docker` when it
is available.
