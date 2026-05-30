# R Environment Setup

## When R Is (and Isn't) Needed

You can work without invoking R for: editing `.R` files, writing roxygen2 docs, writing vignettes, reviewing logic.

You need R when: running tests, checking the package, rebuilding docs, or verifying new code runs.

## What Is Pre-installed

The package requires R >= 4.3.0, while CI is pinned to R 4.5 on Ubuntu. Local macOS development should also use R 4.5 for now because some compiled Bioconductor packages can be incompatible with newer framework R builds.

Use the repository's `renv.lock` and the active R 4.5 library when possible. On Carl's macOS machines, `dev/run-r-4.5.sh` can be used as a convenience wrapper if the default framework R is newer than CI.

## Installing Missing Packages

```bash
# Prefer apt (faster, no compilation)
sudo apt install -y r-bioc-<pkgname>   # e.g., r-bioc-genomicranges
sudo apt install -y r-cran-<pkgname>   # e.g., r-cran-zoo

# Fallback
Rscript -e "install.packages('pkgname')"
Rscript -e "BiocManager::install('pkgname')"
```

Avoid installing packages outside the project library unless you are deliberately fixing the development environment.

## Common Commands

```bash
Rscript -e "devtools::test()"       # run all tests
Rscript -e "devtools::check()"      # full R CMD check
Rscript -e "devtools::document()"   # rebuild man/ from roxygen2
Rscript -e "BiocCheck::BiocCheck()" # Bioconductor-specific checks
```

For CI-equivalent local checks, use R 4.5 and run `R CMD check --no-manual --as-cran` after `devtools::document()` and `devtools::test()`.
