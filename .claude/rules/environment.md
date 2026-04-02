# R Environment Setup

## When R Is (and Isn't) Needed

You can work without invoking R for: editing `.R` files, writing roxygen2 docs, writing vignettes, reviewing logic.

You need R when: running tests, checking the package, rebuilding docs, or verifying new code runs.

## What Is Pre-installed

Ubuntu/Debian environment. R 4.3.3 at `/usr/bin/R`.

- **Bioconductor core:** `GenomicRanges`, `IRanges`, `SummarizedExperiment`, `S4Vectors`, `GenomeInfoDb`, `Rsamtools`, `Biostrings`, `BSgenome`, `rtracklayer`, `BiocGenerics` — via `apt` (`r-bioc-*`)
- **CRAN:** `zoo`, `ggplot2`, `dplyr`, `tidyr`, `devtools`, `testthat`, `knitr`, `rmarkdown`, `ggrepel`, `patchwork` — in `/usr/local/lib/R/site-library/`
- **`BiocManager`** available for additional Bioconductor packages

## Installing Missing Packages

```bash
# Prefer apt (faster, no compilation)
sudo apt install -y r-bioc-<pkgname>   # e.g., r-bioc-genomicranges
sudo apt install -y r-cran-<pkgname>   # e.g., r-cran-zoo

# Fallback
Rscript -e "install.packages('pkgname')"
Rscript -e "BiocManager::install('pkgname')"
```

`sudo` is available without a password.

## Common Commands

```bash
Rscript -e "devtools::test()"       # run all tests
Rscript -e "devtools::check()"      # full R CMD check
Rscript -e "devtools::document()"   # rebuild man/ from roxygen2
Rscript -e "BiocCheck::BiocCheck()" # Bioconductor-specific checks
```

Do not give up if R seems unavailable — R 4.3.3 is at `/usr/bin/R` and all `comma` dependencies are pre-installed.
