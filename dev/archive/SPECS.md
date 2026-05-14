# comma v1.0 — Implementation Specs

Agent-implementable specifications for each v1.0 blocker. Pick up a task, follow the spec, run the validation, done.

---

## Spec B1: R CMD check — zero errors/warnings

### Context
The package must pass `R CMD check --as-cran` with zero errors and zero warnings before Bioconductor submission. CI currently runs on R 4.5 (pinned due to S4Vectors C API break in R 4.6.0).

### Steps
1. Run `R CMD check --as-cran .` from the package root directory
2. Read the output. Categorize each issue as ERROR, WARNING, or NOTE
3. Fix all ERRORs and WARNINGs
4. NOTES are acceptable but minimize them (Bioconductor reviewers notice)
5. Re-run until clean

### Common R CMD check issues to watch for
- Missing `\value` in Rd files (stub docs)
- Non-ASCII characters in R code or docs (×, →, etc.)
- Unused imports in NAMESPACE
- Package version format (must be x.y.z)
- Examples that fail (use `\donttest{}` for slow/optional-dep examples)
- `library()` or `require()` calls in examples (use `::` instead)

### Validation
```bash
R CMD check --as-cran .
# Must show: 0 errors, 0 warnings
```

### Files to modify
- `R/*.R` — fix code issues
- `man/*.Rd` — fix doc issues (regenerate with `devtools::document()`)
- `DESCRIPTION` — fix dependency issues

---

## Spec B2: BiocCheck — zero errors

### Context
`BiocCheck::BiocCheck()` runs additional Bioconductor-specific checks beyond R CMD check.

### Steps
1. Install BiocCheck: `BiocManager::install("BiocCheck")`
2. Run: `BiocCheck::BiocCheck(".")`
3. Fix all ERRORs
4. WARNINGs and NOTEs are acceptable but should be minimized
5. Re-run until zero errors

### Common BiocCheck issues
- `biocViews` not properly specified
- Missing `LICENSE` or `LICENSE.md` (MIT needs file)
- Vignette builder not in `VignetteBuilder` field
- Package imports itself
- T/F instead of TRUE/FALSE
- `1:x` instead of `seq_len(x)` or `seq_along(x)`
- `system()` calls
- `install.packages()` in examples/vignettes

### Validation
```r
BiocCheck::BiocCheck(".")
# Must show: 0 errors
```

---

## Spec B3: Test file naming convention

### Context
`.claude/rules/testing.md` specifies `test-functionName.R` naming. Several test files use snake_case instead. This is cosmetic but important for Bioconductor review consistency.

### Renaming map

| Current file | New file(s) | Notes |
|---|---|---|
| `test-coverageAnalysis.R` | `test-coverageDepth.R` + `test-varianceByDepth.R` | Split into two files, one per function |
| `test-enrichment.R` | `test-enrichMethylation.R` | Rename only |
| `test-find_motif_sites.R` | `test-findMotifSites.R` | Rename only |
| `test-load_annotation.R` | `test-loadAnnotation.R` | Rename only |
| `test-m_values.R` | `test-mValues.R` | Rename only |
| `test-plot_distribution.R` | `test-plot_methylation_distribution.R` | Rename only |
| `test-parse_dorado.R` | Keep as-is (internal function, not exported) | Internal helpers don't need to match export name |
| `test-parse_megalodon.R` | Keep as-is (internal function, not exported) | Same |

### Steps
1. For each rename: `git mv tests/testthat/old.R tests/testthat/new.R`
2. For the split: create `test-coverageDepth.R` with coverageDepth tests from `test-coverageAnalysis.R`, create `test-varianceByDepth.R` with varianceByDepth tests
3. Update the test file table in `.claude/rules/testing.md`
4. Run `devtools::test()` to verify all 938 tests still pass

### Validation
```bash
Rscript -e 'library(testthat); library(comma); test_dir("tests/testthat")'
# Must show: FAIL 0, PASS >= 938
```

---

## Spec B4: Version reconciliation

### Context
- `DESCRIPTION` says `Version: 0.8.0.9000`
- `AGENTS.md` says `Version: 0.9.0.9000 dev`
- `NEWS.md` has entries for 0.8.0 and 0.9.x features
- The README roadmap stops at 0.5.0

### Decision needed from Carl
The NEWS.md suggests v0.8.0 and v0.9.x features are already implemented (mod_context, enrichment, KEGG offline). But DESCRIPTION was never bumped. Options:
1. **Bump DESCRIPTION to 0.9.2.9000** — reflects actual feature state
2. **Keep at 0.8.0.9000** — current dev version, bump only when ready for CRAN/Bioc

### Steps (after Carl decides)
1. Update `DESCRIPTION` Version field
2. Update `AGENTS.md` version line
3. Update `README.Rmd` roadmap table to include all versions through current
4. Run `devtools::document()` to regenerate man pages

### Validation
```bash
grep "^Version:" DESCRIPTION  # Should match AGENTS.md
```

---

## Spec B5: README update

### Context
README.Rmd is outdated: says "300 sites" for example data (should be 588), lists "seven plot_* functions" (should be eight), and roadmap stops at 0.5.0.

### Changes needed

1. **Line ~102:** Change "300 sites" to "588 sites" (or use `nrow(comma_example_data)` dynamically)
2. **Features section:** Change "Seven `plot_*()` functions" to "Eight `plot_*()` functions" and add `plot_tss_profile()` to the list
3. **Roadmap table:** Update to include all versions through current:
   - 0.6.0: `mValues()`, `motifs()` accessor, M-value PCA
   - 0.7.x: `plot_tss_profile()`, `diffMethyl(method="limma"|"quasi_f")`
   - 0.8.0: `mod_context` rowData column, `modContexts()`, `diffMethyl()` loops by mod_context
   - 0.9.x: `enrichMethylation()`, KEGG offline path
4. **Add enrichment workflow section** after Step 4 (Differential Methylation):
   ```r
   # Step 5 — Enrichment Analysis
   cd_dm <- diffMethyl(cd_dm, formula = ~ condition, mod_type = "6mA")
   cd_dm <- annotateSites(cd_dm, keep = "overlap")
   enr <- enrichMethylation(cd_dm, ont = "BP", gene_role = "target")
   ```
5. **Update citation** if preprint is updated

### Validation
```bash
Rscript -e 'rmarkdown::render("README.Rmd")'
# Must knit without errors
```

---

## Spec B6: Vignette — enrichment workflow

### Context
The enrichment pipeline (`enrichMethylation()`, `buildKEGGTermGene()`, `buildKEGGGeneIDMap()`) is the most complex part of comma and has no vignette coverage. Users need a worked example.

### Options
1. **Extend getting-started vignette** — add a "Step 5: Enrichment Analysis" section
2. **Create a new vignette** — `vignettes/enrichment-analysis.Rmd`

Recommendation: extend the existing vignette. It already covers steps 1–4. Adding step 5 keeps the narrative flow.

### Content outline

```
## Step 5 — Enrichment Analysis

### GO enrichment (ORA)
- Run `enrichMethylation()` with `ont = "BP"`
- Explain gene_role semantics: "target" (genes where DM sites overlap) vs "regulator" (genes whose products bind near DM sites)
- Show results extraction and visualization

### KEGG enrichment (offline path)
- Why offline: KEGG API rate limits
- Build term2gene mapping: `buildKEGGTermGene("eco")`
- Build gene ID map: `buildKEGGGeneIDMap("eco", OrgDb = org.EcK12.eg.db)`
- Pass to `enrichMethylation(kegg_term2gene = ..., kegg_term2name = ...)`
- Show results

### GSEA mode
- When to use GSEA vs ORA
- Run `enrichMethylation(method = "gsea")`
- Show gene score metric and interpretation
```

### Steps
1. Edit `vignettes/getting-started.Rmd`
2. Add enrichment section after differential methylation
3. Use `\donttest{}` for chunks requiring optional dependencies (clusterProfiler, KEGGREST)
4. Knit to verify: `rmarkdown::render("vignettes/getting-started.Rmd")`
5. Update `.claude/rules/documentation.md` vignette list

### Validation
```bash
Rscript -e 'rmarkdown::render("vignettes/getting-started.Rmd")'
# Must knit without errors
```

---

## Spec P1: Verify bundled data < 5 MB

### Steps
1. Check total size of `data/` and `inst/extdata/`:
   ```bash
   du -sh data/ inst/extdata/
   ```
2. If > 5 MB, identify large files and consider:
   - Moving to ExperimentHub (Bioconductor's cloud data hosting)
   - Compressing further
   - Removing unnecessary files

### Validation
```bash
du -sh data/ inst/extdata/
# Must be < 5 MB combined
```

---

## Spec P2: Register Zenodo DOI

### Steps
1. Go to https://zenodo.org/
2. Sign in with GitHub
3. Find the carl-stone/CoMMA repository
4. Enable the Zenodo webhook
5. Create a release on GitHub to trigger DOI minting
6. Add the DOI badge to README.Rmd

### Validation
README.Rmd should have a DOI badge.

---

## Spec P3: Review biocViews

### Current
```
biocViews: Sequencing, Epigenetics, Coverage, DifferentialMethylation
```

### Recommended additions
```
biocViews: Sequencing, Epigenetics, Coverage, DifferentialMethylation,
    GenomeAnnotation, DataImport, Visualization, Microbiome
```

### Steps
1. Edit `DESCRIPTION`
2. Verify all biocViews terms are valid: https://bioconductor.org/packages/release/BiocViews.html
3. Run `BiocCheck::BiocCheck()` to verify

### Validation
```r
BiocCheck::BiocCheck(".")
# No errors about biocViews
```
