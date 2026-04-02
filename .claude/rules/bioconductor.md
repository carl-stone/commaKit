---
paths:
  - "DESCRIPTION"
  - "NAMESPACE"
  - "NEWS.md"
---

# Bioconductor

## Requirements Checklist

| Requirement | Status |
|---|---|
| Individual package imports (`dplyr`, not `tidyverse`) | ✅ Done |
| S4 classes with proper `validity()` methods | ✅ Done |
| `show()` methods for all S4 classes | ✅ Done |
| Package-level `?comma` documentation page | ✅ Done |
| At least two vignettes | ✅ Done |
| `NEWS.md` with version history | ✅ Done |
| `biocViews` declared | ✅ Done (Sequencing, Epigenetics, Coverage, DifferentialMethylation, GenomeAnnotation, DataImport, Visualization) |
| `R CMD check --as-cran` zero errors/warnings | ⏳ Verify |
| `BiocCheck::BiocCheck()` zero errors | ⏳ Run and fix |
| Bundled data < 5 MB total | ⏳ Verify |
| Zenodo DOI | ⏳ Register before submission |
| Version bumped to 1.0.0 | ⏳ Pending |

## Out of Scope for v1.0

Do not implement without explicit discussion:

- Multi-species comparative methylomics
- Integration with transcriptomics (RNA-seq correlation)
- Motif discovery
- Phage/plasmid methylation analysis
- Shiny interactive browser
- Python or command-line interface
- Genome browser track export beyond BED (bigWig, etc.)

## Submission Checklist (when ready)

1. Update README — still reflects v0.3.0; needs full workflow + visualization showcase
2. Run `BiocCheck::BiocCheck()` and address all errors/warnings
3. Verify `R CMD check --as-cran` passes with zero errors and zero warnings
4. Register Zenodo DOI before submission
5. Verify bundled data < 5 MB total (`data/` + `inst/extdata/`)
6. Bump version to 1.0.0 in DESCRIPTION and add NEWS.md entry
7. Submit at https://contributions.bioconductor.org/
