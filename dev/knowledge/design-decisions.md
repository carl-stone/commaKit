# Design Decisions — Why We Made These Choices

**Last updated:** 2026-05-30
**Maintained by:** commaBot

This document records significant architectural and API decisions. When you're tempted to change something, check here first to understand why it was designed this way.

---

## D-001: diffMethyl loops by mod_context, not mod_type

**Decision:** Differential methylation tests each `mod_context` (mod_type x motif combination) independently, not each `mod_type`.

**Rationale:** Different modification contexts are biologically distinct. 6mA at GATC (Dam methyltransferase) and 6mA at ACCACC (Cellulomonas-specific) are not the same signal. Pooling them would create spurious associations.

**Example:**
- `6mA:GATC` — 393 sites, tested separately
- `5mC:CCWGG` — 195 sites, tested separately
- A site that's `6mA:GATC` is NOT tested in the `5mC:CCWGG` model

**Consequence:** Multiple testing correction is genome-wide across all mod_contexts, not per-mod_type.

**Do not change without:** Understanding the statistical implications and consulting Carl.

---

## D-002: Effect sizes on beta scale (0-1), not M-value scale

**Decision:** `dm_delta_beta` is always reported on the beta scale (0-1), even when the backend uses M-values internally (limma method).

**Rationale:** Biologists think in methylation percentages. A delta_beta of 0.2 means "20 percentage points difference" — intuitive. M-values are log-odds and harder to interpret.

**Consequence:** Limma backend transforms back to beta scale after fitting on M-values.

**Do not change:** User-facing results must always be on beta scale.

---

## D-003: Multiple testing correction is genome-wide

**Decision:** All p-values from `diffMethyl()` are corrected together across all mod_contexts in a single call to `p.adjust()`.

**Rationale:** Each site is a hypothesis test. The family of tests is all sites tested in that `diffMethyl()` call. Correcting per-mod_context would artificially lower the bar for rare modifications.

**Consequence:** If you have 393 6mA sites and 195 5mC sites, all 588 p-values are corrected together.

**Do not change without:** Statistical justification and consulting Carl.

---

## D-004: annotateSites uses list-columns (CharacterList/IntegerList/NumericList)

**Decision:** `annotateSites()` stores all overlapping/nearby features per site as list-columns, not a single match per site.

**Rationale:** Bacterial genomes are densely annotated. A single site can overlap a gene, a promoter, a TF binding site, and an operon. Forcing a single match loses information.

**Consequence:**
- `feature_types`, `feature_names`, `rel_position`, `frac_position` are all list-columns
- Intergenic sites have length-0 entries (test with `lengths(col) == 0`)
- Never revert to `!duplicated()` or `distanceToNearest()` single-match

**Do not change:** This is core to how annotation works. The test suite enforces list-column types.

---

## D-005: enrichment supports gene_role = target/regulator/both

**Decision:** `enrichMethylation()` distinguishes between target genes (where DM sites overlap gene bodies) and regulator genes (whose products bind near DM sites).

**Rationale:** Methylation can affect gene expression (target) or protein binding (regulator). These are different biological questions with different background universes.

**Consequence:**
- `gene_role = "target"` — background is all genes in annotation
- `gene_role = "regulator"` — background is only regulators of that type (e.g., sigma factors)
- `gene_role = "both"` — runs both analyses separately

**Do not change:** The semantics are intentional and documented.

---

## D-006: KEGG uses offline path with 2 bulk API calls

**Decision:** KEGG enrichment should use pre-built term-to-gene mappings via `buildKEGGTermGene()` and `buildKEGGGeneIDMap()`, not fire one HTTP request per pathway.

**Rationale:** KEGG REST API has rate limits. A typical organism has 100+ pathways. One request per pathway exceeds limits and crashes the analysis.

**Implementation:**
1. `buildKEGGTermGene(organism)` — fires 2 API calls total, caches to RDS
2. `buildKEGGGeneIDMap(organism)` — fires 1 API call, caches to RDS
3. `enrichMethylation(kegg_term2gene = ...)` — uses offline mapping, no network

**Do not change:** The offline path is the recommended approach. Online path exists but is fragile.

---

## D-007: commaData extends RangedSummarizedExperiment

**Decision:** The core data container is an S4 class extending `RangedSummarizedExperiment` (migrated from `SummarizedExperiment` in Schema v2, PR #100).

**Rationale:** Idiomatic for Bioconductor. `RangedSummarizedExperiment` stores genomic positions as `GRanges` in `rowRanges()`, enabling native `findOverlaps()`-based operations. Users familiar with DESeq2, edgeR, or other Bioc packages will recognize the pattern.

**Consequence:**
- Genomic positions stored in `rowRanges()` as 1-bp `GRanges`
- `methylation()` and `coverage()` are assay accessors
- `siteInfo()` provides backward-compatible flat DataFrame access
- Genome info from `Seqinfo` (per `seqinfo(object)`)
- Annotation and motif sites stored in `metadata()`

**Do not change:** Core to the package architecture.

---

## D-008: Modification-type agnostic by design

**Decision:** Every function must work for 6mA, 5mC, 4mC, and any future modification type. No hardcoding for specific types.

**Rationale:** The field is evolving. New modifications (e.g., phosphorothioate) may become relevant. The package should not assume a fixed set.

**Consequence:**
- No `if (mod_type == "6mA")` special cases in core logic
- `modTypes(object)` returns what's present, not a fixed list
- Statistical backends must handle any mod_type string

**Do not change:** This is a core design principle.

---

## D-009: Genome size from Seqinfo, never hardcoded

**Decision:** Any function that needs chromosome length must get it from `seqinfo(object)` (via `seqlengths()`), never hardcode values.

**Rationale:** Different organisms have different genomes. Hardcoding E. coli values would break for Helicobacter or Mycobacterium.

**Consequence:**
- `slidingWindow()` uses `seqlengths()` for boundary handling
- `plot_genome_track()` uses `seqlengths()` for axis limits
- Tests use synthetic 100 kb chromosome, not real values

**Do not change:** This is a core design principle.

---

## D-010: Parameter naming convention

**Decision:**
- Analysis functions: `verbNoun()` camelCase (e.g., `annotateSites`, `diffMethyl`)
- Plot functions: `plot_noun()` snake_case (e.g., `plot_volcano`, `plot_pca`)
- Arguments: snake_case (e.g., `mod_type`, `min_coverage`)
- Internal functions: `.` prefix (e.g., `.parseBetaValues`)

**Rationale:** Follows established R conventions. Analysis functions are actions, plots are nouns, arguments are descriptive.

**Do not change:** Consistency matters for user experience.

---

## D-011: mod_context derived on demand, not stored

**Decision:** `mod_context` is computed from `mod_type` and `motif` columns via `modContexts()` accessor, not stored as a persistent column in `mcols(rowRanges())`.

**Rationale:** Storing a derived column creates a denormalization risk — if `mod_type` or `motif` changes, the stored `mod_context` could become stale. Deriving on demand guarantees consistency.

**Consequence:** Use `modContexts(object)` to get unique mod_context strings, not `mcols(rowRanges(object))$mod_context`.

**Do not change without:** Understanding the Schema v2 migration rationale (issue #95).

---

## D-012: Alignment by genomic position, not string keys

**Decision:** All alignment between sites (e.g., merging per-sample data in the constructor) uses `findOverlaps()` on `GRanges` with `mod_type`/`motif` matching, not row names or string keys.

**Rationale:** String-key alignment is fragile (separator collisions, formatting inconsistencies). `GRanges`-based alignment is mathematically correct and leverages Bioconductor infrastructure.

**Consequence:** Assay matrices have no row names. `site_key` is a computed display column only. Its display format is stable and always includes chromosome: `chrom:position:strand:mod_type:motif`.

**Do not change without:** Understanding the no-rownames alignment design (PR #117, #105).

---

## D-013: Assay layers are additive, named, and defaulted by role

**Decision:** Assay layers are stored as named `SummarizedExperiment` assays, with a small registry in `metadata(object)$assay_provenance` and default role mapping in `metadata(object)$assay_defaults`.

**Rationale:** This follows the spirit of DESeq2 and Seurat-style layered assays: raw observations stay present, derived layers are added under explicit names, and a role can point at the current/default layer without deleting older versions. This lets multiple transformations or analysis runs coexist for audit and comparison.

**Consequence:**
- Core raw assays (`methylation`, `coverage`, `mod_counts`, `canonical_counts`) remain stable.
- Derived assays must use unique names unless overwrite is explicit.
- `assayLayers()` is the public registry view for assay names, types, provenance, parent assays, and default roles.
- Internal layer writers should use `.addAssayLayer()` instead of mutating assay metadata ad hoc.

**Do not change without:** Preserving raw-layer immutability and multi-version derived-layer behavior.

---

## How to Add New Decisions

When making a significant design or API decision:

1. Add an entry with format `D-NNN: Title`
2. Include: Decision, Rationale, Consequence, "Do not change without"
3. Cross-reference related decisions if applicable
4. Update this document, not just code comments

This document is the source of truth for "why did we do it this way."
