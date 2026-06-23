# commaKit

commaKit is an R/Bioconductor package for comparative microbial methylomics.
This glossary captures project-specific language used in parser, data model, and analysis decisions.

## Language

**modkit bedMethyl input**:
A modkit pileup bedMethyl file accepted by commaKit's modkit parser. commaKit should support both all-tab bedMethyl and modkit mixed-delimiter bedMethyl when the implementation remains small and well-tested; otherwise it should require all-tab input with a clear error.
_Avoid_: modkit TSV, generic BED

**modkit count fields**:
The authoritative modkit bedMethyl values for commaKit coverage and methylation beta calculations: `Nvalid_cov`, `Nmod`, `Ncanonical`, and `Nother_mod`. `fraction_modified` is derived from these fields and should not override them.
_Avoid_: modkit fraction as source of truth

**site**:
A 1-bp reference/genomic position represented by a commaKit row after caller evidence has been mapped and aggregated. Do not use `site` for an unmapped base offset within a sequencing read.
_Avoid_: read offset, raw base call

**read position**:
A base offset within an individual sequencing read before or during alignment to a reference. A read position may fail to map to any reference site, for example when it is inserted or soft-clipped.
_Avoid_: genomic site

**modified-base call**:
Caller evidence that a read position carries a specific base modification, usually with probability or count support. Modified-base calls become commaKit site evidence only after they map to reference sites.
_Avoid_: site

**durable knowledge**:
Repository documentation that future agents or maintainers use as current project memory, including `dev/knowledge/`, `CONTEXT.md`, and ADRs. Durable knowledge should change in the same branch as code or behavior changes that make it stale.
_Avoid_: stale project notes, separate cleanup later
