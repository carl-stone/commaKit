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

**commaKit adjusted p-value**:
The backend-independent multiple-testing-adjusted p-value reported by commaKit for a `diffMethyl()` result. It belongs to the full testing family defined by the commaKit analysis call, not to a backend-specific correction statistic.
_Avoid_: q-value, backend q-value, methylKit qvalue

**backend differential statistic**:
A statistic produced by a differential methylation backend that commaKit may preserve as backend-specific evidence. Backend differential statistics use backend-specific result columns, such as `dm_methylkit_qvalue`, and do not define commaKit's shared result columns unless explicitly promoted by the commaKit contract.
_Avoid_: commaKit result column, canonical adjusted p-value

**durable knowledge**:
Repository documentation that future maintainers use as current project memory, including `dev/knowledge/`, `CONTEXT.md`, and ADRs. Durable knowledge should change in the same branch as code or behavior changes that make it stale.
_Avoid_: stale project notes, separate cleanup later
