# 001: Validate counts before integer conversion

commaKit follows the DESeq2 pattern for observed count data.

Before conversion to integer storage, every required count value from modkit
must be present, numeric, finite, non-negative, within R's integer range, and
exactly equal to its rounded value. Invalid values cause an actionable error;
they are never silently rounded or truncated. Valid whole-number values stored
as doubles are converted to integer storage with an informative message.

Class validity applies the same exact whole-number rule to non-missing assay
values without an epsilon tolerance. `NA` is allowed in assembled assays only
where commaKit uses it to represent a site absent from a sample. Rounding is
permitted only in a source-specific import or derivation path whose count
semantics explicitly justify it.
