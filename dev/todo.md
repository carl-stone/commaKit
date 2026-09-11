# TODOs

- `R/commaData_class.R`: `.VALID_MOD_TYPES` should not be hardcoded to the three modification types currently supported; provide a way to configure the allowed types.
- `R/parse_modkit.R` and `R/commaData_class.R`: implement binding decision 001 for exact count validation before integer conversion.
- `R/parse_modkit.R`: expand `.MODKIT_CODE_MAP` to include the full current modkit modified-base code map, including 5hmC.
- Add `rlang` as a dependency and use it consistently for argument validation and user-facing errors.
- Re-evaluate the complete modBAM checking flow to ensure the checks are sensible, correctly ordered, and actionable.
- `R/accessors.R`: remove the `commaData` method for `genome()`; use `genomeSizes()` for chromosome lengths and preserve the standard Bioconductor meaning of `genome()`.
- `R/accessors.R`: replace the `commaData` method for `annotation()` with a package-specific `featureAnnotation()` accessor that clearly returns the stored genomic feature `GRanges`.
