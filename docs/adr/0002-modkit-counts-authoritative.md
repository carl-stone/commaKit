# Treat modkit counts as authoritative

commaKit should treat modkit bedMethyl count fields as authoritative for coverage and beta calculations. `fraction_modified` is derived from `Nmod / Nvalid_cov` and should not make otherwise droppable zero-coverage rows fail import before coverage filtering.

## Consequences

For retained rows with positive coverage, beta should be computed from counts rather than trusting the percentage column. Rows with zero valid coverage may be dropped by the coverage filter even when their derived fraction is missing or undefined.
