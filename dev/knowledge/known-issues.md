# Active Limitations

This file contains limitations that are not represented adequately by current
tests or user documentation. Resolved issues belong in Git history, not here.

## Real modkit files

The parser has focused tests for required fields, tab delimiters, count-derived
beta values, zero coverage, and malformed rows. It has limited validation with
real caller output across modkit versions and organisms.

## Performance at scale

Runtime and memory behavior have not been characterized systematically on full
bacterial methylomes. Synthetic package data are useful for correctness tests
but do not establish practical scaling behavior.

## Differential methylation with sparse groups

Backends differ in how they represent zero variance, perfect separation, and
insufficient observations. User-facing behavior should remain explicit and
consistent where a shared contract is scientifically defensible.

## External annotation services

KEGG-backed workflows depend on an external service when mappings are not
provided locally. Network availability, identifiers, and upstream API changes
can still cause failures; cached mappings are preferred for reproducible work.
