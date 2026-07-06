# Keep commaKit adjusted p-values backend-independent

commaKit's `dm_padj` column means the backend-independent adjusted p-value for
the full testing family in a `diffMethyl()` call. Backend-specific statistics,
such as methylKit q-values, may be preserved as backend-specific per-site result
columns and mirrored with the active result layer, but they must not redefine
`dm_padj`.

## Considered Options

- Use methylKit q-values directly as `dm_padj` for methylKit-backed results.
- Ignore backend-specific adjusted statistics entirely.
- Keep `dm_padj` as commaKit's canonical adjusted p-value and preserve
  backend-specific statistics under backend-specific names, such as
  `dm_methylkit_qvalue`.

## Consequences

`dm_padj` stays comparable across backends and across multiple
`mod_context`s. Users who need methylKit-specific q-values can still inspect
them through `dm_methylkit_qvalue`, but the column name makes clear that they
are backend-specific evidence rather than commaKit's canonical adjusted
p-values. Backend-specific columns are present only in result layers whose
backend produced them; non-methylKit result layers should not carry an all-`NA`
`dm_methylkit_qvalue` column. Because commaKit has not had a real release with
the current methylKit result-layer contract, there is no compatibility promise
to backfill methylKit q-values into older saved objects that did not preserve
them.
