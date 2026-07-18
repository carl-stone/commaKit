# Load only explicit assay and result layer registries

commaKit supports explicit named layer registries as its assay and differential
methylation compatibility surface. Assay metadata is read from the named
`metadata(object)$assay_provenance` record list and
`metadata(object)$assay_defaults` role map. Differential methylation metadata is
read from the named `diffMethyl_results` and `diffMethyl_result_layers`
registries plus the explicit `diffMethyl_default_result` name.

The active differential result remains mirrored into bare `dm_*` row metadata
columns for existing downstream workflows. The mirror is not a standalone
result registry.

## Supported legacy shapes

| Shape | Load expectation |
|---|---|
| Named assay provenance records | Load records by assay name; missing record fields remain `NA`. |
| Named assay default role map | Mark only the mapped assays as defaults. |
| `default_for` in an assay provenance record | Fill a missing role-map entry for that record. |
| Named differential result tables | Load only as named result data when a corresponding named layer record exists. |
| Named differential result-layer records | List and resolve the named layer; an absent result table produces the existing aligned-table error. |
| Explicit differential default name | Use only this name for default resolution; no name is inferred. |
| Bare active `dm_*` row metadata | Preserve and use as the mirror of an active named result; do not reconstruct a layer from it. |

## Unsupported inference

`assayLayers()` does not infer provenance, roles, types, sources, or defaults
from assay names. Result helpers do not reconstruct layers from
`diffMethyl_result_cols`, `diffMethyl_params`, a sole result table, or bare
`dm_*` columns without named registries.

## Consequences

Incomplete legacy objects remain loadable as `commaData` objects and can expose
their explicitly named layers, but missing metadata or defaults remains visible
instead of being silently fabricated. New `diffMethyl()` calls write the named
registries, the explicit active name, and the bare mirror only.
