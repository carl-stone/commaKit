# Distinguish Dorado read-level skips from call-level drops

commaKit documentation should reserve `site` for mapped reference/genomic positions and distinguish read positions from modified-base calls during Dorado MM/ML parsing. Direct Dorado parsing may skip entire reads that cannot be parsed, but individual modified-base calls whose read positions cannot be mapped to reference sites are dropped without implying the whole read was discarded.

## Consequences

User-facing troubleshooting text should say that malformed CIGAR strings or unusable MM/ML tags can cause read-level skips, while inserted or soft-clipped modification evidence can cause call-level drops. Avoid saying that unmappable calls cause the whole read to be skipped.
