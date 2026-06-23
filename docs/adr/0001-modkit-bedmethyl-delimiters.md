# Support safe modkit bedMethyl delimiter variants

commaKit's modkit parser should support both all-tab bedMethyl and modkit mixed-delimiter bedMethyl when that support remains small and well-tested; otherwise it should require all-tab input with a clear error. This keeps commaKit compatible with common modkit output shapes without allowing whitespace-collapsing parsing to silently shift missing fields into the wrong columns.

## Considered Options

- Require all-tab bedMethyl only.
- Support both all-tab and mixed-delimiter bedMethyl.
- Continue using whitespace-collapsing parsing.

## Consequences

Parser code must preserve enough field structure to detect blank required fields instead of repairing them by column shifting. Tests should cover both accepted delimiter shapes, or the parser should clearly reject the unsupported shape.
