# Developing commaKit

commaKit development and CI target R 4.6.1 and Bioconductor 3.23. Package
dependencies are declared in `DESCRIPTION`.

Common commands:

```r
devtools::load_all()
devtools::test()
devtools::document()
devtools::check()
```

Format R code with Air:

```bash
air format .
```

GitHub Issues track concrete work. The maintained developer notes are:

- `ROADMAP.md`: current scientific and engineering priorities
- `knowledge/architecture.md`: non-obvious package contracts
- `knowledge/design-decisions.md`: rationale behind important behavior
- `knowledge/known-issues.md`: active limitations not captured adequately by tests

User documentation belongs in the README, vignettes, and roxygen comments.
Generated `man/*.Rd` files and `NAMESPACE` should be updated with
`devtools::document()` rather than edited directly.
