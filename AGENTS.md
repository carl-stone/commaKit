# Instructions for coding agents

## Developing commaKit

commaKit development and CI target R 4.6.1 and Bioconductor 3.23. Package dependencies are declared in `DESCRIPTION`. Dependencies should only include packages needed for end users to run the package, not for development.

Common commands:

```r
devtools::load_all()    # load package into R session
devtools::test()        # run tests
devtools::document()    # update documentation
devtools::check()       # run checks
```

Format R code with Air:

```bash
air format .
```

GitHub Issues track concrete work. The maintained developer notes are:

- `knowledge/architecture.md`: non-obvious package contracts
- `knowledge/design-decisions.md`: rationale behind important behavior

User documentation belongs in the README, vignettes, and roxygen comments.
Generated `man/*.Rd` files and `NAMESPACE` should be updated with
`devtools::document()` rather than edited directly.

## Maintain `NEWS.md` while you work

Document all user-facing package changes as a bullet in `NEWS.md` in the same commit where the change happens. The bullet should describe the consequence for the user.

change implementation
→ change tests
→ change documentation
→ add NEWS bullet
→ commit
