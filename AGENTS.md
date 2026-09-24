# Instructions for coding agents

## Developing commaKit

commaKit development and CI target R 4.6.1 and Bioconductor 3.23. Declare
required package dependencies in `DESCRIPTION` (`Imports`), and optional
features, tests, and vignette dependencies in `Suggests`. Maintainer-only tools
such as `devtools` do not belong in `Imports`.

Use native R for the fast edit/test loop (install dependencies as described in
`dev/environment.md`):

```bash
./dev/run test-native
./dev/run test-native parsers  # only matching test files
```

In RStudio, use `devtools::load_all()` for interactive work and
`devtools::document()` when documentation changes. The container is the
reproducible built-package check (`fast` is an optional earlier gate; `pr`
includes the same steps):

```bash
./dev/run validate fast
./dev/run validate pr
```

See `dev/environment.md` for environment and validation details. Run `pr`
before opening a pull request.

GitHub Issues track concrete work. The maintained developer notes are:

- `dev/knowledge/architecture.md`: non-obvious package contracts
- `dev/knowledge/design-decisions.md`: provisional design ideas that have not been confirmed
- `dev/knowledge/decisions/`: numbered, binding package design decisions
- `dev/todo.md`: concise bullet points for work identified during code audits. Entries should be simple, concise, and actionable. They should contain file paths and line numbers only if they are relevant.

User documentation belongs in the README, vignettes, and roxygen comments.
Generated `man/*.Rd` files and `NAMESPACE` should be updated with
`devtools::document()` rather than edited directly.

Don't duplicate documentation between README and vignettes.

## Maintain `NEWS.md` while you work

Document all user-facing package changes as a bullet in `NEWS.md` in the same commit where the change happens. The bullet should describe the consequence for the user.

change implementation
→ change tests
→ change documentation
→ add NEWS bullet
→ commit
