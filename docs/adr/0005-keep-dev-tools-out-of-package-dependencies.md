# Keep dev tools out of package dependencies

commaKit should remain Bioconductor-package-first: `DESCRIPTION` is the
authoritative installation and check contract for end users and package
automation, not a catch-all development environment manifest. Developer-only
tools such as `lintr`, `styler`, `devtools`, `roxygen2`, `pkgdown`,
`BiocCheck`, `languageserver`, and `httpgd` should stay out of `DESCRIPTION`
unless package code, tests, examples, or vignettes actually require them.

## Consequences

Contributor tooling should be installed through the devcontainer, pre-commit
setup, dedicated CI workflow package lists, or other development-environment
mechanisms. The Bioconductor-typical path is to keep `renv.lock` focused on the
check-ready package environment, including `Suggests`, rather than adding a
custom gate-tool manifest that future snapshots can accidentally drop. The
devcontainer should bootstrap `renv`, restore the lockfile, and install
contributor/editor tools once after restore. CI should remain split:
`R-CMD-check` follows the package-first `DESCRIPTION` contract, while tooling
gates such as style checks and `BiocCheck` install their tools explicitly.
Dependencies used only by archival scripts under `inst/scripts` should not
expand either dependency contract. Interactive R startup should activate `renv`
without auto-attaching development packages.
