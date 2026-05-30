.MOD_TYPE_PALETTE <- c(
    "6mA" = "#e41a1c",
    "5mC" = "#377eb8",
    "4mC" = "#4daf4a"
)

.modTypePalette <- function(mod_type, warn_unknown = TRUE) {
    mod_type <- unique(as.character(mod_type))
    mod_type <- mod_type[!is.na(mod_type)]

    unknown <- setdiff(mod_type, names(.MOD_TYPE_PALETTE))
    if (length(unknown) > 0L) {
        if (isTRUE(warn_unknown)) {
            warning(
                "No palette color defined for mod_type value(s): ",
                paste(unknown, collapse = ", "),
                ". Using grey50 fallback.",
                call. = FALSE
            )
        }
        fallback <- stats::setNames(rep("grey50", length(unknown)), unknown)
        return(c(.MOD_TYPE_PALETTE, fallback))
    }

    .MOD_TYPE_PALETTE
}
