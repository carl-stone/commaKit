#!/usr/bin/env Rscript
# Test timing utility.
# Runs each test file individually and reports per-file timing.
# Usage: Rscript dev/test-timing.R
# Output: test-timing-report.csv and console summary

library(testthat)

if (requireNamespace("pkgload", quietly = TRUE)) {
  pkgload::load_all(".", export_all = FALSE, helpers = FALSE, quiet = TRUE)
} else {
  library(commaKit)
}

test_path <- "tests/testthat"
test_files <- list.files(
  test_path,
  pattern = "^test-.*\\.R$",
  full.names = TRUE
)
test_files <- sort(test_files)

cat("Running", length(test_files), "test files with timing...\n\n")

results <- vector("list", length(test_files))

run_test_file <- function(path) {
  test_file_formals <- names(formals(testthat::test_file))
  supports_dots <- "..." %in% test_file_formals
  args <- list(
    path = path,
    reporter = testthat::SummaryReporter$new()
  )

  add_arg <- function(name, value) {
    if (name %in% test_file_formals || supports_dots) {
      args[[name]] <<- value
    }
  }

  add_arg("package", "commaKit")
  add_arg("load_package", "none")
  add_arg("load_helpers", TRUE)
  add_arg("stop_on_failure", TRUE)
  add_arg("stop_on_warning", FALSE)

  do.call(testthat::test_file, args)
}

for (i in seq_along(test_files)) {
  f <- test_files[i]
  elapsed <- system.time({
    run_test_file(f)
  })["elapsed"]

  results[[i]] <- data.frame(
    file = basename(f),
    duration_sec = round(as.numeric(elapsed), 2),
    stringsAsFactors = FALSE
  )

  cat(sprintf("  %-50s %7.2fs\n", basename(f), as.numeric(elapsed)))
}

df <- do.call(rbind, results)
df <- df[order(-df$duration_sec), ]

cat("\n", paste(rep("=", 65), collapse = ""), "\n", sep = "")
cat("Test Timing Summary\n")
cat(paste(rep("=", 65), collapse = ""), "\n\n", sep = "")
cat(sprintf("  Total files:  %d\n", nrow(df)))
cat(sprintf("  Total time:   %.2fs\n", sum(df$duration_sec)))
cat(sprintf("  Average:      %.2fs per file\n", mean(df$duration_sec)))
cat(sprintf(
  "  Fastest:      %.2fs (%s)\n",
  min(df$duration_sec),
  df$file[which.min(df$duration_sec)]
))
cat(sprintf(
  "  Slowest:      %.2fs (%s)\n",
  max(df$duration_sec),
  df$file[which.max(df$duration_sec)]
))

cat("\n  Top 10 slowest:\n")
cat(sprintf("    %-50s %8s\n", "File", "Time(s)"))
cat("    ", paste(rep("-", 60), collapse = ""), "\n", sep = "")
for (i in seq_len(min(10, nrow(df)))) {
  cat(sprintf("    %-50s %8.2fs\n", df$file[i], df$duration_sec[i]))
}

write.csv(df, "test-timing-report.csv", row.names = FALSE)
cat(sprintf("\n  Report saved to: test-timing-report.csv\n"))
