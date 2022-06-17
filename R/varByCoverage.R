function(dataset) {
  f_all_coverage <- rep(NA, max(dataset$Coverage_Sample))
  for (cov in 1:max(dataset$Coverage_Sample)) {
    if (cov %in% dataset$Coverage_Sample) {
      f_all_coverage[cov] <- TRUE
    } else {
      f_all_coverage[cov] <- FALSE
    }
  }
  #Create vector containing all coverage levels present
  f_coverage_levels <- rep(NA, sum(f_all_coverage))
  cov_i <- 1
  for (cov in 1:length(f_all_coverage)) {
    if (f_all_coverage[cov] == TRUE) {
      f_coverage_levels[cov_i] <- cov
      cov_i <- cov_i + 1
    }
  }
  #Calculate variance of all coverage levels
  f_coverage_dist <- as.data.frame(cbind(f_coverage_levels, rep(NA, length(f_coverage_levels))))
  names(f_coverage_dist)[1] <- "coverage"
  names(f_coverage_dist)[2] <- "variance"
  for (cov in 1:nrow(f_coverage_dist)) {
    temp <- dataset[dataset$Coverage_Sample == f_coverage_dist[cov,1],]
    f_coverage_dist[cov,2] <- var(temp["Percent_Methyl_Sample"]-temp["Ancestor_Mean"])
  }
  return(f_coverage_dist)
}
