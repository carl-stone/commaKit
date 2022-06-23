calculateMethylSiteDepth <- function(df, position_col, cov_col, w_size, calc_log2 = FALSE) {
  out_df <- tibble(position = integer(),
                   coverage = double())
  for (i in 1:ceiling(max(df[[position_col]])/w_size)) {
    out_df[i,1] <- i*w_size
    out_df[i,2] <- mean(df[between(df[[position_col]], (i-1)*w_size-1, i*w_size), cov_col][[1]])
  }
  if (calc_log2 == TRUE) {
    out_df$log2_coverage <- log2(out_df$coverage)
  }
  return(out_df)
}
