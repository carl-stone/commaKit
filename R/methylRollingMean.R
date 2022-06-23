methylRollingMean <- function(df, position_col, methyl_col, w_size=1000, genome_size=4641652, verbose=FALSE) {
  tstart <- Sys.time()
  # Take the first w_size of sites from the beginning of the chromosome and add them to the end so it wraps around
  df <- dplyr::bind_cols(position = df[[position_col]], methyl = df[[methyl_col]])
  nsites <- nrow(df)
  beginning_sites <- df %>%
    dplyr::filter(position <= w_size) %>%
    dplyr::mutate(position = position + genome_size)
  df <- dplyr::bind_rows(df, beginning_sites)
  df <- data.matrix(df)
  # TODO sort df by position
  out_df <- data.frame(position = as.numeric(rep(NA, nsites)),
                   mean_methyl = as.double(rep(NA, nsites)))
  out_df <- data.matrix(out_df)

  for (i in 1:nsites) {
    if(verbose && i %% 1000 == 0) {
      print(i)
    }
    out_df[[i,'position']] <- df[[i,'position']]
    out_df[[i,'mean_methyl']] <- mean(df[df[,1] >= df[[i,1]] & df[,1] < (df[[i,1]]+w_size), 2])
  }

  out_df <- as_tibble(out_df)
  if (verbose) print(Sys.time()-tstart)
  return(out_df)
}
