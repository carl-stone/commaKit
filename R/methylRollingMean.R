function(df, position_col, methyl_col, w_size, genome_size=4641652, verbose = FALSE) {
  tstart <- Sys.time()
  out_df <- tibble(position = integer(),
                   mean_methyl = double())
  # Take the first w_size of sites from the beginning of the chromosome and add them to the end so it wraps around
  df <- bind_cols(position = df[[position_col]], methyl = df[[methyl_col]])
  nsites <- nrow(df)
  beginning_sites <- df %>% filter(position <= w_size) %>% mutate(position = position + genome_size)
  df <- bind_rows(df, beginning_sites)
  out_df <- tibble(position = as.numeric(rep(NA, nsites)),
                   mean_methyl = as.double(rep(NA, nsites)))
  for (i in 1:nsites) {
    if(i %% 1000 == 0 & verbose) {
      print(i)
    }
    out_df[[i,'position']] <- df[[i,'position']]
    out_df[[i,'mean_methyl']] <- mean(dplyr::filter(df, position>=df[[i,'position']]&position<=(df[[i,'position']]+w_size))[['methyl']])
  }
  if (verbose) print(Sys.time()-tstart)
  return(out_df)
}
