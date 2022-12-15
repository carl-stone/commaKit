#' Title
#'
#' @param df A df with columns position_col and methyl_col.
#' @param position_col A string.
#' @param methyl_col A string.
#' @param w_size An integer for sliding window size.
#' @param genome_size An integer, by default the E. coli K-12 MG1655 genome size.
#' @param method A string "exact" (default) or "fast"
#'
#' @return A dataframe.
#' @export
#'
#' @examples
#' sliding_window_methylation <- methylRollingMedian(WT_average, position_col = 'Position', methyl_col = 'beta', w_size = 10000, method = "exact")
methylRollingMedian <- function(df,
                                position_col,
                                methyl_col,
                                w_size,
                                genome_size = 4641652,
                                method = "exact") {
  # Combine input into one df
  df <- bind_cols(position = df[[position_col]], methyl = df[[methyl_col]])

  if (method == "exact") {
    # Create empty df of every genomic position + extra for the end to wrap around
    all_pos <- data.frame(position = seq(1, genome_size + w_size, 1),
                          methyl = NA)
    all_pos[all_pos$position %in% df$position,]$methyl <- df$methyl
    # Take the first w_size of sites from the beginning of the chromosome and add them to the end
    all_pos[all_pos$position > genome_size, "methyl"] <- all_pos[all_pos$position <= w_size, "methyl"]
    # Calculate median sliding window
    all_pos$med_methyl <- zoo::rollapply(all_pos$methyl, w_size + 1,
                                         median, na.rm = TRUE, partial = TRUE,
                                         align = "center")
    out_df <- all_pos[1:genome_size,]
  }

  if (method == "fast") {
    nsites <- nrow(df)
    beginning_sites <- df %>%
      filter(position <= w_size) %>%
      mutate(position = position + genome_size)
    df <- bind_rows(df, beginning_sites)
    df <- data.matrix(df)
    # TODO sort df by position
    out_df <- tibble(position = as.numeric(rep(NA, nsites)),
                     mean_methyl = as.double(rep(NA, nsites)))
    out_df <- data.matrix(out_df)
    for (i in 1:nsites) {
      if (i %% 1000 == 0) {
        print(i)
      }
      out_df[[i,'position']] <- df[[i,'position']]
      out_df[[i,'mean_methyl']] <- median(df[df[,1] >= df[[i,1]] & df[,1] < (df[[i,1]] + w_size), 2])
    }
    out_df <- as_tibble(out_df)
  }
  return(out_df)
}
