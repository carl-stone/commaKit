#' Title
#'
#' @param methyl_df A dataframe.
#' @param meta_df A dataframe.
#' @param location A string.
#' @param size An integer.
#' @param long A boolean.
#'
#' @return A dataframe.
#' @export
#'
#' @examples
#' WT_all_TSS <- annotateTSS(WT_average_annotated, genome_sites, location = 'Position', size = 500)
annotateTSS <-
  function(methyl_df, meta_df, location, size, long = TRUE) {
    meta_df <- meta_df %>%
      dplyr::filter(Type == 'Transcription-Units')
    for (position in methyl_df[[location]]) {
      if (nrow(meta_df[(
        meta_df$Strand == '+' &
        meta_df$Left - size <= position &
        meta_df$Left + size >= position
      ) |
      (
        meta_df$Strand == '-' &
        meta_df$Right - size <= position &
        meta_df$Right + size >= position
      ), ]) == 0) {
        methyl_df[methyl_df[[location]] == position, 'NoTSS'] <- 'X'
        next
      }
      SenseTU_at_position <- meta_df[meta_df$Strand == '+' &
                                       meta_df$Left - size <= position &
                                       meta_df$Left + size >= position, ]
      AntisenseTU_at_position <- meta_df[meta_df$Strand == '-' &
                                           meta_df$Right - size <= position &
                                           meta_df$Right + size >= position, ]
      for (i in 1:nrow(SenseTU_at_position)) {
        if (nrow(SenseTU_at_position) == 0) {
          next
        }
        methyl_df[methyl_df[[location]] == position, paste0('RelPos_+', i)] <-
          position - SenseTU_at_position[i, 'Left']
      }
      for (i in 1:nrow(AntisenseTU_at_position)) {
        if (nrow(AntisenseTU_at_position) == 0) {
          next
        }
        methyl_df[methyl_df[[location]] == position, paste0('RelPos_-', i)] <-
          AntisenseTU_at_position[i, 'Right'] - position
      }
    }
    if (long) {
      methyl_df <-
        methyl_df %>% pivot_longer(
          cols = starts_with('RelPos'),
          names_to = 'TSS_strand',
          names_pattern = 'RelPos_(.)[0-9]*',
          values_to = 'RelPos'
        )
      methyl_df <-
        methyl_df %>% distinct(Position, RelPos, .keep_all = TRUE)
      methyl_df <- methyl_df %>% dplyr::filter(!is.na(RelPos))
      return(methyl_df)
    } else {
      return(methyl_df)
    }
  }
