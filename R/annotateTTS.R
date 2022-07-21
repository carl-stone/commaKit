annotateTTS <- function(methyl_df, meta_df, location, size, long = TRUE, annotatedOnly = FALSE) {
  if (annotatedOnly == TRUE) {
    meta_df <- meta_df %>%
      filter(Type == 'Rho-Independent-Terminators' |
               Type == 'Rho-Dependent-Terminators')
    for (position in methyl_df[[location]]) {
      if (nrow(meta_df[(meta_df$Strand == '+' &
                        meta_df$Left-size <= position &
                        meta_df$Left+size >= position) |
                       (meta_df$Strand == '-' &
                        meta_df$Right-size <= position &
                        meta_df$Right+size >= position),]) == 0) {
        methyl_df[methyl_df[[location]] == position, 'NoTTS'] <- 'X'
        next
      }
      SenseTU_at_position <- meta_df[meta_df$Strand == '+' &
                                       meta_df$Left-size <= position &
                                       meta_df$Left+size >= position,]
      AntisenseTU_at_position <- meta_df[meta_df$Strand == '-' &
                                           meta_df$Right-size <= position &
                                           meta_df$Right+size >= position,]
      for (i in 1:nrow(SenseTU_at_position)) {
        if(nrow(SenseTU_at_position) == 0) {
          next
        }
        methyl_df[methyl_df[[location]] == position, paste0('RelPos_+',i)] <- position-SenseTU_at_position[i,'Left']
      }
      for (i in 1:nrow(AntisenseTU_at_position)) {
        if(nrow(AntisenseTU_at_position) == 0) {
          next
        }
        methyl_df[methyl_df[[location]] == position, paste0('RelPos_-',i)] <- AntisenseTU_at_position[i,'Right']-position
      }
    }
  } else if (annotatedOnly == FALSE) {
    meta_df <- meta_df %>%
      filter(Type == 'Transcription-Units')
    for (position in methyl_df[[location]]) {
      if (nrow(meta_df[(meta_df$Strand == '+' &
                        meta_df$Right-size <= position &
                        meta_df$Right+size >= position) |
                       (meta_df$Strand == '-' &
                        meta_df$Left-size <= position &
                        meta_df$Left+size >= position),]) == 0) {
        methyl_df[methyl_df[[location]] == position, 'NoTTS'] <- 'X'
        next
      }
      SenseTU_at_position <- meta_df[meta_df$Strand == '+' &
                                       meta_df$Right-size <= position &
                                       meta_df$Right+size >= position,]
      AntisenseTU_at_position <- meta_df[meta_df$Strand == '-' &
                                           meta_df$Left-size <= position &
                                           meta_df$Left+size >= position,]
      for (i in 1:nrow(SenseTU_at_position)) {
        if(nrow(SenseTU_at_position) == 0) {
          next
        }
        methyl_df[methyl_df[[location]] == position, paste0('RelPos_+',i)] <- position-SenseTU_at_position[i,'Right']
      }
      for (i in 1:nrow(AntisenseTU_at_position)) {
        if(nrow(AntisenseTU_at_position) == 0) {
          next
        }
        methyl_df[methyl_df[[location]] == position, paste0('RelPos_-',i)] <- AntisenseTU_at_position[i,'Left']-position
      }
    }
  }
  if (long == TRUE) {
    methyl_df <- methyl_df %>% pivot_longer(cols = starts_with('RelPos'),
                                            names_to = 'TTS_strand',
                                            names_pattern = 'RelPos_(.)[0-9]*',
                                            values_to = 'RelPos')
    methyl_df <- methyl_df %>% distinct(Position, RelPos, .keep_all = TRUE)
    methyl_df <- methyl_df %>% filter(!is.na(RelPos))
    return(methyl_df)
  } else {
    return(methyl_df)
  }
}
