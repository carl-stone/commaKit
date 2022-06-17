# Functions for things

### Annotate methylation sites
# Use genome_sites file, which has all genome annotations together for K-12 and
#   columns 'Type' (category of feature like promoter), 'Site' (name of feature),
#   'Left' and 'Right' (start and end of feature)
# Input a file like AvWT which has a column 'start' (position of methylated site)
# Other metadata columns do not matter
annotateMethylSites <- function(methyl_df, meta_df, location) {
  for (position in methyl_df[[location]]) {
    if (nrow(meta_df[meta_df$Left <= position &
                     meta_df$Right >= position,]) == 0){
      methyl_df[methyl_df[[location]] == position, 'No_Feature'] <- '1'
      next
    }
    sites_at_position <- meta_df[meta_df$Left <= position &
                                   meta_df$Right >= position,]
    for (i in 1:nrow(sites_at_position)) {
      methyl_df[methyl_df[[location]] == position, toString(sites_at_position[i,'Type'])] <- sites_at_position[i,'Site']
    }
  }
  return(methyl_df)
}

### Annotate methylation sites relative to TSS
# Similar to the above, but only using Transcription-Units
# Use genome_sites file, which has all genome annotations together for K-12 and
#   columns 'Type' (category of feature like promoter), 'Site' (name of feature),
#   'Left' and 'Right' (start and end of feature)
# Input a file like AvWT which has a column 'start' (position of methylated site)
# Other metadata columns do not matter
annotateTSS <- function(methyl_df, meta_df, location, size) {
  meta_df <- meta_df %>% 
    filter(Type == 'Transcription-Units')
  for (position in methyl_df[[location]]) {
    if (nrow(meta_df[(meta_df$Strand == '+' &
                      meta_df$Left-size <= position &
                      meta_df$Left+size >= position) |
                     (meta_df$Strand == '-' &
                      meta_df$Right-size <= position &
                      meta_df$Right+size >= position),]) == 0) {
      methyl_df[methyl_df[[location]] == position, 'NoTSS'] <- 'X'
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
  return(methyl_df)
}

# Define function to calculate variance by coverage
# Input df, 
# Outputs dataframe with 2 columns, coverage and variance
# Calculate variance of change in evolved methylation sites by coverage
# Generate df with coverage (min to max) and variance at that coverage
# Determine which coverages are (and are not) present in dataset
var_by_coverage <- function(dataset) {
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

# Function to calculate sequencing depth across non-overlapping windows of arbitrary size
# Input dataframe, column with position, column with coverages, and window size
# Returns df with right position of each window and average depth within window
calculate_methyl_site_depth <- function(df, position_col, cov_col, w_size, calc_log2 = FALSE) {
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

# Calculate rolling mean of methylation values to show whole chromosome at once
# Default using MG1655 genome which is 4641652 long, but can put in other genome size
methylRollingMean <- function(df, position_col, methyl_col, w_size, genome_size=4641652) {
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
    if(i %% 1000 == 0) {
      print(i)
    }
    out_df[[i,'position']] <- df[[i,'position']]
    out_df[[i,'mean_methyl']] <- mean(dplyr::filter(df, position>=df[[i,'position']]&position<=(df[[i,'position']]+w_size))[['methyl']])
  }
  print(Sys.time()-tstart)
  return(out_df)
}
