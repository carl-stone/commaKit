#' Annotate methylation sites
#'
#' Takes a dataframe with a numeric vector for methylated sites of interest (for example, all adenines within GATC sites) and a dataframe with the genomic locations of features of interest (like Promoters, Sigma factor binding sites, genes, etc.) and determines which features (can be more than one) are present at each methylation site.
#'
#' @param methyl_df A dataframe with a column listing every methylated site of interest, for example adenines within GATC sites.
#' @param meta_df A dataframe containing columns Type (feature type e.g. Promoter), Site (name of site e.g. thrLp), Left, and Right (for start and end of each feature).
#' @param location Column within methyl_df containing the genomic position of the methylation sites
#'
#' @return Dataframe methyl_df with columns added for genomic features present in each methylation site.
#' @export
#'
#' @name annotateMethylSites
#'
#' @examples
#' df <- methyl_df
#' meta <- genome_sites
#' df_annotated <- annotateMethylSites(df, meta, location='Position')
annotateMethylSites <- function(methyl_df, meta_df, location) {
  for (position in methyl_df[[location]]) {
    if (nrow(meta_df[meta_df$Left <= position &
                     meta_df$Right >= position, ]) == 0) {
      methyl_df[methyl_df[[location]] == position, 'No_Feature'] <- '1'
      next
    }
    sites_at_position <- meta_df[meta_df$Left <= position &
                                   meta_df$Right >= position, ]
    for (i in 1:nrow(sites_at_position)) {
      methyl_df[methyl_df[[location]] == position,
                toString(sites_at_position[i, 'Type'])] <-
        sites_at_position[i, 'Site']
    }
  }
  return(methyl_df)
}
