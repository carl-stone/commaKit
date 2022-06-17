function(methyl_df, meta_df, location) {
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
