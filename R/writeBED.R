writeBED <- function() {
  bed_file <- tibble(
    chrom = rep("U00096.3", length(WT_average$Position)),
    chromStart = WT_average$Position,
    chromEnd = WT_average$Position,
    name = WT_average$Position,
    score = round(WT_average$beta*100),
    strand = WT_average$Strand,
    thickStart = WT_average$Position,
    thickEnd = WT_average$Position,
    itemRGB = NA
  )
  # Deciles
  bed_file[bed_file$score <= 100, "itemRGB"] <- "250,0,0"
  bed_file[bed_file$score <= 95, "itemRGB"] <- "222,0,28"
  bed_file[bed_file$score <= 90, "itemRGB"] <- "167,0,85"
  bed_file[bed_file$score <= 80, "itemRGB"] <- "125,0,128"
  bed_file[bed_file$score <= 60, "itemRGB"] <- "0,0,255"


  cat("track name=methylation description=\"6mA methylation\" itemRgb=\"On\"\n",
      file = "/Users/carlstone/Library/CloudStorage/Box-Box/Behringer_Lab_Box_Drive/Projects/LongTermExpEvo/MethylationProject/Rscripts/WT_BED_UCSC.txt")

  write.table(bed_file,
              file = "/Users/carlstone/Library/CloudStorage/Box-Box/Behringer_Lab_Box_Drive/Projects/LongTermExpEvo/MethylationProject/Rscripts/WT_BED_UCSC.txt",
              sep = "\t",
              row.names = FALSE,
              col.names = FALSE,
              append = TRUE,
              quote = FALSE)

  bed_file$chrom <- rep("gi|556503834|ref|NC_000913.3|", length(bed_file$chrom))

  cat("track name=methylation description=\"6mA methylation\" itemRgb=\"On\"\n",
      file = "/Users/carlstone/Library/CloudStorage/Box-Box/Behringer_Lab_Box_Drive/Projects/LongTermExpEvo/MethylationProject/Rscripts/WT_BED_IGV.bed")

  write.table(bed_file,
              file = "/Users/carlstone/Library/CloudStorage/Box-Box/Behringer_Lab_Box_Drive/Projects/LongTermExpEvo/MethylationProject/Rscripts/WT_BED_IGV.bed",
              sep = "\t",
              row.names = FALSE,
              col.names = FALSE,
              append = TRUE,
              quote = FALSE)
}
