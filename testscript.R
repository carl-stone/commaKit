library(tidyverse)
library(devtools)
load_all()

df <- read.table(file="./data/WT_6mA_Mg.txt", header=T)
meta <- read_tsv('./all_site_annotations.txt')

df_annotated <- annotateMethylSites(df, meta, location='Position')
