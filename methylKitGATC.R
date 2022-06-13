# Adapting methylKit to our data
# Made copy of AllSamples.6mA.Megalodon.Split_By_Strand.Header.noblanks.txt for the analysis
# Renamed it to mgAllSamples.txt for simplicity

library(tidyverse)
library(methylKit)
library(cowplot)
library(ggrepel)
library(readr)
library(vegan)
source('Projects/LongTermExpEvo/MethylationProject/Rscripts/functions.R')

# Laptop
# setwd("/Users/carlstone/Documents/Vanderbilt/Behringer_Lab/Methylation")
# Mac
setwd("/Users/carlstone/Library/CloudStorage/Box-Box/Behringer_Lab_Box_Drive")



# I don't actually use this here but useful for later
mgAllSamples <- as_tibble(read.table(file='Projects/LongTermExpEvo/MethylationProject/ComparativeAnalysis/RawDataFiles/SortedFiles/MethylPercent_Coverage_All_NoBlank.txt', sep='\t', header=T))
mgAllSamples <- mgAllSamples %>% 
  dplyr::mutate(across(c(cov_WT1:methyl_129.3), as.numeric)) %>% 
  drop_na()
mgAllSamplesLong <- mgAllSamples %>% 
  pivot_longer(
    cols = cov_WT1:methyl_129.3,
    names_to = c(".value", "sample"),
    names_sep = "_",
    values_drop_na = TRUE
  )

mgFileList <- list(
  mgWT1 <- 'Projects/LongTermExpEvo/MethylationProject/ComparativeAnalysis/RawDataFiles/SortedFiles/WT1.6mA.Megalodon.Filtered.sort.txt',
  mgWT2 <- 'Projects/LongTermExpEvo/MethylationProject/ComparativeAnalysis/RawDataFiles/SortedFiles/WT2.6mA.Megalodon.Filtered.sort.txt',
  mgWT3 <- 'Projects/LongTermExpEvo/MethylationProject/ComparativeAnalysis/RawDataFiles/SortedFiles/WT3.6mA.Megalodon.Filtered.sort.txt',
  mg113.1 <- 'Projects/LongTermExpEvo/MethylationProject/ComparativeAnalysis/RawDataFiles/SortedFiles/113_1.6mA.Megalodon.Filtered.sort.txt',
  mg113.2 <- 'Projects/LongTermExpEvo/MethylationProject/ComparativeAnalysis/RawDataFiles/SortedFiles/113_2.6mA.Megalodon.Filtered.sort.txt',
  mg125.6 <- 'Projects/LongTermExpEvo/MethylationProject/ComparativeAnalysis/RawDataFiles/SortedFiles/125_6.6mA.Megalodon.Filtered.sort.txt',
  mg125.7 <- 'Projects/LongTermExpEvo/MethylationProject/ComparativeAnalysis/RawDataFiles/SortedFiles/125_7.6mA.Megalodon.Filtered.sort.txt',
  mg126.2 <- 'Projects/LongTermExpEvo/MethylationProject/ComparativeAnalysis/RawDataFiles/SortedFiles/126_2.6mA.Megalodon.Filtered.sort.txt',
  mg126.4 <- 'Projects/LongTermExpEvo/MethylationProject/ComparativeAnalysis/RawDataFiles/SortedFiles/126_4.6mA.Megalodon.Filtered.sort.txt',
  mg129.1 <- 'Projects/LongTermExpEvo/MethylationProject/ComparativeAnalysis/RawDataFiles/SortedFiles/129_1.6mA.Megalodon.Filtered.sort.txt',
  mg129.3 <- 'Projects/LongTermExpEvo/MethylationProject/ComparativeAnalysis/RawDataFiles/SortedFiles/129_3.6mA.Megalodon.Filtered.sort.txt'
)

mgDataList <- list(
  mgWT1Data = read.table(mgWT1, sep='\t', header=T, na.strings='NA'),
  mgWT2Data = read.table(mgWT2, sep='\t', header=T, na.strings='NA'),
  mgWT3Data = read.table(mgWT3, sep='\t', header=T, na.strings='NA'),
  mg113.1 = read.table(mg113.1, sep='\t', header=T, na.strings='NA'),
  mg113.2 = read.table(mg113.2, sep='\t', header=T, na.strings='NA'),
  mg125.6 = read.table(mg125.6, sep='\t', header=T, na.strings='NA'),
  mg125.7 = read.table(mg125.7, sep='\t', header=T, na.strings='NA'),
  mg126.2 = read.table(mg126.2, sep='\t', header=T, na.strings='NA'),
  mg126.4 = read.table(mg126.4, sep='\t', header=T, na.strings='NA'),
  mg129.1 = read.table(mg129.1, sep='\t', header=T, na.strings='NA'),
  mg129.3 = read.table(mg129.3, sep='\t', header=T, na.strings='NA')
)
# Read in table of mutated sites
# Don't do these next two commands if you want to remove mutated GATC sites after diff meth analysis
siteMutations <- read_tsv(file = 'Projects/LongTermExpEvo/MethylationProject/GATC_Mutation_Data/SummaryData/GATCsiteSummary.txt')
# This is an incredibly janky way to make a list of positions of mutated GATC sites
# Fix this later because it's not really exact but it's close enough to work
mutatedGATCs <- rep(siteMutations$position, each=7) + seq(-3,3)

# Drop NA rows, delete rows with GATC mutations, save as new tsv files, and make list of new file names
mgNewFileList <- vector(mode="list", length = length(mgFileList))
for (i in 1:length(mgDataList)) {
  mgDataList[[i]] <- as.data.frame(mgDataList[[i]]) %>% drop_na()
  #mgDataList[[i]] <- as.data.frame(mgDataList[[i]]) %>% filter(!(Position %in% mutatedGATCs))
  write.table(mgDataList[i], file=paste0('/Users/carlstone/Documents/Methylation', names(mgDataList)[i], '.txt'), quote=F, sep='\t', row.names=F)
  mgNewFileList[i] <- paste0('/Users/carlstone/Documents/Methylation', names(mgDataList)[i], '.txt')
}


# Selects all WT together
# mgAllSamplesLong %>% dplyr::filter(grepl('WT', sample))
# Use filter when you have to extract a sample from the long dataset
# mgWT1 <- mgAllSamplesLong %>% dplyr::filter(sample == 'WT1')
# mgWT2 <- mgAllSamplesLong %>% dplyr::filter(sample == 'WT2')
# mgWT3 <- mgAllSamplesLong %>% dplyr::filter(sample == 'WT3')
# mg113.1 <- mgAllSamplesLong %>% dplyr::filter(sample == '113.1')
# mg113.2 <- mgAllSamplesLong %>% dplyr::filter(sample == '113.2')
# mg125.6 <- mgAllSamplesLong %>% dplyr::filter(sample == '125.6')
# mg125.7 <- mgAllSamplesLong %>% dplyr::filter(sample == '125.7')
# mg126.2 <- mgAllSamplesLong %>% dplyr::filter(sample == '126.2')
# mg126.4 <- mgAllSamplesLong %>% dplyr::filter(sample == '126.4')
# mg129.1 <- mgAllSamplesLong %>% dplyr::filter(sample == '129.1')
# mg129.3 <- mgAllSamplesLong %>% dplyr::filter(sample == '129.3')

# Convert coverage and percent to number of methylated and unmethylated reads
# mgAllSamplesLong$nMethyl = mgAllSamplesLong$methyl * mgAllSamplesLong$cov
# mgAllSamplesLong$nUnmethyl = (1 - mgAllSamplesLong$methyl) * mgAllSamplesLong$cov


meta_cols <- list(fraction = TRUE,
                  chr.col = 1,
                  start.col = 2,
                  end.col = 2,
                  coverage.col = 4,
                  strand.col = 3,
                  freqC.col = 5)

mgMethylObject <- methRead(mgNewFileList, 
                           sample.id = list('WT1', 'WT2', 'WT3', 'mg113.1', 'mg113.2',
                                                        'mg125.6', 'mg125.7', 'mg126.2', 'mg126.4',
                                                        'mg129.1', 'mg129.3'),
                           assembly = 'NC_000913.3',
                           treatment = c(0,0,0,2,2,1,1,2,2,1,1),
                           context = 'GATC',
                           mincov=10,
                           pipeline = meta_cols
                           )



# Normalize coverage by the median coverage
mgMethylObject <- normalizeCoverage(mgMethylObject, method = 'median')


# Plot histogram of methylation per base for WT2
getMethylationStats(mgMethylObject[[2]],plot=FALSE,both.strands=FALSE)
# Now plot coverage
getCoverageStats(mgMethylObject[[1]],plot=FALSE,both.strands=FALSE)

# Combine all samples into one dataframe
# numC's is the number of methylated reads
# numT's is the number of unmethylated reads
allMethyl <- unite(mgMethylObject)
percent_allMethyl <- cbind(allMethyl$start, percMethylation(allMethyl))
colnames(percent_allMethyl)[1] <- 'position'
# write_tsv(as.data.frame(percent_allMethyl), file='/Users/carlstone/Box/Behringer_Lab_Box_Drive/Projects/LongTermExpEvo/MethylationProject/ComparativeAnalysis/DifferentialMethylation/percentMethylationMatrixNormalizedMedian_noMutantGATCs.txt')

##################### Start here for future analysis ###########################
percent_allMethyl <- read_tsv(file='Projects/LongTermExpEvo/MethylationProject/ComparativeAnalysis/DifferentialMethylation/percentMethylationMatrixNormalizedMedian_noMutantGATCs.txt')
percent_allMethyl$position <- percent_allMethyl$position + 1

getCorrelation(allMethyl,plot=FALSE)
clusterSamples(allMethyl, dist="correlation", method="ward", plot=TRUE)

# Screeplot
PCASamples(allMethyl, screeplot=T)
# Output prcomp object
pca_prcomp <- PCASamples(allMethyl, scale=F, center=T, obj.return=T)
# Plot first 2 PCs
PCASamples(allMethyl, center=T, scale=F, transpose=T, sd.filter=T, filterByQuantile = T)

# PCA with percent methyl matrix
pca_percent <- prcomp(t(percMethylation(allMethyl)))
library(ggfortify)
autoplot(pca_percent)

# Try NMDS because why the fuck not
sites <- as.data.frame(percent_allMethyl)$position
percent_allMethyl_t <- as_tibble(percent_allMethyl[,2:ncol(percent_allMethyl)])
percent_allMethyl_t <- t(percent_allMethyl_t)
methyl_distmat <- vegdist(percent_allMethyl_t, method='euclidean')
methyl_distmat <- as.matrix(methyl_distmat, labels=T)
NMDS <- metaMDS(methyl_distmat,
                distance = 'euclidean')
goodness(NMDS)
stressplot(NMDS)
screeplot(NMDS)
plot(NMDS, 'sites')
orditorp(NMDS, 'sites')
groups_NMDS <- rep(NA, 11)
groups_NMDS[1:3] <- 'Ancestor'
groups_NMDS[c(4,5,8,9)] <- 'MMR-'
groups_NMDS[c(6,7,10,11)] <- 'WT'
#ggplot it
NMDS_scores <- as.data.frame(scores(NMDS))
NMDS_scores$Sample <- rownames(NMDS_scores)
NMDS_scores$group <- groups_NMDS
plot_NMDS <- ggplot(NMDS_scores, aes(x=NMDS1, y=NMDS2, color=group)) +
  geom_point(size=2) +
  theme_bw() +
  theme(legend.text = element_text(size = 14),
        legend.position = c(0.75, 0.75),
        legend.background = element_rect(fill='white', color='black')) +
  labs(color = 'Group') +
  scale_color_manual(values=c('#619CFF', 'red', 'black')) +
  coord_fixed(ratio=1.5625)
ggsave(filename = 'NMDS.png',
       path = 'Projects/LongTermExpEvo/MethylationProject/Manuscript/Figures',
       plot = plot_NMDS)

# Calculate differential methylation
diffMethyl = calculateDiffMeth(allMethyl)
# Run differential methylation pairwise to compare each evolved to ancestor
methylAncvWT <- reorganize(allMethyl,
                           sample.ids = c('WT1', 'WT2', 'WT3', 'mg125.6', 'mg125.7', 'mg129.1', 'mg129.3'),
                           treatment = c(0,0,0,1,1,1,1))
methylAncvMMR <- reorganize(allMethyl,
                            sample.ids = c('WT1', 'WT2', 'WT3', 'mg113.1', 'mg113.2', 'mg126.2', 'mg126.4'),
                            treatment = c(0,0,0,2,2,2,2))
diffMethylAvWT <- calculateDiffMeth(methylAncvWT)
diffMethylAvMMR <- calculateDiffMeth(methylAncvMMR)
# Convert to tibbles for easy graphing
diffMethylAvWT_tibble <- as_tibble(diffMethylAvWT)
diffMethylAvMMR_tibble <- as_tibble(diffMethylAvMMR)

# Read Gwyn's table that has metadata
allSites <- read_tsv(file='/Users/carlstone/Box/Behringer_Lab_Box_Drive/Projects/LongTermExpEvo/MethylationProject/ComparativeAnalysis/RawDataFiles/WT1_2.Megalodon.Filtered.txt')
# Add whether each position is within a gene and which gene
diffMethylAvWT_tibble_metadata <- left_join(diffMethylAvWT_tibble, allSites[,c('location_types', 'position', 'genic_location')],
                                            by = c("start" = "position"))
diffMethylAvMMR_tibble_metadata <- left_join(diffMethylAvMMR_tibble, allSites[,c('location_types', 'position', 'genic_location')],
                                            by = c("start" = "position"))
# Save tsvs
#write_delim(diffMethylAvWT_tibble_metadata, file='/Users/carlstone/Box/Behringer_Lab_Box_Drive/Projects/LongTermExpEvo/MethylationProject/ComparativeAnalysis/DifferentialMethylation/diffMethAncVsWT.txt', delim='\t')
#write_delim(diffMethylAvMMR_tibble_metadata, file='/Users/carlstone/Box/Behringer_Lab_Box_Drive/Projects/LongTermExpEvo/MethylationProject/ComparativeAnalysis/DifferentialMethylation/diffMethAncVsMMR.txt', delim='\t')

# Read tsvs to change graphs later
#diffMethylAvWT_tibble_metadata <- read_tsv(file='Projects/LongTermExpEvo/MethylationProject/ComparativeAnalysis/DifferentialMethylation/diffMethAncVsWT.txt')
#diffMethylAvMMR_tibble_metadata <- read_tsv(file='Projects/LongTermExpEvo/MethylationProject/ComparativeAnalysis/DifferentialMethylation/diffMethAncVsMMR.txt')

#diffMethylAvMMR_tibble_metadata <- diffMethylAvMMR_tibble_metadata %>% filter(!(start %in% mutatedGATCs))

# tsvs with mutated GATCs removed
diffMethylAvWT_tibble_metadata <- read_tsv(file='Projects/LongTermExpEvo/MethylationProject/ComparativeAnalysis/DifferentialMethylation/diffMethAncVsWT_noMutantGATCs.txt')
diffMethylAvMMR_tibble_metadata <- read_tsv(file='Projects/LongTermExpEvo/MethylationProject/ComparativeAnalysis/DifferentialMethylation/diffMethAncVsMMR_noMutantGATCs.txt')

# get hyper methylated bases
myDiff25p.hyper=getMethylDiff(diffMethyl,difference=25,qvalue=0.01,type="hyper")
# get hypo methylated bases
# Doesn't work rn
myDiff25p.hypo=getMethylDiff(diffMethyl,difference=1,qvalue=0.5,type="hypo")
# get all differentially methylated bases
myDiff25p=getMethylDiff(diffMethyl,difference=25,qvalue=0.01)
dmAvWT10p <- getMethylDiff(diffMethylAvWT, difference=10, qvalue=0.01)
dmAvMMR10p <- getMethylDiff(diffMethylAvMMR, difference=10, qvalue=0.01)

# Order by methylation difference
diffMethyl[order(diffMethyl$meth.diff),]

### Volcano plots
# Determine whether GATC site was mutated or not
siteMutations <- read_tsv(file = 'Projects/LongTermExpEvo/MethylationProject/GATC_Mutation_Data/SummaryData/GATCsiteSummary.txt')
diffMethylAvWT_tibble_metadata$mutated <- 'NO'
diffMethylAvMMR_tibble_metadata$mutated <- 'NO'
mutatedGATCs <- rep(siteMutations$position, each=7) + seq(-3,3)
diffMethylAvMMR_tibble_metadata[diffMethylAvMMR_tibble_metadata$start %in% mutatedGATCs,]$mutated <- 'YES'
# Define % difference cutoff and q value cutoff for both of the following volcano plots
percentCutoff <- 10
qvalueCutoff <- 0.05
# Classify genes as diffmeth up or down for later coloring
diffMethylAvWT_tibble_metadata$diffexpressed <- "NO"
diffMethylAvWT_tibble_metadata$diffexpressed[diffMethylAvWT_tibble_metadata$meth.diff > percentCutoff & 
                                      diffMethylAvWT_tibble_metadata$qvalue < qvalueCutoff] <- "UP"
diffMethylAvWT_tibble_metadata$diffexpressed[diffMethylAvWT_tibble_metadata$meth.diff < -percentCutoff & 
                                      diffMethylAvWT_tibble_metadata$qvalue < qvalueCutoff] <- "DOWN"
diffMethylAvMMR_tibble_metadata$diffexpressed <- "NO"
diffMethylAvMMR_tibble_metadata$diffexpressed[diffMethylAvMMR_tibble_metadata$meth.diff > percentCutoff & 
                                      diffMethylAvMMR_tibble_metadata$qvalue < qvalueCutoff] <- "UP"
diffMethylAvMMR_tibble_metadata$diffexpressed[diffMethylAvMMR_tibble_metadata$meth.diff < -percentCutoff & 
                                      diffMethylAvMMR_tibble_metadata$qvalue < qvalueCutoff] <- "DOWN"
# Label differentially methylated sites
diffMethylAvWT_tibble_metadata$dmlabel <- NA
diffMethylAvWT_tibble_metadata$dmlabel[diffMethylAvWT_tibble_metadata$diffexpressed != "NO"] <- diffMethylAvWT_tibble_metadata$genic_location[diffMethylAvWT_tibble_metadata$diffexpressed != "NO"]
diffMethylAvMMR_tibble_metadata$dmlabel <- NA
diffMethylAvMMR_tibble_metadata$dmlabel[diffMethylAvMMR_tibble_metadata$diffexpressed != "NO"] <- diffMethylAvMMR_tibble_metadata$genic_location[diffMethylAvMMR_tibble_metadata$diffexpressed != "NO"]


# Define colors for each group
mycolors <- c('blue', 'red', 'black')
names(mycolors) <- c("DOWN", "UP", "NO")
# Define shapes for mutated/not mutated
myshapes <- c(15,16)
names(myshapes) <- c("YES", "NO")
# Plot the plots
WT_volcano <- ggplot(data=diffMethylAvWT_tibble_metadata, aes(x=meth.diff, y=-log10(qvalue), col=diffexpressed)) +
  geom_point(aes(shape = mutated)) +
  scale_shape_manual(values = myshapes) +
  theme_minimal() +
  theme(legend.position = 'none') +
  scale_x_continuous(limits = c(-65, 65), breaks = seq(-75, 75, 25)) +
#  scale_x_continuous(limits = c(-52.5, 52.5), breaks = seq(-75, 75, 25)) +
  scale_y_continuous(limits = c(NA, 250)) +
#  scale_y_continuous(limits = c(NA, 200)) +
  geom_vline(xintercept = c(-percentCutoff, percentCutoff), col='red') +
  geom_hline(yintercept = -log10(qvalueCutoff), col='red') +
  scale_colour_manual(values=mycolors) +
  labs(x = 'Percent Methylation Difference')
WT_volcano_inset <- ggplot(data=diffMethylAvWT_tibble_metadata, aes(x=meth.diff, y=-log10(qvalue), col=diffexpressed, label=dmlabel)) +
#  geom_point(aes(shape = mutated)) +
  scale_shape_manual(values = myshapes) +
  theme_minimal() +
  theme(legend.position = 'none', plot.background = element_rect(fill = 'white', colour = 'black')) +
  scale_y_continuous(limits = c(NA, NA)) +
  geom_vline(xintercept = c(-percentCutoff, percentCutoff), col='red') +
  geom_hline(yintercept = -log10(qvalueCutoff), col='red') +
  scale_colour_manual(values=mycolors) +
  labs(x = 'Percent Methylation Difference') +
  geom_text_repel() +
  coord_fixed(ratio = 2)
# Add inset volcano plot to WT using cowplot ggdraw
WT_volcano_combined <- ggdraw(WT_volcano) +
  draw_plot(WT_volcano_inset,
            x = 0.25,
            y = 0.4,
            width = 0.7, 
            height = 0.7)
MMR_volcano <- ggplot(data=diffMethylAvMMR_tibble_metadata, aes(x=meth.diff, y=-log10(qvalue), col=diffexpressed, label=dmlabel)) +
  geom_point(aes(shape = mutated)) +
  scale_shape_manual(values = myshapes) +
  theme_minimal() +
  theme(legend.position = "none") +
  theme(axis.title.y = element_blank()) +
#  scale_x_continuous(limits = c(-85, 20), breaks = seq(-75, 75, 25)) + 
  scale_y_continuous(limits = c(NA, 250)) +
#  scale_y_continuous(limits = c(NA, 200)) +
  geom_vline(xintercept = c(-percentCutoff, percentCutoff), col='red') +
  geom_hline(yintercept = -log10(qvalueCutoff), col='red') +
  scale_colour_manual(values=mycolors) +
  labs(x = 'Percent Methylation Difference', col='Differential Methylation') +
  geom_text_repel()
diffMethylFigure <- plot_grid(WT_volcano_combined, MMR_volcano, labels = 'AUTO')

MMR_volcano_small <- ggplot(data=diffMethylAvMMR_tibble_metadata, aes(x=meth.diff, y=-log10(qvalue), col=diffexpressed)) +
  geom_point() +
  scale_shape_manual(values = myshapes) +
  theme_minimal() +
  theme(legend.position = 'none') +
  scale_x_continuous(limits = c(-25,25)) +
  scale_y_continuous(limits = c(NA, 25)) +
  geom_vline(xintercept = c(-percentCutoff, percentCutoff), col='red') +
  geom_hline(yintercept = -log10(qvalueCutoff), col='red') +
  scale_colour_manual(values=mycolors) +
  labs(x = 'Percent Methylation Difference') +
  coord_fixed(ratio = 2) +
  annotate('text', 
           x=c(-19,19), y=c(20,20), 
           label=c('141','28'),
           size=c(7,7),
           color=c('blue','red'))
WT_volcano_small <- ggplot(data=diffMethylAvWT_tibble_metadata, aes(x=meth.diff, y=-log10(qvalue), col=diffexpressed)) +
  geom_point() +
  theme_minimal() +
  theme(legend.position = 'none') +
  scale_x_continuous(limits = c(-25,25)) +
  scale_y_continuous(limits = c(NA, 25)) +
  geom_vline(xintercept = c(-percentCutoff, percentCutoff), col='red') +
  geom_hline(yintercept = -log10(qvalueCutoff), col='red') +
  scale_colour_manual(values=mycolors) +
  labs(x = '') +
  coord_fixed(ratio=2) +
  annotate('text', 
           x=c(-19,19), y=c(20,20), 
           label=c('16','22'),
           size=c(7,7),
           color=c('blue','red'))
smallPlots_noMutants <- plot_grid(WT_volcano_small, MMR_volcano_small, labels = c("MMR+", "MMR-"))

# Save it
ggsave(file = paste0('diffMethylFigure_',Sys.Date(),'.pdf'),
       path = '/Users/carlstone/Box/Behringer_Lab_Box_Drive/Projects/LongTermExpEvo/MethylationProject/ComparativeAnalysis/DifferentialMethylation',
       plot = diffMethylFigure,
       device = 'pdf')

# Save small plot with no mutants
ggsave(file = paste0('diffMethyl_noMutants_', Sys.Date(),'.pdf'),
       path = '/Users/carlstone/Box/Behringer_Lab_Box_Drive/Projects/LongTermExpEvo/MethylationProject/ComparativeAnalysis/DifferentialMethylation',
       plot = smallPlots_noMutants,
       device = 'pdf')

# Combined NMDS and volcano plots
vertical_volcano <- plot_grid(NULL,WT_volcano_small,NULL, MMR_volcano_small,
                              rel_heights = c(0.15,1,-0.4,1),
                              label_x = 0,
                              label_y = 0.855,
                              labels = c('','A','','B'),
                              ncol = 1)
cow_NMDS <- plot_grid(NULL, plot_NMDS,
                      rel_heights = c(0.2,1.85),
                      ncol = 1,
                      label_y = 0.92,
                      labels = c('','C'))
all_plots <- plot_grid(vertical_volcano, NULL, cow_NMDS,
                       rel_widths = c(2,0.1,4),
                       labels = c('','', ''),
                       nrow = 1)

ggsave(file = 'combined_dM_data.png',
       path = 'Projects/LongTermExpEvo/MethylationProject/ComparativeAnalysis/DifferentialMethylation',
       plot = all_plots,
       width = 7, height = 6)

################# Compare methylation at different features between ############
################# ancestor and evolved #########################################
# First read percent_allMethyl above
# Optional step below, average WT into one value
percent_allMethyl <- percent_allMethyl %>% 
  mutate(WT = rowMeans(select(., starts_with('WT')), na.rm=T)) %>% 
  select(!matches('WT[123]'))

genome_sites <- read_tsv(file='Projects/LongTermExpEvo/MethylationProject/Ecoli_Annotation/all_site_annotations.txt')
all_sites_annotated <- annotateMethylSites(percent_allMethyl, genome_sites, location = 'position')

# Elongate on site annotations to make it easier to summarize
# and remove useless rows (where feature_name = NA)
all_annotated_long <- all_sites_annotated %>% 
  pivot_longer(
    cols = 'Transcription-Units':'Origin-of-Replication',
    names_to = 'feature_type',
    values_to = 'feature_name'
  ) %>% 
  filter(!is.na(feature_name)) %>% 
  pivot_longer(cols = 'mg113.1':'WT',
               names_to = 'sample',
               values_to = 'methylation')

# Add group and genotype to each sample, then elongate those
MMRplus_genotypes <- list('mg125.6','mg125.7','mg129.1','mg129.3')
MMRminus_genotypes <- list('mg113.1','mg113.2','mg126.2','mg126.4')

all_annotated_long$group = NA
all_annotated_long$genotype = NA

all_annotated_long <- all_annotated_long %>% 
  mutate(group = if_else(sample == 'WT', 'ancestor', 'evolved'),
         genotype = if_else(sample %in% MMRminus_genotypes, 'MMR-', 'MMR+'),
         groupxgenotype = if_else(sample == 'WT', 'ancestor', paste(group, genotype, sep='_')))

all_annotated_long %>% 
  filter(feature_type != 'Promoters' & feature_type != 'Transcription-Units') %>%  
  filter(str_detect(feature_type, 'Sigma', negate = TRUE)) %>% 
  kruskal.test(methylation ~ sample, data = .)

pwWilcoxTest <- all_annotated_long %>% 
  filter(feature_type != 'Promoters' & feature_type != 'Transcription-Units') %>% 
  filter(str_detect(feature_type, 'Sigma', negate = TRUE)) %>% 
  {pairwise.wilcox.test(x = .$methylation,
                        g = .$sample,
                        p.adjust.method = 'BH')}

all_features_plot <- all_annotated_long %>% 
  filter(feature_type != 'Promoters' &
           feature_type != 'Transcription-Units') %>% 
#  filter(feature_type == 'DNA-Binding-Sites') %>% 
  ggplot(aes(x = feature_type, y = methylation)) +
  geom_boxplot(notch = T) +
  facet_wrap(~groupxgenotype)

sigma_avg_methylation_plot <- all_annotated_long %>% 
  filter(str_detect(feature_type, 'Sigma')) %>% 
  ggplot(aes(x = feature_type, y = methylation, color = groupxgenotype)) +
  geom_boxplot(outlier.size = 0.5,
               outlier.alpha = 0.3,
               notch = T) +
  theme_classic()


# Make a boxplot of methylation values by feature type
# First define pairwise comparisons of interest
comparison_list <- rev(list(c('Sigma24','Sigma38'),
                            c('Sigma24','Sigma70'),
                            c('Sigma28','Sigma38'),
                            c('Sigma28','Sigma70')))

#comparison_list <- rev(list(c('Sigma24','Sigma28'),
#                        c('Sigma24','Sigma32'),
#                        c('Sigma24','Sigma38'),
#                        c('Sigma24','Sigma70'),
#                        c('Sigma28','Sigma32'),
#                        c('Sigma28','Sigma38'),
#                        c('Sigma28','Sigma70'),
#                        c('Sigma32','Sigma38'),
#                        c('Sigma32','Sigma70'),
#                        c('Sigma38','Sigma70')))

feature_avg_methylation_plot <- WT_annotated_long %>% 
  filter(feature_type != 'Promoters' &
           feature_type != 'Transcription-Units') %>% 
  filter(str_detect(feature_type, 'Sigma', negate = TRUE)) %>% 
  ggplot(aes(x = feature_type, y = Mg_MethylPercent_Mean)) +
  geom_boxplot(outlier.size = 0.5,
               outlier.alpha = 0.3,
               notch = T) +
  theme_classic()

sigma_avg_methylation_plot <- WT_annotated_long %>% 
  filter(str_detect(feature_type, 'Sigma')) %>% 
  ggplot(aes(x = feature_type, y = Mg_MethylPercent_Mean)) +
  geom_boxplot(outlier.size = 0.5,
               outlier.alpha = 0.3,
               notch = T) +
  theme_classic() +
  stat_compare_means(method = 'wilcox.test',
                     comparisons = comparison_list,
                     label = 'p.signif')


ggsave(plot = feature_avg_methylation_plot,
       filename = 'WT_all_features_methylation_plot.png',
       path = 'Projects/LongTermExpEvo/MethylationProject/Manuscript/Figures',
       width = 22,
       height = 7)

ggsave(plot = sigma_avg_methylation_plot,
       filename = 'WT_sigma_methylation_plot.png',
       path = 'Projects/LongTermExpEvo/MethylationProject/Manuscript/Figures',
       width = 7,
       height = 7)









