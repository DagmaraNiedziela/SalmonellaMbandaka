library(readxl) 

# FASTQC results ####
merged_fastqc_summaries <- read_excel("fastqc results/merged_fastqc_summaries.xlsx", sheet = 1)
View(merged_fastqc_summaries)
colnames(merged_fastqc_summaries)

library(tidyr)
pivoted_fastqc_summaries <- pivot_wider(merged_fastqc_summaries, id_cols = Sample, 
                                        names_from = Statistic, values_from = Result)
View(pivoted_fastqc_summaries)

write.csv(pivoted_fastqc_summaries, "Pivoted_fastqc_summaries.csv") 

# BacMet EXP vs PRE #### 

BacMet_EXP_n <- read_excel("BacMet_results_separate.xlsx", sheet = "EXP")
BacMet_EXP_n
nrow(BacMet_EXP_n) #1882 

BacMet_PRE <- read_excel("BacMet_results_separate.xlsx", sheet = "PRE") 
nrow(BacMet_PRE) #1891 
1882+1891 #3773, previous no of rows when joined by everything 

EXP_vs_PRE <- full_join(BacMet_EXP_n, BacMet_PRE, by = c("Gene_Alias", "File"), suffix = c("_EXP", "_PRE"))
EXP_vs_PRE
View(EXP_vs_PRE)
nrow(EXP_vs_PRE) # 2614
?full_join
#by = "Gene_Alias", 
write.csv(EXP_vs_PRE, "EXP_vs_PRE.csv") 

EXP_vs_PRE_inner <- inner_join(BacMet_EXP_n, BacMet_PRE, by = c("Gene_Alias", "File"), suffix = c("_EXP", "_PRE"))
nrow(EXP_vs_PRE_inner) #1416 
write.csv(EXP_vs_PRE_inner, "EXP_vs_PRE_inner.csv") 

# Get total length into the EXP file ####

gene_legend <- read_excel("BacMet_results_separate.xlsx", sheet = "gene legend2") 
View(gene_legend) 
gene_legend <- gene_legend %>% select(Gene_Alias, `Total length`)

BacMet_EXP_n_length <- full_join(BacMet_EXP_n, gene_legend, by = "Gene_Alias") 
View(BacMet_EXP_n_length) 

BacMet_EXP_n_length$`Total length` <- as.numeric(as.character(BacMet_EXP_n_length$`Total length`))
BacMet_EXP_n_length <- BacMet_EXP_n_length %>% mutate(coverage=`Match length`/`Total length`*100)

BacMet_EXP_n_length_reduced <- BacMet_EXP_n_length %>% filter(coverage > 85) %>% filter(`Percent identity` > 95)
write_csv(BacMet_EXP_n_length_reduced, "BacMet_EXP_n_length_reduced.csv")
nrow(BacMet_EXP_n_length_reduced) #1128 
# THIS IS THE BEST OPTION 

colnames(BacMet_EXP_n_length)
smaller <- BacMet_EXP_n_length %>% filter(coverage < 90)
View(smaller) 

# ** THIS figure BacMet heatmap on reduced file ####

BacMet_EXP_n_length_reduced_cols <- BacMet_EXP_n_length_reduced %>% select(File, Gene_Alias, `Percent identity`)
BacMet_EXP_reduced_for_heatmap <- BacMet_EXP_n_length_reduced_cols %>% pivot_wider(names_from = Gene_Alias, values_from = `Percent identity`)

View(BacMet_EXP_reduced_for_heatmap) 

BacMet_EXP_reduced_for_heatmap_df <- apply(BacMet_EXP_reduced_for_heatmap,2,as.character)
write.csv(BacMet_EXP_reduced_for_heatmap_df, "BacMet_EXP_reduced_for_heatmap_df.csv") 
# Manually fix 

BacMet_EXP_reduced_for_heatmap_df_cp <- read.csv("BacMet_EXP_reduced_for_heatmap_df2.csv", na.strings = "NA", row.names = "File")
head(BacMet_EXP_reduced_for_heatmap_df_cp)

# Two colours for 0 and 1 
library(RColorBrewer) 
RColorBrewer::brewer.pal.info
RColorBrewer::display.brewer.all()
bacmet_cols <- colorRampPalette(brewer.pal(9, "RdPu") )(9) # or PiYG, last bracket was 100
bacmet_cols
bacmet_cols2 <- c("#FA9FB5","#F768A1","#DD3497","#AE017E","#7A0177","#49006A")
bacmet_cols3 <- c("#FFF7F3","#FCC5C0","#F768A1","#AE017E","#49006A")
?colorRampPalette

?brewer.pal
cols <- brewer.pal(11, "RdYlGn")
cols2 <- brewer.pal(11, "PiYG")
cols3 <- brewer.pal(11, "YlOrRd")

# Italicise column names 
newcolnames <- lapply(
  colnames(BacMet_EXP_reduced_for_heatmap_df_cp),
  function(x) bquote(italic(.(x))))
newcolnames

library(pheatmap)
?pheatmap
BacMet_reduced_heatmap <- pheatmap(BacMet_EXP_reduced_for_heatmap_df_cp, col = bacmet_cols2, na_col = "grey", 
                                      cluster_cols = FALSE, cluster_rows = TRUE, angle_col = 45, breaks = c(95,96,97,98,99,100), labels_col = as.expression(newcolnames)) 
BacMet_reduced_heatmap 

BacMet_reduced_heatmap2 <- pheatmap(BacMet_EXP_reduced_for_heatmap_df_cp, col = bacmet_cols3, na_col = "grey", 
                                   cluster_cols = FALSE, cluster_rows = TRUE, angle_col = 45, breaks = c(95,96,97,98,99,100), labels_col = as.expression(newcolnames))

BacMet_reduced_heatmap2
#cellwidth = 10, fontsize = 12

tiff("BacMet_reduced_heatmap_new2.tiff", height = 30, width = 30, units = 'cm', 
     compression = "lzw", res = 600) 
BacMet_reduced_heatmap2
dev.off() 

# And try Cairo 
library(Cairo)
Cairo(file="BacMet_reduced_heatmap_new2.png", 
      type="png",
      units="in", 
      width=10, 
      height=12, 
      pointsize=12, 
      dpi=300)
BacMet_reduced_heatmap2
dev.off()

# BacMet results #### 

library(readxl) 
BacMet_EXP <- read_excel("BacMet_results_separate.xlsx", sheet = 1) 
BacMet_EXP 
colnames(BacMet_EXP) 
nrow(BacMet_EXP) #1882
#  [1] "File"             "Query"            "Subject"          "Gene_Alias"       "Gene"             "Description"      "Location"        
# [8] "Organism"         "Compounds"        "Percent identity" "Match length"     "E-value"          "Score per length"

library(tidyverse) 
BacMet_EXP_Salmonella <- BacMet_EXP %>% filter_at((vars(starts_with("Salmonella"))), all_vars(. > 0)) 
BacMet_EXP_Salmonella <- BacMet_EXP %>% filter(grepl("Salmonella", BacMet_EXP$Organism)) 
nrow(BacMet_EXP_Salmonella) #806 
View(BacMet_EXP_Salmonella)

BacMet_EXP_identity <- BacMet_EXP %>% filter(`Percent identity` > 95) 

nrow(BacMet_EXP_identity) #1188
BacMet_EXP_identity 
View(BacMet_EXP_identity)

# For now make a heatmap with all of the variables 
library(tidyr)
BacMet_EXP_columns <- BacMet_EXP %>% select(File, Gene_Alias, `Percent identity`)
BacMet_EXP_for_heatmap <- BacMet_EXP_columns %>% pivot_wider(names_from = Gene_Alias, values_from = `Percent identity`)
BacMet_EXP_columns 
View(BacMet_EXP_for_heatmap)

View(BacMet_EXP_columns)
BacMet_EXP_columns_internals <- BacMet_EXP_columns[83:1599,] 
BacMet_EXP_columns_internals_heatm <- BacMet_EXP_columns_internals %>% pivot_wider(names_from = Gene_Alias, values_from = `Percent identity`)
View(BacMet_EXP_columns_internals_heatm)
BacMet_EXP_columns_internals_heatm_df <- apply(BacMet_EXP_columns_internals_heatm,2,as.character)
write.csv(BacMet_EXP_columns_internals_heatm_df, "BacMet_EXP_columns_internals_heatm.csv") 

BacMet_EXP_Salmonella_internals <- BacMet_EXP_Salmonella[33:704,]
BacMet_EXP_Salmonella_internals_columns <- BacMet_EXP_Salmonella_internals %>% select(File, Gene_Alias, `Percent identity`)
BacMet_EXP_Salmonella_for_heatmap <- BacMet_EXP_Salmonella_internals_columns %>% pivot_wider(names_from = Gene_Alias, values_from = `Percent identity`) 

View(BacMet_EXP_Salmonella_for_heatmap)
BacMet_EXP_Salmonella_for_heatmap_df <- as.data.frame(BacMet_EXP_Salmonella_for_heatmap)
BacMet_EXP_Salmonella_for_heatmap_df2 <- apply(BacMet_EXP_Salmonella_for_heatmap_df,2,as.character)
write.csv(BacMet_EXP_Salmonella_for_heatmap_df2, "BacMet_EXP_Salmonella_internals_for_heatmap.csv") 

BacMet_EXP_Salmonella_columns <- BacMet_EXP_Salmonella %>% select(File, Gene_Alias, `Percent identity`)
BacMet_EXP_Salmonella_for_heatmap2 <- BacMet_EXP_Salmonella_columns %>% pivot_wider(names_from = Gene_Alias, values_from = `Percent identity`) 
BacMet_EXP_Salmonella_for_heatmap_df3 <- apply(BacMet_EXP_Salmonella_for_heatmap2,2,as.character)
write.csv(BacMet_EXP_Salmonella_for_heatmap_df3, "BacMet_EXP_Salmonella_all_for_heatmap.csv") 

# ** Heatmap ##### 

# Salmonella, internals only 
# CAREFUL! This file is no longer correct, I would have to make a new one 
BacMet_EXP_Salmonella_for_heatmap_cp <- read.csv("BacMet_EXP_Salmonella_internals_for_heatmap3.csv", na.strings = "NULL", row.names = "File")
head(BacMet_EXP_Salmonella_for_heatmap_cp)

#virulence_mbandaka_cp[is.na(virulence_mbandaka_cp)] <- as.double("NA")

# Two colours for 0 and 1 
library(RColorBrewer) 
RColorBrewer::brewer.pal.info
RColorBrewer::display.brewer.all()
bacmet_cols <- colorRampPalette(brewer.pal(9, "RdPu") )(100) # or PiYG
cols <- brewer.pal(11, "RdYlGn")
cols2 <- brewer.pal(11, "PiYG")
cols3 <- brewer.pal(11, "YlOrRd")

library(pheatmap)
BacMet_Salmonella_heatmap <- pheatmap(BacMet_EXP_Salmonella_for_heatmap_cp, col = bacmet_cols, na_col = "grey", 
                         cluster_cols = FALSE, cluster_rows = FALSE) 
BacMet_Salmonella_heatmap

#cellwidth = 10, fontsize = 12

tiff("BacMet_Salmonella_heatmap.tiff", height = 30, width = 30, units = 'cm', 
     compression = "lzw", res = 600) 
BacMet_Salmonella_heatmap
dev.off() 

# And try Cairo 
library(Cairo)
Cairo(file="BacMet_Salmonella_heatmap.png", 
      type="png",
      units="in", 
      width=10, 
      height=10, 
      pointsize=12, 
      dpi=72)
BacMet_Salmonella_heatmap
dev.off()

# Salmonella, internals only - columns sorted alphabetically 
BacMet_EXP_Salmonella_for_heatmap_cp2 <- read.csv("BacMet_EXP_Salmonella_internals_for_heatmap2.csv", na.strings = "NULL", row.names = "File")
head(BacMet_EXP_Salmonella_for_heatmap_cp2)

#virulence_mbandaka_cp[is.na(virulence_mbandaka_cp)] <- as.double("NA")

# Two colours for 0 and 1 
library(RColorBrewer) 
RColorBrewer::brewer.pal.info
RColorBrewer::display.brewer.all()
bacmet_cols <- colorRampPalette(brewer.pal(9, "RdPu") )(100) # or PiYG
cols <- brewer.pal(11, "RdYlGn")
cols2 <- brewer.pal(11, "PiYG")
cols3 <- brewer.pal(11, "YlOrRd")

library(pheatmap)
BacMet_Salmonella_heatmap2 <- pheatmap(BacMet_EXP_Salmonella_for_heatmap_cp2, col = bacmet_cols, na_col = "grey", 
                                      cluster_cols = FALSE, cluster_rows = FALSE) 
BacMet_Salmonella_heatmap2

#cellwidth = 10, fontsize = 12

tiff("BacMet_Salmonella_heatmap2.tiff", height = 30, width = 30, units = 'cm', 
     compression = "lzw", res = 600) 
BacMet_Salmonella_heatmap2
dev.off() 

# And try Cairo 
library(Cairo)
Cairo(file="BacMet_Salmonella_heatmap2.png", 
      type="png",
      units="in", 
      width=10, 
      height=10, 
      pointsize=12, 
      dpi=72)
BacMet_Salmonella_heatmap2
dev.off() 

# I don't really like the alphabetical way 

# *** Salmonella, all strains #####
BacMet_EXP_Salmonella_for_heatmap_all_cp <- read.csv("BacMet_EXP_Salmonella_all_for_heatmap2.csv", na.strings = "NULL", row.names = "File")
head(BacMet_EXP_Salmonella_for_heatmap_all_cp)

#virulence_mbandaka_cp[is.na(virulence_mbandaka_cp)] <- as.double("NA")

# Two colours for 0 and 1 
library(RColorBrewer) 
RColorBrewer::brewer.pal.info
RColorBrewer::display.brewer.all()
bacmet_cols <- colorRampPalette(brewer.pal(9, "RdPu") )(100) # or PiYG
cols <- brewer.pal(11, "RdYlGn")
cols2 <- brewer.pal(11, "PiYG")
cols3 <- brewer.pal(11, "YlOrRd")

library(pheatmap)
BacMet_Salmonella_heatmap_all <- pheatmap(BacMet_EXP_Salmonella_for_heatmap_all_cp, col = bacmet_cols, na_col = "grey", 
                                      cluster_cols = FALSE, cluster_rows = FALSE, angle_col = 45) 
?pheatmap
BacMet_Salmonella_heatmap_all

#cellwidth = 10, fontsize = 12

tiff("BacMet_Salmonella_heatmap_all.tiff", height = 30, width = 30, units = 'cm', 
     compression = "lzw", res = 600) 
BacMet_Salmonella_heatmap_all
dev.off() 

# And try Cairo 
# This is a good size image with clear enough font, it takes A4 
Cairo(file="BacMet_Salmonella_heatmap_all.png", 
      type="png",
      units="in", 
      width=10, 
      height=15, 
      pointsize=12, 
      dpi=300)
BacMet_Salmonella_heatmap_all
dev.off()


# Abricate results #### 

library(dplyr)

# ** VFDB ####

virulence_results <- read_excel("vfdb_results.xlsx", sheet = 2)
head(virulence_results)
nrow(virulence_results) #10268

library(janitor)
virulence_results <- virulence_results %>% clean_names() 
colnames(virulence_results)

virulence_results_reduced <- virulence_results %>% 
  filter(percent_coverage > 80) %>% filter(percent_identity > 95)
virulence_results_reduced 
nrow(virulence_results_reduced) #6159 
virulence_results_reduced_df <- as.data.frame(virulence_results_reduced)
head(virulence_results_reduced_df)
virulence_results_reduced_df <- virulence_results_reduced_df[,-1]
rownames(virulence_results_reduced_df) <- NULL

write.csv(virulence_results_reduced_df, "vfdb_reduced.csv") 

# ** Plasmid #### 

plasmid_results <- read.csv("Analyses results no SPADES/MLST and abricate results/abricate_plasmid.csv")
head(plasmid_results)
nrow(plasmid_results) #258 
colnames(plasmid_results)

plasmid_results_reduced <- plasmid_results %>% 
  filter(`X.COVERAGE`> 80) %>% filter(`X.IDENTITY`> 95)
plasmid_results_reduced 
nrow(plasmid_results_reduced) #178
plasmid_results_reduced_df <- as.data.frame(plasmid_results_reduced)
head(plasmid_results_reduced_df)
 
write.csv(plasmid_results_reduced_df, "plasmid_reduced.csv", row.names = FALSE) 
write.table(plasmid_results_reduced_df, "plasmid_reduced.tab", sep = "\t", quote = FALSE, row.names = FALSE) 

# ** Argannot ##### 

argannot_results <- read.csv("Analyses results no SPADES/MLST and abricate results/abricate_argannot_amr.csv")
head(argannot_results)
nrow(argannot_results) #275 
colnames(argannot_results)

argannot_results_reduced <- argannot_results %>% 
  filter(`X.COVERAGE`> 80) %>% filter(`X.IDENTITY`> 95)
argannot_results_reduced 
nrow(argannot_results_reduced) #100
argannot_results_reduced_df <- as.data.frame(argannot_results_reduced)
head(argannot_results_reduced_df)

write.csv(argannot_results_reduced_df, "argannot_reduced.csv", row.names = FALSE) 
write.table(argannot_results_reduced_df, "argannot_reduced.tab", sep = "\t", quote = FALSE, row.names = FALSE) 

# ** Card #### 

card_results <- read.csv("Analyses results no SPADES/MLST and abricate results/abricate_card_amr.csv")
head(card_results)
nrow(card_results) #3345 
colnames(card_results)

card_results <- card_results %>% clean_names() 
colnames(card_results)

card_results_reduced <- card_results %>% 
  filter(x_coverage > 80) %>% filter(x_identity > 95)
card_results_reduced 
nrow(card_results_reduced) #577
card_results_reduced_df <- as.data.frame(card_results_reduced)
head(card_results_reduced_df)

write.csv(card_results_reduced_df, "card_reduced.csv", row.names = FALSE) 
write.table(card_results_reduced_df, "card_reduced.tab", sep = "\t", quote = FALSE, row.names = FALSE) 
