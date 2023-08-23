# Original heatmap - yes or no ####

library(pheatmap)
library(readxl) 

virulence_mbandaka <- read_excel("abricate_resfinder_all.xlsx", sheet = 3)
# Make sure all NAs are 0 
virulence_mbandaka

??pheatmap
virulence_mbandaka_cp <- virulence_mbandaka 
virulence_mbandaka_cp <- as.data.frame(virulence_mbandaka_cp)
rownames(virulence_mbandaka_cp) <- virulence_mbandaka_cp$Sample
virulence_mbandaka_cp <- virulence_mbandaka_cp[,-1] # Remove term.name as a column 
View(virulence_mbandaka_cp) 
rownames(virulence_mbandaka_cp) 
virulence_mbandaka_cp[is.na(virulence_mbandaka_cp)] <- as.double("NA")

# Two colours for 0 and 1 
library(RColorBrewer) 
RColorBrewer::brewer.pal.info
RColorBrewer::display.brewer.all()
virulence_cols <- colorRampPalette(brewer.pal(11, "RdYlGn") )(100) # or PiYG
cols <- brewer.pal(11, "RdYlGn")
cols2 <- brewer.pal(11, "PiYG")
cols3 <- brewer.pal(11, "YlOrRd")

vfdb_heatmap <- pheatmap(virulence_mbandaka_cp, col = cols, na_col = "grey", 
                         cluster_cols = FALSE, cluster_rows = FALSE, angle_col = 315)

cellwidth = 10, fontsize = 12

tiff("vfdb_heatmap.tiff", height = 30, width = 60, units = 'cm', 
     compression = "lzw", res = 600) 
vfdb_heatmap
dev.off() 

ncol(virulence_mbandaka_cp) #94

virulence1 <- virulence_mbandaka_cp[,1:50]
head(virulence1)

virulence2 <- virulence_mbandaka_cp[,50:94]
head(virulence2)

vfdb_heatmap1 <- pheatmap(virulence1, col = cols, na_col = "grey", 
                         cluster_cols = FALSE, cluster_rows = FALSE, angle_col = 315)

tiff("vfdb_heatmap1.tiff", height = 30, width = 30, units = 'cm', 
     compression = "lzw", res = 600) 
vfdb_heatmap1
dev.off() 

vfdb_heatmap2 <- pheatmap(virulence2, col = cols, na_col = "grey", 
                          cluster_cols = FALSE, cluster_rows = FALSE, angle_col = 315)

tiff("vfdb_heatmap2.tiff", height = 30, width = 30, units = 'cm', 
     compression = "lzw", res = 600) 
vfdb_heatmap2
dev.off() 


# AMR results reduced - % identity #### 


View(card_results_reduced) 

colnames(card_results_reduced) 
#[1] "x_file"       "sequence"     "start"        "end"          "strand"       "gene"         "coverage"     "coverage_map"
#[9] "gaps"         "x_coverage"   "x_identity"   "database"     "accession"    "product"      "resistance"  

library(tidyr)
card_results_reduced_wide <- card_results_reduced %>% select(x_file,gene,x_identity) %>% 
  pivot_wider(names_from = gene, values_from = x_identity)
View(card_results_reduced_wide) 
card_results_reduced_wide_df <- apply(card_results_reduced_wide,2,as.character)
write.csv(card_results_reduced_wide_df, "card_results_reduced_wide.csv")

card_for_heatmap <- read.csv("card_results_reduced_wide.csv", na.strings = "NULL", row.names = "x_file")
head(card_for_heatmap)
# These have the external strains, vfdb results do not 

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
newcolnames_card <- lapply(
  colnames(card_for_heatmap),
  function(x) bquote(italic(.(x))))
newcolnames_card

library(pheatmap)
?pheatmap
card_identity_heatmap <- pheatmap(card_for_heatmap, col = bacmet_cols2, na_col = "grey", 
                                  cluster_cols = FALSE, cluster_rows = TRUE, angle_col = 45, breaks = c(95,96,97,98,99,100), labels_col = as.expression(newcolnames_card)) 
card_identity_heatmap 

card_identity_heatmap2 <- pheatmap(card_for_heatmap, col = bacmet_cols3, na_col = "grey", 
                                   cluster_cols = FALSE, cluster_rows = TRUE, angle_col = 45, breaks = c(95,96,97,98,99,100), labels_col = as.expression(newcolnames_card))

card_identity_heatmap2
#cellwidth = 10, fontsize = 12

tiff("card_reduced_heatmap_identity.tiff", height = 30, width = 15, units = 'cm', 
     compression = "lzw", res = 600) 
card_identity_heatmap
dev.off() 

# And try Cairo 
library(Cairo)
Cairo(file="card_reduced_heatmap_identity.png", 
      type="png",
      units="in", 
      width=6, 
      height=12, 
      pointsize=12, 
      dpi=300)
card_identity_heatmap
dev.off()

# Virulence results reduced - % identity #### 
View(virulence_results_reduced) 

colnames(virulence_results_reduced)

library(tidyr)
virulence_results_reduced_wide <- virulence_results_reduced %>% select(number_file,gene,percent_identity) %>% 
  pivot_wider(names_from = gene, values_from = percent_identity)
View(virulence_results_reduced_wide) 
virulence_results_reduced_wide_df <- apply(virulence_results_reduced_wide,2,as.character)
write.csv(virulence_results_reduced_wide_df, "virulence_results_reduced_wide.csv")

vfdb_for_heatmap <- read.csv("virulence_results_reduced_wide.csv", na.strings = "NULL", row.names = "number_file")
head(vfdb_for_heatmap)

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
newcolnames_vfdb <- lapply(
  colnames(vfdb_for_heatmap),
  function(x) bquote(italic(.(x))))
newcolnames_vfdb

library(pheatmap)
?pheatmap
vfdb_identity_heatmap <- pheatmap(vfdb_for_heatmap, col = bacmet_cols2, na_col = "grey", 
                                   cluster_cols = FALSE, cluster_rows = TRUE, angle_col = 45, breaks = c(95,96,97,98,99,100), labels_col = as.expression(newcolnames_vfdb)) 
vfdb_identity_heatmap 

vfdb_identity_heatmap2 <- pheatmap(vfdb_for_heatmap, col = bacmet_cols3, na_col = "grey", 
                                    cluster_cols = FALSE, cluster_rows = TRUE, angle_col = 45, breaks = c(95,96,97,98,99,100), labels_col = as.expression(newcolnames_vfdb))

vfdb_identity_heatmap2
#cellwidth = 10, fontsize = 12

tiff("vfdb_reduced_heatmap_identity2.tiff", height = 30, width = 60, units = 'cm', 
     compression = "lzw", res = 600) 
vfdb_identity_heatmap
dev.off() 

# And try Cairo 
library(Cairo)
Cairo(file="vfdb_reduced_heatmap_identity2.png", 
      type="png",
      units="in", 
      width=22, 
      height=12, 
      pointsize=12, 
      dpi=300)
vfdb_identity_heatmap
dev.off()

