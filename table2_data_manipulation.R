
# Import files 
library(readxl) 
library(dplyr)
BacMet_gene_legend <- read_excel("BacMet_results_separate.xlsx", sheet = "gene_legend2") 
View(BacMet_gene_legend)
colnames(BacMet_gene_legend)

BacMet_gene_legend_columns <- BacMet_gene_legend %>% select(Gene_Alias, Gene_name, Compounds)

BacMet_gene_counts <- read_excel("BacMet_EXP_reduced_for_heatmap_final.xlsx", sheet = "gene_counts") 
View(BacMet_gene_counts)
colnames(BacMet_gene_counts)
BacMet_gene_counts <- BacMet_gene_counts[,1:3]

BacMet_table2 <- inner_join(BacMet_gene_counts, BacMet_gene_legend_columns, by = "Gene_Alias")
View(BacMet_table2)
write.csv(BacMet_table2, "BacMet_table2.csv")
