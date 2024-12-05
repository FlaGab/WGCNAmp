library(openxlsx)
library(dplyr)
library(tidyr)
library(ggplot2)
library(stringr)

GeneNames = read.xlsx("in/GeneNames.xlsx")

for (i in dir("out/DEGs", pattern = "DEGs_")){
  df = read.xlsx(paste0("DEGs/",i))
  df$Gene_name = coalesce(df$Gene_name, df$Gene_ID)
  assign(i, df)
}

GO_experiment = read.csv("out/GO/ivory_ALL_GO_AGI_simplified_enrichment_results.csv")
#print eriched GO for selection
GO_experiment$Description

AGIs = GO_experiment[9,"geneID"] %>% 
  str_split(pattern = "/") %>% 
  unlist()
AGIs

AGIs = read.xlsx("out/Oct_30_18.28_2024corType_pearson_mCH0.25_pt19deepSplit2_wgcna_wgcna_km_gene_to_module.xlsx") %>%
  filter(module == "yellow") %>%
  select(gene) %>% 
  as.vector() %>%
  unlist()

query = tibble(Gene_ID = AGIs) %>% left_join(GeneNames[,1:2])
query$Gene_name = coalesce(query$Gene_name, query$Gene_ID)
OUT=left_join(query, DEGs_2Xx4X.vs.2Xx2X_STAR.xlsx[,c(1,2,7)]) %>% 
  left_join(DEGs_4Xx2X.vs.2Xx2X_STAR.xlsx[,c(1,2,7)], suffix = c("-Pat.Exc","-Mat.Exc"), by = join_by("Gene_ID", "Gene_name")) %>%
  left_join(DEGs_tt8x2X.vs.2Xx2X_STAR.xlsx[,c(1,2,7)], by = join_by("Gene_ID", "Gene_name")) %>%
  left_join(DEGs_tt8x4X.vs.2Xx2X_STAR.xlsx[,c(1,2,7)], suffix = c("-tt8x2X","tt8x4X"), by = join_by("Gene_ID", "Gene_name"))

OUT[is.na(OUT)] <- 0
arrange(OUT, desc(OUT$`logFC-tt8x2X`))

# Convert data to long format
data_long <- pivot_longer(OUT, cols = colnames(OUT[,c(3:6)]), names_to = "Condition", values_to = "logFC")

# Determine the positions of vertical separator
x_breaks <- seq(0.5, length(unique(data_long$Gene_name))+ 0.5)

# Plot
ggplot(data_long, aes(x = Gene_name, y = logFC, fill = Condition)) + 
  geom_bar(stat = "identity", position = "dodge") +
  labs(title = "Differential Expression",
       x = "Gene Name",
       y = "Log2 Fold Change") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "top",
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank()) +
  scale_fill_manual(values = c("#FF9999", "#66CCFF", "#cce5cc", "#329932"), 
                    labels = c("4xX2x vs 2xX2X", "2xX4x vs 2xX2X", "tt8X2x vs 2xX2X", "tt8X4x vs 2xX2X")) +
  geom_vline(xintercept = x_breaks, color = "gray", linewidth = 1) +
  geom_hline(yintercept = 0, color = "black", linewidth = 0.5)

