source("src/GO_env_setup.R")

#import module-associated genes
gene_list <- list()
gene_list$bkg <- read.xlsx("out/Oct_30_18.28_2024corType_pearson_mCH0.25_pt19deepSplit2_wgcna_wgcna_km_gene_to_module.xlsx", sheet = 1)
gene_list$ivory <- read.xlsx("out/Oct_30_18.28_2024corType_pearson_mCH0.25_pt19deepSplit2_wgcna_wgcna_km_gene_to_module.xlsx", sheet = 2)
gene_list$darkgrey <- read.xlsx("out/Oct_30_18.28_2024corType_pearson_mCH0.25_pt19deepSplit2_wgcna_wgcna_km_gene_to_module.xlsx", sheet = 3)
gene_list$darkorange2 <- read.xlsx("out/Oct_30_18.28_2024corType_pearson_mCH0.25_pt19deepSplit2_wgcna_wgcna_km_gene_to_module.xlsx", sheet = 4)
gene_list$blue <- read.xlsx("out/Oct_30_18.28_2024corType_pearson_mCH0.25_pt19deepSplit2_wgcna_wgcna_km_gene_to_module.xlsx", sheet = 5)
gene_list$darkmagenta <- read.xlsx("out/Oct_30_18.28_2024corType_pearson_mCH0.25_pt19deepSplit2_wgcna_wgcna_km_gene_to_module.xlsx", sheet = 6)

#set parameters
keyType <- "TAIR"
ont <- "ALL"     # Options are "BP", "MF", "CC", and "ALL"
pAdjustMethod <- "BH"    # Benjamini-Hochberg adjustment method
pvalueCutoff <- 0.05
qvalueCutoff <- 0.2
universe <- gene_list[['bkg']]$gene

#parameters for simplification
simCutoff <- 0.7
simBy <- "p.adjust"
simSelect_fun <- min
simMeasure <- "Wang"
simSemData <- NULL

exp_list = names(gene_list[-1])

#exp_list = c("4xX2x.vs.2xX2x_DOWN", "4xX2x.vs.2xX4x_DOWN")

for (i in exp_list){
  
  #select your genes list
  title = paste0(i,"_",ont)
  print(paste("Processing:", title))
  gene_vec = as.vector(gene_list[[i]]$gene)
  
  # Perform GO enrichment analysis
  go_enrichment <- enrichGO(
    gene = gene_vec,
    OrgDb = org.At.tair.db,
    keyType = keyType,
    ont = ont,    
    pAdjustMethod = pAdjustMethod,
    universe = universe,
    pvalueCutoff = pvalueCutoff,
    qvalueCutoff = qvalueCutoff,
    readable = FALSE       # Convert gene IDs to gene names
  )
  
  # View the results
  #head(go_enrichment)
  
  # Check if enrichment results are valid
  if (is.null(go_enrichment) || nrow(go_enrichment) == 0) {
    print(paste("No enrichment results for:", title))
    next
  }
  
  #calculate pairwise similarity between terms
  go_enrichment_sim = simplify(x = go_enrichment,
                               cutoff = simCutoff,
                               by = simBy,
                               select_fun = min,
                               measure = simMeasure,
                               semData = simSemData)

  # Specify the file path for the PDF file
  pdf_file <- paste0("out/GO/", title, "_GO_plots.pdf")
  
  # Open the PDF device with a higher resolution
  pdf(pdf_file, width = 8, height = 12, onefile = TRUE, family = "Helvetica", pointsize = 12)
  
  # Visualization - Dotplot
  dp = dotplot(go_enrichment_sim, showCategory = 20, title = title)
  
  # Visualization - Barplot
  bp = barplot(go_enrichment_sim, showCategory = 20, title = title)
  
  # Visualization - category-gene-network plot
  np = cnetplot(go_enrichment_sim, showCategory = 5, title = title)
  
  print(dp)
  print(bp)
  print(np)
  # Close the PDF device
  dev.off()
  
  # Save the results to a file
  csv_file <- paste0("out/GO/", title, "_GO_AGI_simplified_enrichment_results.csv")
  print(paste("Saving results to:", csv_file))
  write.csv(as.data.frame(go_enrichment_sim), file = csv_file)
}

# Visualization - Dotplot
png_file_dot <- paste0("out/GO/", title, "_GO_dot.png")
print(paste("Saving dotplot to:", png_file_dot))
png(png_file_dot, width = 800, height = 600, )
dotplot(go_enrichment, showCategory = 20, title = title)
dev.off()

# Visualization - Barplot
png_file_bar <- paste0("out/GO/", title, "_GO_bar.png")
print(paste("Saving barplot to:", png_file_bar))
png(png_file_bar, width = 800, height = 600)
barplot(go_enrichment, showCategory = 20, title = title)
dev.off()

# Visualization - cnetplot (category-gene-network plot)
png_file_net <- paste0("out/GO", title, "_GO_net.png")
print(paste("Saving cnetplot to:", png_file_net))
png(png_file_net, width = 800, height = 600)
cnetplot(go_enrichment, showCategory = 5, title = title)
dev.off()

