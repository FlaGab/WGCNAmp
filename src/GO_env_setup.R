#GO enrichmant with clusterprofiler evironment setup

if (!("clusterProfiler" %in% installed.packages())) {
  # Install this package if it isn't installed yet
  BiocManager::install("clusterProfiler", update = FALSE)
}

if (!("org.At.tair.db" %in% installed.packages())) {
  # Install this package if it isn't installed yet
  BiocManager::install("org.At.tair.db", update = FALSE)
}

if (!("GO.db" %in% installed.packages())) {
  # Install this package if it isn't installed yet
  BiocManager::install("GO.db", update = FALSE)
}

if (!("enrichplot" %in% installed.packages())) {
  # Install this package if it isn't installed yet
  BiocManager::install("enrichplot", update = FALSE)
}

library(clusterProfiler)
library(org.At.tair.db)
library(AnnotationDbi)
library(GO.db)
library(enrichplot)

GO_enrichment <- function(gma,
                          modules,
                          session_name = session_name,
                          orgDb = org.At.tair.db,
                          keyType = "TAIR",
                          ont = "ALL",     # Options are "BP", "MF", "CC", and "ALL"
                          pAdjustMethod = "BH",    # Benjamini-Hochberg adjustment method
                          pvalueCutoff = 0.05,
                          qvalueCutoff = 0.2,
                          simCutoff = 0.7,
                          simBy = "p.adjust",
                          simSelect_fun = min,
                          simMeasure = "Wang",
                          simSemData = NULL){
  
  for (i in modules){
    gene_lists <- lapply(modules, function(i) gma[gma$module == i, 1])
    names(gene_lists) <- modules
  }
  
  #set background genes
  universe <- gma$gene
  
  for (i in modules){
    
    #select your genes list
    title <- paste0(i,"_",ont)
    print(paste("Processing:", title))
    gene_vec <- gene_lists[[i]]
    
    # Perform GO enrichment analysis
    go_enrichment <- enrichGO(
      gene = gene_vec,
      OrgDb = orgDb,
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
    go_enrichment_sim = clusterProfiler::simplify(x = go_enrichment,
                                 cutoff = simCutoff,
                                 by = simBy,
                                 select_fun = min,
                                 measure = simMeasure,
                                 semData = simSemData)
    
    #head(go_enrichment_sim)
    
    # Specify the file path for the PDF file
    pdf_file <- paste0(session_name, "/out/GO/", title, "_GO_plots.pdf")
    
    # Open the PDF device with a higher resolution
    pdf(pdf_file, width = 8, height = 12, onefile = TRUE, family = "Helvetica", pointsize = 12)
    
    # Visualization - Dotplot
    dp = dotplot(go_enrichment_sim, showCategory = 20, title = title)
    
    # Visualization - Barplot
    bp = barplot(go_enrichment_sim, showCategory = 20, title = title)
    
    # Visualization - category-gene-network plot
    #np = cnetplot(go_enrichment_sim, showCategory = 5, title = title, circular = T)
    
    print(dp)
    print(bp)
    #print(np)
    # Close the PDF device
    dev.off()
    
    # Save the results to a file
    csv_file <- paste0(session_name, "/out/GO/", title, "_GO_AGI_simplified_enrichment_results.csv")
    print(paste("Saving results to:", csv_file))
    write.csv(as.data.frame(go_enrichment_sim), file = csv_file)
  }
  
  # Visualization - Dotplot
  #png_file_dot <- paste0(session_name, "/out/GO/", title, "_GO_dot.png")
  #print(paste("Saving dotplot to:", png_file_dot))
  #png(png_file_dot, width = 800, height = 600, )
  #dotplot(go_enrichment, showCategory = 20, title = title)
  #dev.off()
  
  # Visualization - Barplot
  #png_file_bar <- paste0(session_name, "/out/GO/", title, "_GO_bar.png")
  #print(paste("Saving barplot to:", png_file_bar))
  #png(png_file_bar, width = 800, height = 600)
  #barplot(go_enrichment, showCategory = 20, title = title)
  #dev.off()
  
  # Visualization - cnetplot (category-gene-network plot)
  #png_file_net <- paste0(session_name, "/out/GO", title, "_GO_net.png")
  #print(paste("Saving cnetplot to:", png_file_net))
  #png(png_file_net, width = 800, height = 600)
  #cnetplot(go_enrichment_sim, showCategory = 5, title = title, circular = T)
  #dev.off()
  
}

KEGG_enrichment <- function(gma,
                          modules,
                          session_name = session_name,
                          keyType = "kegg",#one of "kegg", 'ncbi-geneid', 'ncbi-proteinid' and 'uniprot'
                          org = "ath", #supported organism listed in 'https://www.genome.jp/kegg/catalog/org_list.html'
                          pAdjustMethod = "BH",    # Benjamini-Hochberg adjustment method
                          pvalueCutoff = 0.05,
                          qvalueCutoff = 0.2,
                          simCutoff = 0.7,
                          simBy = "p.adjust",
                          simSelect_fun = min,
                          simMeasure = "Wang",
                          simSemData = NULL){

  
  for (i in modules){
    gene_lists <- lapply(modules, function(i) gma[gma$module == i, 1])
    names(gene_lists) <- modules
  }
  
  #set background genes
  universe <- gma$gene
  
  for (i in modules){
    
    #select your genes list
    title <- paste0(i,"_KEGG")
    print(paste("Processing:", title))
    gene_vec <- gene_lists[[i]]
    
    # Perform GO enrichment analysis
    go_enrichment <- enrichKEGG(
      gene = gene_vec,
      organism = org,
      keyType = keyType,
      pAdjustMethod = pAdjustMethod,
      universe = universe,
      pvalueCutoff = pvalueCutoff,
      qvalueCutoff = qvalueCutoff
    )
    
    # View the results
    #head(go_enrichment)
    
    # Check if enrichment results are valid
    if (is.null(go_enrichment) || nrow(go_enrichment) == 0) {
      print(paste("No enrichment results for:", title))
      next
    }
    
    # Specify the file path for the PDF file
    pdf_file <- paste0(session_name, "/out/GO/", title, "_plots.pdf")
    
    # Open the PDF device with a higher resolution
    pdf(pdf_file, width = 8, height = 12, onefile = TRUE, family = "Helvetica", pointsize = 12)
    
    # Visualization - Dotplot
    dp = dotplot(go_enrichment, showCategory = 20, title = title)
    
    # Visualization - Barplot
    bp = barplot(go_enrichment, showCategory = 20, title = title)
    
    # Visualization - category-gene-network plot
    #np = cnetplot(go_enrichment_sim, showCategory = 5, title = title, circular = T)
    
    print(dp)
    print(bp)
    #print(np)
    # Close the PDF device
    dev.off()
    
    # Save the results to a file
    csv_file <- paste0(session_name, "/out/GO/", title, "_enrichment_results.csv")
    print(paste("Saving results to:", csv_file))
    write.csv(as.data.frame(go_enrichment), file = csv_file)
  }
}
