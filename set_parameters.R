#parameters setting

session_name <- "TT8_RNAseq"

#input files
count.matrix <- "in/count_matrix.csv"
gene.length.column <- "Length"
metadata.filename <- "in/metadata.csv"
traits.filename <- "in/traits.csv"

#filtering parameters
excluded_samples <- c("E2") #samples to be excluded from the analysis
filtering_method <- "2"

n <- 10 #number of reads
m <- 3 #number of samples (used in method "2")
p <- 0.75 #proportion of samples (used in method "3")

#network parameters
MaxBS <- 20000 #only for blockwise analysis
network_type <- "signed"
TOMtype <- "signed"
mergeCutHeight <- 0.3
corType <- "pearson"
deepSplit <- 2
minModSize <- 30
rs <- 42

#categorical module-trait correlation parameters
#NULL automatically selects factor, character, and logical metadata columns.
categorical_variables <- NULL

#Optional named vector defining the reference level for selected variables.
#Variables not listed here use their existing first factor level.
categorical_reference_levels <- c(
  mat.gen = "Col-0",
  mat.ploidy = "2x",
  pat.ploidy = "2x",
  viability = "viable"
)

#Generate pairwise contrasts and, optionally, level-versus-all contrasts.
categorical_include_pairwise <- TRUE
categorical_include_vs_all <- TRUE
categorical_min_count <- 1
categorical_drop_redundant <- TRUE

#NULL retains every non-redundant generated contrast. To select or exclude
#contrasts, use the exact names printed by 4_module-trait_cor.qmd.
categorical_contrasts <- NULL
categorical_exclude_contrasts <- NULL
