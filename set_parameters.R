#parameters setting

session_name <- "TT8_RNAseq"

#input files
count.matrix <- "in/count_matrix.csv"
gene.length.column <- "Length"
metadata.filename <- "in/metadata.csv"
traits.filename <- "in/traits.csv"

#filtering parameters
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