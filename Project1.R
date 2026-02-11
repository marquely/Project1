# Project 1 template

# important -- this makes sure our runs are consistent and reproducible
set.seed(0)

file <- "TCGA_breast_cancer_LumA_vs_Basal_PAM50.tsv"
first10 <- c('NAT1','BIRC5','BAG1','BCL2','BLVRA','CCNB1','CCNE1','CDC6','CDC20','CDH3')
nfold=5
  
header <- scan(file, nlines = 1, sep="\t", what = character())
data <- read.table(file, skip = 2, header = FALSE, sep = "\t", quote = "", check.names=FALSE)

header[1] <- "gene_id"
names(data) <- header

header2 <- scan(file, skip = 1, nlines = 1, sep="\t", what = character())

LumA <- data[data$gene_id %in% first10,header2=='Luminal A']
Basal <- data[data$gene_id %in% first10,header2=='Basal-like']


# define function cross_valid so we can rerun the cross validataion with various parameters
#cross_validation <- function (nfold, alg="centroid") {

  # split each cancer type samples into nfold groups
  LumA_groups <- split(sample(colnames(LumA)), 1+(seq_along(colnames(LumA)) %% nfold))
  Basal_groups <- split(sample(colnames(Basal)), 1+(seq_along(colnames(Basal)) %% nfold))
  
  result <- array()
 
  # iterate from 1 to nfold groups -- to choose test group
  for (test_group in 1:nfold) {
   
    # return all samples in the chosen test group
    testLumA <- LumA[,colnames(LumA) %in% unlist(LumA_groups[test_group])]
    testBasal <- Basal[,colnames(Basal) %in% unlist(Basal_groups[test_group])]
   
    # return all samples *not* in the chosen test group 
    trainingLumA <- LumA[,!(colnames(LumA) %in% unlist(LumA_groups[test_group]))]
    trainingBasal <- Basal[,!(colnames(Basal) %in% unlist(Basal_groups[test_group]))]
   
    # compute centroid for each cancer type -- mean for each gene based on all samples
    # note -- rows are gene
    centroidLumA <- rowMeans(trainingLumA)
    centroidBasal <- rowMeans(trainingBasal)
   
    # For each sample in the test set decide whether it will be classified
    # distance from centroid Lum A: sum(abs(x-centroidLumA))
    # distance from centroid Basal: sum(abs(x-centroidBasal))
    # distance is a sum of distances over all genes 
    # misclassification if when the distance is greater from centroid associated with known result
    misclassifiedLumA <- sum(sapply(testLumA, function(x) { sum(abs(x-centroidLumA))>sum(abs(x-centroidBasal)) }))
    misclassifiedBasal <- sum(sapply(testBasal, function(x) { sum(abs(x-centroidLumA))<sum(abs(x-centroidBasal)) }))
    
    result[test_group] <- (misclassifiedLumA+misclassifiedBasal)/(ncol(testLumA)+ncol(testBasal))
 }
 
 c(mean(result), sd(result))
#}

x<-data.frame(three=cross_validation(3), five=cross_validation(5), ten=cross_validation(10))
rownames(x) <- c('mean','sd')
x

