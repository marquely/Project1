# important -- this makes sure our runs are consistent and reproducible
set.seed(0)

file <- "TCGA_breast_cancer_ERstatus_allGenes.txt"
nfold <- 5
sd_threashold <- 4
top_num <- 10
  
header <- scan(file, nlines = 1, sep="\t", what = character())
data <- read.table(file, skip = 2, header = FALSE, sep = "\t", quote = "", check.names=FALSE)

header[1] <- "gene_id"
names(data) <- header

header2 <- scan(file, skip = 1, nlines = 1, sep="\t", what = character())

# cleanup - remove genes with sd < 1
# compute sd for each gene
data_sd<-sapply(seq(nrow(data)), function(x) { as.numeric(sd(data[x,-1])) })

# add gene names to the sd list
data_sd_names<-cbind(data.frame(data_sd),data[,1])

# create an "include" list of all those genes where sd > threshold
include_list <- data_sd_names[data_sd_names[,1]>sd_threashold,2]
sample_classes <- header2[-1]
Positive <- data[data$gene_id %in% include_list,header2=='Positive']
Negative <- data[data$gene_id %in% include_list,header2=='Negative']

# define function cross_valid so we can rerun the cross validation with various parameters
cross_validation <- function (nfold, alg="centroid") {

  Positive_groups <- split(sample(colnames(Positive)), 1+(seq_along(colnames(Positive)) %% nfold))
  Negative_groups <- split(sample(colnames(Negative)), 1+(seq_along(colnames(Negative)) %% nfold))
  
  result <- array()
  
  for (test_group in 1:nfold) {
    
    testA <- Positive[,colnames(Positive) %in% unlist(Positive_groups[test_group])]
    testB <- Negative[,colnames(Negative) %in% unlist(Negative_groups[test_group])]
    
    trainingA <- Positive[,!(colnames(Positive) %in% unlist(Positive_groups[test_group]))]
    trainingB <- Negative[,!(colnames(Negative) %in% unlist(Negative_groups[test_group]))]
    
    # Feature selection -- 
    
    # compute t-statistic for each row
    training_t_stat<-data.frame(sapply(seq(nrow(trainingA)), function(x) { abs(as.numeric(t.test(trainingA[x,], trainingB[x,])$statistic)) }))
    
    # add gene id column
    training_t_stat_geneid<-cbind(training_t_stat,rownames(trainingA))
    colnames(training_t_stat_geneid) <- c('t','id')
    
    # pick top 50 based on t-statistic
    selected_genes <- head(training_t_stat_geneid[order(-training_t_stat_geneid$t),],n=top_num)[,2]
    
    # narrow down the list of genes based on t-statistic
    testA <- testA[rownames(testA) %in% selected_genes,]
    testB <- testB[rownames(testB) %in% selected_genes,]
    trainingA <- trainingA[rownames(trainingA) %in% selected_genes,]
    trainingB <- trainingB[rownames(trainingB) %in% selected_genes,]
    
    # Centroid model
    if (alg=="centroid") {
    centroidA <- rowMeans(trainingA)
    centroidB <- rowMeans(trainingB)
    
    misclassifiedA <- sum(sapply(testA, function(x) { sqrt(sum((x-centroidA)^2))-sqrt(sum((x-centroidB)^2))>0 }))
    misclassifiedB <- sum(sapply(testB, function(x) { sqrt(sum((x-centroidA)^2))-sqrt(sum((x-centroidB)^2))<0 }))
    
    result[test_group] <- (misclassifiedA+misclassifiedB)/(ncol(testA)+ncol(testB))
    }
    
    # General linear model
    if (alg=="GLM") {
      # Build training data frame 
      trainingC <- rbind( cbind(data.frame(t(trainingA)), cancer = 1),   cbind(data.frame(t(trainingB)), cancer = 0) )
      
      # Build test data set 
      test_C <- rbind( data.frame(t(testA)), 
                        data.frame(t(testB)))
      
      # Provide true labels 
      truth <- c(rep(1, ncol(testA)), rep(0, ncol(testB)))
      
      # Create logistic regression model
      model <- glm(cancer ~ ., data = trainingC, family = binomial)
      
      # Prediction model using regression
      p <- predict(model, newdata = test_C, type = "response")
      preds <- ifelse(p < 0.5, 0, 1)
      
      # Calculate false positive and false negative predictions
      misclassifiedC <- sum(preds[truth == 1] == 0)
      misclassifiedD <- sum(preds[truth == 0] == 1)
    
    # Calculate fold-level error rate
    result[test_group] <- (misclassifiedC + misclassifiedD) / (ncol(testA) + ncol(testB))
    }
  }
  return(list(mean = mean(result, na.rm = TRUE), sd = sd(result, na.rm = TRUE)))
}

centroid_res <- cross_validation(5, alg="centroid")
GLM_res <- cross_validation(5, alg="GLM")

comparison <- data.frame(
  Algorithm = c("Centroid", "GLM"),
  MeanError = c(centroid_res$mean, GLM_res$mean),
  SD = c(centroid_res$sd, GLM_res$sd)
)

knitr::kable(comparison)
