# -----------------------------------------------------------------------
#
#   Code to reproduce results from Franks et al. (2018)
#   Feature Specific Quantile Normalization Enables Cross-Platform
#   Classification of Molecular Subtypes using Gene Expression Data
#
#   Contact: Jennifer Franks - jennifer.m.franks.gr@dartmouth.edu
#
#   Microarray and RNA-seq data from TCGA CRC datasets
#
# -----------------------------------------------------------------------


readPCLtoMatrix <- function(filename){
    
    # reads PCL file, returns data frame
    # assumes first column contains unique identifiers
    
    data <- read.delim(filename, header=FALSE, skip = 2, stringsAsFactors=FALSE)
    hl <- readLines(filename, 1)
    hl <- strsplit(hl, "\t")
    colnames(data) <- hl[[1]]
    data <- data[,c(-2,-3)]
    rownames(data) <- data[,1]
    data <- data[,-1]
    
    return(data)
    
}


library(RCurl)
library(preprocessCore)
library(DMwR)
library(scales)
library(plotrix)
library(devtools)
library(TDM)
library(preprocessCore)
library(Biobase)
library(huge)
library(utils)
library(caret)
library(devtools)
library(kernlab)


#set the working directory
setwd("Replace me!")
array <- readPCLtoMatrix("tcgacrc_microarray.imputed.collapsed.PCL")

rnaseq1 <- read.delim("tcgacrc_expression.tsv", header = T)
sub <- read.delim("cms_labels_temp.txt")

rownames(rnaseq1) <- rnaseq1$feature
rnaseq1<- rnaseq1[ , !(names(rnaseq1) == "feature")]

#change to m by n 
array <- t(array)
rnaseq1 <- t(rnaseq1)


#only get the samples in common
comm_samples <- list()
comm_samples$array <- rownames(array)
comm_samples$seq <- rownames(rnaseq1)
comm_samples$sub <- sub$Sample[!is.na(sub$PAM50)]
common.symbols <- Reduce(intersect, comm_samples)

array <- array[which(rownames(array) %in% common.symbols),]
array <- array[order(rownames(array)),]
rnaseq1 <- rnaseq1[which(rownames(rnaseq1) %in% common.symbols),]
rnaseq1 <- rnaseq1[order(rownames(rnaseq1)),]


comm_genes <- list()
comm_genes$array <- colnames(array)
comm_genes$seq <- colnames(rnaseq1)
common.symbols2 <- Reduce(intersect, comm_genes)

#get matrices to match by common rownames and colnames 
array <- array[,which(colnames(array) %in% common.symbols2)]
array <- array[,order(colnames(array))]
rnaseq <- rnaseq1[,which(colnames(rnaseq1) %in% common.symbols2)]
rnaseq <- rnaseq1[,order(colnames(rnaseq1))]


array.key <- list()
array.key$id <- rownames(array)
for (i in 1:length(array.key$id)){
  array.key$subtype[i] <- sub$PAM50[which(sub$Sample == array.key$id[i])]
}

rnaseq.key <- list()
rnaseq.key$id <- rownames(rnaseq1)
for (i in 1:length(rnaseq.key$id)){
  rnaseq.key$subtype[i] <- sub$PAM50[which(sub$Sample == rnaseq.key$id[i])]
}


# median center the microarray dataset by gene
rowmed <- apply(array,2,median)
array <- sweep(array,2,rowmed,"-")


#scale the dataset by gene/feature -----------------

#Takes M x N matrix!
scale_bygene <- function(to_scale){
  data.scaled <- matrix(0, nrow = nrow(to_scale), ncol = ncol(to_scale))
  for (i in 1:ncol(to_scale)){
    gene.to.scale <- to_scale[,i]
    result <- rescale(gene.to.scale, c(0,1))
    data.scaled[,i] <- result
  }
  
  rownames(data.scaled) = rownames(to_scale)
  colnames(data.scaled) = colnames(to_scale)
  
  return(data.scaled)
}

array_scaled <- scale_bygene(array)
rnaseq_scaled <- scale_bygene(rnaseq)

colnames(array_scaled) = colnames(array)
colnames(rnaseq_scaled) = colnames(rnaseq)

# FSQN code ---------------
fsqn <- function(to_normalize, target_dist){
  data.qn <- matrix(0, nrow = nrow(to_normalize), ncol = ncol(to_normalize))
  
  for (i in 1:ncol(to_normalize)){
    gene.to.normalize <- to_normalize[,i]
    target.gene.dist <- target_dist[,i]
    result <- normalize.quantiles.use.target(x = as.matrix(gene.to.normalize), 
                                             target = target.gene.dist, 
                                             copy = TRUE)
    data.qn[,i] <- result
  }
  
  rownames(data.qn) = rownames(to_normalize)
  colnames(data.qn) = colnames(to_normalize)
  
  return(data.qn)
}

rnaseq_fsqn_scaled <- fsqn(rnaseq, array_scaled)
rnaseq_fsqn <- fsqn(rnaseq, array)


# TDM / NPN / QN  --------------

#TDM 
rnaseq_tdm <- tdm_transform(ref_data = data.table(cbind(gene=colnames(array), t(array))), 
                            target_data = data.table(cbind(gene = colnames(rnaseq), t(rnaseq))))

rnaseq_tdm_scaled <- tdm_transform(ref_data = data.table(cbind(gene=colnames(array_scaled), t(array_scaled))), 
                                   target_data = data.table(cbind(gene = colnames(rnaseq), t(rnaseq))))

rnaseq_tdm <- t(data.matrix(rnaseq_tdm))[1:length(rownames(array))+1,]
rnaseq_tdm_scaled <- t(data.matrix(rnaseq_tdm_scaled))[1:length(rownames(array))+1,]
rnaseq_tdm_scaled <- scale_bygene(rnaseq_tdm_scaled)

colnames(rnaseq_tdm_scaled) = colnames(rnaseq)
colnames(rnaseq_tdm) = colnames(rnaseq)


# QN (using target distribution)
rnaseq_qn <- normalize.quantiles.use.target(x = rnaseq, target = as.vector(array), copy = TRUE)
rnaseq_qn_scaled <- normalize.quantiles.use.target(x = rnaseq, target = as.vector(array_scaled), copy = TRUE)

# NPN 
rnaseq_npn <- t(huge.npn(t(rnaseq)))
rnaseq_npn_scaled <- scale_bygene(rnaseq_npn)

# Train the models ---------------------------------------------------------

fitControl <- trainControl(method = "repeatedcv", number = 3, repeats = 10)

set.seed(7)
modelSvm <- train(array,as.factor(array.key$subtype), method = "svmLinear",  trControl = fitControl)
modelSvm_scaled <- train(array_scaled,as.factor(array.key$subtype), method = "svmLinear",  trControl = fitControl)

modelRF <- train(array,as.factor(array.key$subtype), method = "rf", trControl = fitControl)
modelRF_scaled <- train(array_scaled,as.factor(array.key$subtype), method = "rf", trControl = fitControl)

modelGLM <- train(array,as.factor(array.key$subtype), method = "glmnet",  trControl = fitControl)
modelGLM_scaled <- train(array_scaled,as.factor(array.key$subtype), method = "glmnet",  trControl = fitControl)

save.image("FSQN_CRC_models.Rdata")

# test final datasets ---------------------------------------------------

# svm
confusionMatrix(predict(modelSvm$finalModel, array),array.key$subtype)
confusionMatrix(predict(modelSvm$finalModel, rnaseq),rnaseq.key$subtype)
confusionMatrix(predict(modelSvm$finalModel, rnaseq_fsqn),rnaseq.key$subtype)
confusionMatrix(predict(modelSvm$finalModel, rnaseq_qn),rnaseq.key$subtype)
confusionMatrix(predict(modelSvm$finalModel, rnaseq_tdm),rnaseq.key$subtype)
confusionMatrix(predict(modelSvm$finalModel, rnaseq_npn),rnaseq.key$subtype)


confusionMatrix(predict(modelSvm_scaled$finalModel, array_scaled),array.key$subtype)
confusionMatrix(predict(modelSvm_scaled$finalModel, rnaseq_scaled),rnaseq.key$subtype)
confusionMatrix(predict(modelSvm_scaled$finalModel, rnaseq_fsqn_scaled),rnaseq.key$subtype)
confusionMatrix(predict(modelSvm_scaled$finalModel, rnaseq_qn_scaled),rnaseq.key$subtype)
confusionMatrix(predict(modelSvm_scaled$finalModel, rnaseq_tdmscaled),rnaseq.key$subtype)
confusionMatrix(predict(modelSvm_scaled$finalModel, rnaseq_npn_scaled),rnaseq.key$subtype)


# random forest 
confusionMatrix(predict(modelRF$finalModel, array),array.key$subtype)
confusionMatrix(predict(modelRF$finalModel, rnaseq), rnaseq.key$subtype)
confusionMatrix(predict(modelRF$finalModel, rnaseq_fsqn),rnaseq.key$subtype)
confusionMatrix(predict(modelRF$finalModel, rnaseq_qn),rnaseq.key$subtype)
confusionMatrix(predict(modelRF$finalModel, rnaseq_tdm),rnaseq.key$subtype)
confusionMatrix(predict(modelRF$finalModel, rnaseq_npn),rnaseq.key$subtype)


confusionMatrix(predict(modelRF_scaled$finalModel, array_scaled),array.key$subtype)
confusionMatrix(predict(modelRF_scaled$finalModel, rnaseq_scaled),rnaseq.key$subtype)
confusionMatrix(predict(modelRF_scaled$finalModel, rnaseq_fsqn_scaled),rnaseq.key$subtype)
confusionMatrix(predict(modelRF_scaled$finalModel, rnaseq_qn_scaled),rnaseq.key$subtype)
confusionMatrix(predict(modelRF_scaled$finalModel, rnaseq_tdm_scaled),rnaseq.key$subtype)
confusionMatrix(predict(modelRF_scaled$finalModel, rnaseq_npn_scaled),rnaseq.key$subtype)

# glmnet
confusionMatrix(predict(modelGLM$finalModel, array, s = modelGLM$finalModel$lambdaOpt, type = "class"),array.key$subtype)
confusionMatrix(predict(modelGLM$finalModel, rnaseq, s = modelGLM$finalModel$lambdaOpt, type = "class"), rnaseq.key$subtype)
confusionMatrix(predict(modelGLM$finalModel, rnaseq_fsqn, s = modelGLM$finalModel$lambdaOpt, type = "class"),rnaseq.key$subtype)
confusionMatrix(predict(modelGLM$finalModel, rnaseq_qn, s = modelGLM$finalModel$lambdaOpt, type = "class"),rnaseq.key$subtype)
confusionMatrix(predict(modelGLM$finalModel, rnaseq_tdm, s = modelGLM$finalModel$lambdaOpt, type = "class"),rnaseq.key$subtype)
confusionMatrix(predict(modelGLM$finalModel, rnaseq_npn, s = modelGLM$finalModel$lambdaOpt, type = "class"),rnaseq.key$subtype)


confusionMatrix(predict(modelGLM_scaled$finalModel, array_scaled, s = modelGLM_scaled$finalModel$lambdaOpt, type = "class"),array.key$subtype)
confusionMatrix(predict(modelGLM_scaled$finalModel, rnaseq_scaled, s = modelGLM_scaled$finalModel$lambdaOpt, type = "class"),rnaseq.key$subtype)
confusionMatrix(predict(modelGLM_scaled$finalModel, rnaseq_fsqn_scaled, s = modelGLM_scaled$finalModel$lambdaOpt, type = "class"),rnaseq.key$subtype)
confusionMatrix(predict(modelGLM_scaled$finalModel, rnaseq_qn_scaled, s = modelGLM_scaled$finalModel$lambdaOpt, type = "class"),rnaseq.key$subtype)
confusionMatrix(predict(modelGLM_scaled$finalModel, rnaseq_tdm_scaled, s = modelGLM_scaled$finalModel$lambdaOpt, type = "class"),rnaseq.key$subtype)
confusionMatrix(predict(modelGLM_scaled$finalModel, rnaseq_npn_scaled, s = modelGLM_scaled$finalModel$lambdaOpt, type = "class"),rnaseq.key$subtype)

save.image("FSQN_fullcode_CRC_results.Rdata")


