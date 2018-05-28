# -----------------------------------------------------------------------
#
#   Code to reproduce results from Franks et al. (2018)
#   Feature Specific Quantile Normalization Enables Cross-Platform
#   Classification of Molecular Subtypes using Gene Expression Data
#
#   Contact: Jennifer Franks - jennifer.m.franks.gr@dartmouth.edu
#
#   Microarray and RNA-seq data from TCGA BRCA datasets
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
library(devtools)
library(TDM)
library(preprocessCore)
library(Biobase)
library(huge)
library(scales)
library(gridExtra)
library(viridis)
library(scales)
library(plotrix)
library(gPCA)
library(utils)
library(caret)
library(kernlab)
library(glmnet)
library(randomForest)


setwd("Replace me!!")

load("BRCA_UNC_microarray.rda")
array <- t(UNC.AgilentG4502A)
t.rnaseq1 <- readPCLtoMatrix("TCGA_BRCA_RNA-seq_RPKM.PCL") #RPKM recalculated from RSEM

rnaseq1 <- t(t.rnaseq1)
sub <- read.delim("BRCA_547_PAM50_Subtypes.txt")


rnaseq.copy <- rnaseq1
rnaseq1 <- log2(rnaseq1)
is.na(rnaseq1) <- sapply(rnaseq1, is.infinite)
rnaseq1 <- rnaseq1[, colSums(is.na(rnaseq1)) < 50] # only select columns with enough data to impute


#only get the samples in common so we can use true gold standard
comm_samples <- list()
comm_samples$array <- rownames(array)
comm_samples$seq <- rownames(rnaseq1)
comm_samples$sub <- sub$Sample
common.symbols <- Reduce(intersect, comm_samples)

array <- array[which(rownames(array) %in% common.symbols),]
array <- array[order(rownames(array)),]
rnaseq1 <- rnaseq1[which(rownames(rnaseq1) %in% common.symbols),]
rnaseq1 <- rnaseq1[order(rownames(rnaseq1)),]
rnaseq.copy <- rnaseq.copy[which(rownames(rnaseq.copy) %in% common.symbols),]
rnaseq.copy <- rnaseq.copy[order(rownames(rnaseq.copy)),]

comm_genes <- list()
comm_genes$array <- colnames(array)
comm_genes$seq <- colnames(rnaseq1)
common.symbols2 <- Reduce(intersect, comm_genes)

# need to get the matrices to match by rownames and colnames according to what is in common
array <- array[,which(colnames(array) %in% common.symbols2)]
array <- array[,order(colnames(array))]
rnaseq1 <- rnaseq1[,which(colnames(rnaseq1) %in% common.symbols2)]
rnaseq1 <- rnaseq1[,order(colnames(rnaseq1))]
rnaseq.copy <- rnaseq.copy[,which(colnames(rnaseq.copy) %in% common.symbols2)]
rnaseq.copy <- rnaseq.copy[,order(colnames(rnaseq.copy))]


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


# impute missing values for RNA-seq and Microarray
array <- knnImputation(array)
colmed <- apply(array,2,median) #median center by gene (remember transposed!)
array =  sweep(array,2,colmed,"-")

t.rnaseq <- knnImputation(t(rnaseq1), k = 5)
rnaseq <- t(t.rnaseq)


#scale the dataset by feature -----------------

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


#test the batch bias before normalizaing ------------------
together <- rbind(rnaseq, array)
together1 <- rbind(matrix(1, nrow=nrow(rnaseq)), matrix(2, nrow=nrow(array)))
library(gPCA)
bd.results <- gPCA.batchdetect(together, together1)
print(paste("p-value =", bd.results$p.val))
PCplot(bd.results, ug="unguided", type="comp", npcs=3)
ks.test(rnaseq, array)

# then quantile normalize ---------------
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

# and also TDM / NPN / QN it --------------

colnames(array_scaled) = colnames(array)

#TDM 
rnaseq_tdm <- tdm_transform(ref_data = data.table(cbind(gene=colnames(array), t(array))), target_data = data.table(cbind(gene = colnames(rnaseq), t(rnaseq))))
rnaseq_tdm_scaled <- tdm_transform(ref_data = data.table(cbind(gene=colnames(array_scaled), t(array_scaled))), target_data = data.table(cbind(gene = colnames(rnaseq), t(rnaseq))))

rnaseq_tdm <- t(data.matrix(rnaseq_tdm))[2:540,]
rnaseq_tdm_scaled <- scale_bygene(t(data.matrix(rnaseq_tdm_scaled))[2:540,])


# QN
rnaseq_qn <- normalize.quantiles.use.target(x = rnaseq, target = as.vector(array), copy = TRUE)
rnaseq_qn_scaled <- normalize.quantiles.use.target(x = rnaseq, target = as.vector(array_scaled), copy = TRUE)

# NPN 
rnaseq_npn <- t(huge.npn(t(rnaseq)))
rnaseq_npn_scaled <- scale_bygene(rnaseq_npn)

# and test the batch bias with the new normalized datasets --------------
library(gPCA)
together <- rbind(rnaseq_qn, array) # fill in appropriate dataset for rnaseq (and microarray if scaled)
together1 <- rbind(matrix(1, nrow=nrow(rnaseq)), matrix(2, nrow=nrow(array)))
bd.results <- gPCA.batchdetect(together, together1)
print(paste("p-value =", bd.results$p.val))

PCplot(bd.results, ug="unguided", type="comp", npcs=3)

 # create the models ---------------------------------------------------------

fitControl <- trainControl(method = "repeatedcv", number = 3, repeats = 10)

set.seed(7)
modelSvm <- train(array,as.factor(array.key$subtype), method = "svmLinear",  trControl = fitControl)
modelSvm_scaled <- train(array_scaled,as.factor(array.key$subtype), method = "svmLinear",  trControl = fitControl)

set.seed(7)
modelRF <- train(array,as.factor(array.key$subtype), method = "rf",  trControl = fitControl)
modelRF_scaled <- train(array_scaled,as.factor(array.key$subtype), method = "rf",  trControl = fitControl)

set.seed(7)
modelGLM <- train(array,as.factor(array.key$subtype), method = "glmnet",  trControl = fitControl)
modelGLM_scaled <- train(array_scaled,as.factor(array.key$subtype), method = "glmnet",  trControl = fitControl)

# test final datasets ---------------------------------------------------

colnames(rnaseq_tdm_scaled) = colnames(rnaseq)

save.image("FSQN_fullcode_BRCA_results.Rdata")

# SVM unscaled
confusionMatrix(predict(modelSvm$finalModel, array),array.key$subtype)
confusionMatrix(predict(modelSvm$finalModel, rnaseq),rnaseq.key$subtype)
confusionMatrix(predict(modelSvm$finalModel, rnaseq_fsqn),rnaseq.key$subtype)
confusionMatrix(predict(modelSvm$finalModel, rnaseq_qn),rnaseq.key$subtype)
confusionMatrix(predict(modelSvm$finalModel, rnaseq_tdm),rnaseq.key$subtype)
confusionMatrix(predict(modelSvm$finalModel, rnaseq_npn),rnaseq.key$subtype)

# SVM scaled
confusionMatrix(predict(modelSvm_scaled$finalModel, array_scaled),array.key$subtype)
confusionMatrix(predict(modelSvm_scaled$finalModel, rnaseq_scaled),rnaseq.key$subtype)
confusionMatrix(predict(modelSvm_scaled$finalModel, rnaseq_fsqn_scaled),rnaseq.key$subtype)
confusionMatrix(predict(modelSvm_scaled$finalModel, rnaseq_qn_scaled),rnaseq.key$subtype)
confusionMatrix(predict(modelSvm_scaled$finalModel, rnaseq_tdm_scaled),rnaseq.key$subtype)
confusionMatrix(predict(modelSvm_scaled$finalModel, rnaseq_npn_scaled),rnaseq.key$subtype)

# random forest unscaled
confusionMatrix(predict(modelRF$finalModel, array),array.key$subtype)
confusionMatrix(predict(modelRF$finalModel, rnaseq), rnaseq.key$subtype)
confusionMatrix(predict(modelRF$finalModel, rnaseq_fsqn),rnaseq.key$subtype)
confusionMatrix(predict(modelRF$finalModel, rnaseq_qn),rnaseq.key$subtype)
confusionMatrix(predict(modelRF$finalModel, rnaseq_tdm),rnaseq.key$subtype)
confusionMatrix(predict(modelRF$finalModel, rnaseq_npn),rnaseq.key$subtype)

# random forest scaled
confusionMatrix(predict(modelRF_scaled$finalModel, array_scaled),array.key$subtype)
confusionMatrix(predict(modelRF_scaled$finalModel, rnaseq_scaled),rnaseq.key$subtype)
confusionMatrix(predict(modelRF_scaled$finalModel, rnaseq_fsqn_scaled),rnaseq.key$subtype)
confusionMatrix(predict(modelRF_scaled$finalModel, rnaseq_qn_scaled),rnaseq.key$subtype)
confusionMatrix(predict(modelRF_scaled$finalModel, rnaseq_tdm_scaled),rnaseq.key$subtype)
confusionMatrix(predict(modelRF_scaled$finalModel, rnaseq_npn_scaled),rnaseq.key$subtype)

# glmnet unscaled
confusionMatrix(predict(modelGLM$finalModel, array, s = modelGLM$finalModel$lambdaOpt, type = "class"),array.key$subtype)
confusionMatrix(predict(modelGLM$finalModel, rnaseq, s = modelGLM$finalModel$lambdaOpt, type = "class"), rnaseq.key$subtype)
confusionMatrix(predict(modelGLM$finalModel, rnaseq_fsqn, s = modelGLM$finalModel$lambdaOpt, type = "class"),rnaseq.key$subtype)
confusionMatrix(predict(modelGLM$finalModel, rnaseq_qn, s = modelGLM$finalModel$lambdaOpt, type = "class"),rnaseq.key$subtype)
confusionMatrix(predict(modelGLM$finalModel, rnaseq_tdm, s = modelGLM$finalModel$lambdaOpt, type = "class"),rnaseq.key$subtype)
confusionMatrix(predict(modelGLM$finalModel, rnaseq_npn, s = modelGLM$finalModel$lambdaOpt, type = "class"),rnaseq.key$subtype)

# glmnet scaled
confusionMatrix(predict(modelGLM_scaled$finalModel, array_scaled, s = modelGLM_scaled$finalModel$lambdaOpt, type = "class"),array.key$subtype)
confusionMatrix(predict(modelGLM_scaled$finalModel, rnaseq_scaled, s = modelGLM_scaled$finalModel$lambdaOpt, type = "class"),rnaseq.key$subtype)
confusionMatrix(predict(modelGLM_scaled$finalModel, rnaseq_fsqn_scaled, s = modelGLM_scaled$finalModel$lambdaOpt, type = "class"),rnaseq.key$subtype)
confusionMatrix(predict(modelGLM_scaled$finalModel, rnaseq_qn_scaled, s = modelGLM_scaled$finalModel$lambdaOpt, type = "class"),rnaseq.key$subtype)
confusionMatrix(predict(modelGLM_scaled$finalModel, rnaseq_tdm_scaled, s = modelGLM_scaled$finalModel$lambdaOpt, type = "class"),rnaseq.key$subtype)
confusionMatrix(predict(modelGLM_scaled$finalModel, rnaseq_npn_scaled, s = modelGLM_scaled$finalModel$lambdaOpt, type = "class"),rnaseq.key$subtype)




# Figure 3 ----------------------------------------
df = data.frame(Microarray=sample(array, 100000),
                FSQN=sample(rnaseq_fsqn, 100000),
                QN=sample(rnaseq_qn,100000),
                NPN=sample(rnaseq_npn, 100000),
                TDM = sample(rnaseq_tdm, 100000),
                LOG2 = sample(rnaseq, 1000000)
)


df2 = data.frame(Microarray = array[,100], 
                 FSQN=rnaseq_fsqn[,100],
                 QN=rnaseq_qn[,100], 
                 NPN=rnaseq_npn[,100],
                 TDM=rnaseq_tdm[,100],
                 LOG2=rnaseq[,100])

colnames(array_scaled)[100]
df.m <- reshape2::melt(df, id.vars = NULL)
df2.m <- reshape2::melt(df2, id.vars = NULL)

data_summary <- function(x) {
    m <- mean(x)
    ymin <- m-sd(x)
    ymax <- m+sd(x)
    return(c(y=m,ymin=ymin,ymax=ymax))
}


col1 = c("#404040", "#0d0887", "#6a00a8", "#b12a90", "#e16462", "#fca636", "#f0f921")

p1 <- ggplot(df.m, aes(x = variable, y = value)) + geom_violin(aes(fill = variable), scale = "width")+ylim(c(-10,10))+
    theme_classic() + scale_fill_manual(values = col1)+ stat_summary(fun.data=data_summary) + ggtitle("Unscaled")
p2 <- ggplot(df2.m, aes(x = variable, y = value)) + geom_violin(aes(fill = variable), scale = "width")+
    theme_classic() + scale_fill_manual(values = col1)+ stat_summary(fun.data=data_summary) + ggtitle("Unscaled") +ylim(c(-10,10))


grid.arrange(p1,p2)


ggplot(df.m, aes(value)) + geom_density(aes(colour = variable))+
    theme_classic() + scale_fill_manual(values = col1) + ggtitle("Unscaled")


# Figure 4 --------------------------------------------
df = data.frame(Microarray=sample(array_scaled, 100000),
                FSQN=sample(rnaseq_fsqn_scaled, 100000),
                QN=sample(rnaseq_qn_scaled,100000),
                NPN=sample(rnaseq_npn_scaled, 100000),
                TDM = sample(rnaseq_tdm_scaled, 100000),
                LOG2 = sample(rnaseq_scaled, 1000000)
)


df2 = data.frame(Microarray = array_scaled[,100], 
                 FSQN=rnaseq_fsqn_scaled[,100],
                 QN=rnaseq_qn_scaled[,100], 
                 NPN=rnaseq_npn_scaled[,100],
                 TDM=rnaseq_tdm_scaled[,100],
                 LOG2=rnaseq_scaled[,100])


df.m <- reshape2::melt(df, id.vars = NULL)
df2.m <- reshape2::melt(df2, id.vars = NULL)

data_summary <- function(x) {
    m <- mean(x)
    ymin <- m-sd(x)
    ymax <- m+sd(x)
    return(c(y=m,ymin=ymin,ymax=ymax))
}


col1 = c("#404040", "#0d0887",  "#b12a90", "#e16462", "#fca636", "#f0f921")

p1 <- ggplot(df.m, aes(x = variable, y = value)) + geom_violin(aes(fill = variable), scale = "width")+
    theme_classic() + scale_fill_manual(values = col1)+ stat_summary(fun.data=data_summary) + ggtitle("Scaled")
p2 <- ggplot(df2.m, aes(x = variable, y = value)) + geom_violin(aes(fill = variable), scale = "width")+
    theme_classic() + scale_fill_manual(values = col1)+ stat_summary(fun.data=data_summary) +ggtitle("Scaled")
grid.arrange(p1,p2)



#----- Figure 2 ------------

df = data.frame(BRCA1 = array[,"BRCA1"],
                BRCA1_rnaseq = rnaseq[,"BRCA1"],
                CDH1 = array[,"CDH1"],
                CDH1_rnaseq = rnaseq[,"CDH1"],
                ERBB2=array[,"ERBB2"], 
                ERBB2_rnaseq=rnaseq[,"ERBB2"],
                TP53=array[,"TP53"],
                TP52_rnaseq=array[,"TP53"],
                PTEN=array[,"PTEN"],
                PTEN_rnaseq=rnaseq[,"PTEN"])

df.m <- reshape2::melt(df, id.vars = NULL)

data_summary <- function(x) {
    m <- mean(x)
    ymin <- m-sd(x)
    ymax <- m+sd(x)
    return(c(y=m,ymin=ymin,ymax=ymax))
}

col2 = c("#1b0c42", "#1b0c42", "#1b0c42", "#1b0c42", "#1b0c42")

ggplot(df.m, aes(x = variable, y = value)) + geom_violin(aes(fill = variable), scale = "width")+ scale_y_continuous(limits = c(-10, 10))+
    theme_classic() + stat_summary(fun.data=data_summary) + ggtitle("Unscaled")


df = data.frame("Microarray"=sample(array, 100000),
                "RNA-seq" = sample(rnaseq, 100000))

df.m <- reshape2::melt(df, id.vars = NULL)

data_summary <- function(x) {
    m <- mean(x)
    ymin <- m-sd(x)
    ymax <- m+sd(x)
    return(c(y=m,ymin=ymin,ymax=ymax))
}

ggplot(df.m, aes(x = variable, y = value)) + geom_violin(aes(fill = variable), scale = "width")+scale_y_continuous(limits = c(-10, 10))+
    theme_classic() + stat_summary(fun.data=data_summary) + ggtitle("Unscaled")

# SUPPLEMENTARY FIGURES 5 and 6 -------------------

#supp figure 5
colnames(rnaseq_qn) = colnames(rnaseq_tdm) = colnames(array)
df = data.frame(Microarray_gene = array[,"BRCA1"], 
                Log2RPKM=rnaseq[,"BRCA1"], 
                FSQN=rnaseq_fsqn[,"BRCA1"], 
                QN=rnaseq_qn[,"BRCA1"], 
                NPN=rnaseq_npn[,"BRCA1"],
                TDM=rnaseq_tdm[,"BRCA1"])

df2 = data.frame(Microarray_gene = array[,"ERBB2"], 
                 Log2RPKM=rnaseq[,"ERBB2"], 
                 FSQN=rnaseq_fsqn[,"ERBB2"], 
                 QN=rnaseq_qn[,"ERBB2"], 
                 NPN=rnaseq_npn[,"ERBB2"],
                 TDM=rnaseq_tdm[,"ERBB2"])

df3 = data.frame(Microarray_gene = array[,"CDH1"], 
                 Log2RPKM=rnaseq[,"CDH1"], 
                 FSQN=rnaseq_fsqn[,"CDH1"], 
                 QN=rnaseq_qn[,"CDH1"], 
                 NPN=rnaseq_npn[,"CDH1"],
                 TDM=rnaseq_tdm[,"CDH1"])


df.m <- reshape2::melt(df, id.vars = NULL)
df2.m <- reshape2::melt(df2, id.vars = NULL)
df3.m <- reshape2::melt(df3, id.vars = NULL)

data_summary <- function(x) {
    m <- mean(x)
    ymin <- m-sd(x)
    ymax <- m+sd(x)
    return(c(y=m,ymin=ymin,ymax=ymax))
}
library(viridis)
col1 = inferno(10)[2:10]
col1 = c("#1b0c42", "#781c69", "#a52c60", "#cf4446", "#ed6925", "#ffca53")

p1 <- ggplot(df.m, aes(x = variable, y = value)) + geom_violin(aes(fill = variable), scale = "width")+ylab("Expression Values")+
    theme_classic() + scale_fill_manual(values = col1)+ stat_summary(fun.data=data_summary) + ggtitle("Unscaled - ERBB2")
p2 <- ggplot(df2.m, aes(x = variable, y = value)) + geom_violin(aes(fill = variable), scale = "width")+ylab("Expression Values")+
    theme_classic() + scale_fill_manual(values = col1)+ stat_summary(fun.data=data_summary) + ggtitle("Unscaled - BRCA1")
p3 <- ggplot(df3.m, aes(x = variable, y = value)) + geom_violin(aes(fill = variable), scale = "width")+ ylab("Expression Values")+
    theme_classic() + scale_fill_manual(values = col1)+ stat_summary(fun.data=data_summary) + ggtitle("Unscaled - CDH1")
library(gridExtra)
grid.arrange(p1,p2,p3)

#supp figure 6
colnames(rnaseq_qn_scaled) = colnames(rnaseq_tdm_scaled) = colnames(array)
df = data.frame(Microarray_gene = array_scaled[,"BRCA1"], 
                Log2RPKM=rnaseq_scaled[,"BRCA1"], 
                FSQN=rnaseq_fsqn_scaled[,"BRCA1"], 
                QN=rnaseq_qn_scaled[,"BRCA1"], 
                NPN=rnaseq_npn_scaled[,"BRCA1"],
                TDM=rnaseq_tdm_scaled[,"BRCA1"])

df2 = data.frame(Microarray_gene = array_scaled[,"ERBB2"], 
                 Log2RPKM=rnaseq_scaled[,"ERBB2"], 
                 FSQN=rnaseq_fsqn_scaled[,"ERBB2"], 
                 QN=rnaseq_qn_scaled[,"ERBB2"], 
                 NPN=rnaseq_npn_scaled[,"ERBB2"],
                 TDM=rnaseq_tdm_scaled[,"ERBB2"])

df3 = data.frame(Microarray_gene = array_scaled[,"CDH1"], 
                 Log2RPKM=rnaseq_scaled[,"CDH1"], 
                 FSQN=rnaseq_fsqn_scaled[,"CDH1"], 
                 QN=rnaseq_qn_scaled[,"CDH1"], 
                 NPN=rnaseq_npn_scaled[,"CDH1"],
                 TDM=rnaseq_tdm_scaled[,"CDH1"])


df.m <- reshape2::melt(df, id.vars = NULL)
df2.m <- reshape2::melt(df2, id.vars = NULL)
df3.m <- reshape2::melt(df3, id.vars = NULL)

data_summary <- function(x) {
    m <- mean(x)
    ymin <- m-sd(x)
    ymax <- m+sd(x)
    return(c(y=m,ymin=ymin,ymax=ymax))
}

col1 = c("#1b0c42", "#781c69", "#a52c60", "#cf4446", "#ed6925", "#ffca53")

p1 <- ggplot(df.m, aes(x = variable, y = value)) + geom_violin(aes(fill = variable), scale = "width")+ylab("Expression Values")+
    theme_classic() + scale_fill_manual(values = col1)+ stat_summary(fun.data=data_summary) + ggtitle("Scaled - ERBB2")
p2 <- ggplot(df2.m, aes(x = variable, y = value)) + geom_violin(aes(fill = variable), scale = "width")+ylab("Expression Values")+
    theme_classic() + scale_fill_manual(values = col1)+ stat_summary(fun.data=data_summary) + ggtitle("Scaled - BRCA1")
p3 <- ggplot(df3.m, aes(x = variable, y = value)) + geom_violin(aes(fill = variable), scale = "width")+ ylab("Expression Values")+
    theme_classic() + scale_fill_manual(values = col1)+ stat_summary(fun.data=data_summary) + ggtitle("Scaled - CDH1")

grid.arrange(p1,p2,p3)


