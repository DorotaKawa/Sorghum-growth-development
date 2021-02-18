#RNA seq - differential expression

library(edgeR)
library(RColorBrewer)
library(limma)

setwd("/Users/dorota/Dropbox/UC_DAVIS/PROMISE/SRN39 paper/RNAseq/")

########################  Read counts sets
Counts <- read.csv("Counts/20190321_Counts_Genos.csv",
                        header = T,stringsAsFactors = F,row.names = 1)
dim(Counts)

########################  Filtering counts with 0 in either dataset
filterZero <- T
if(filterZero){
  cat ("Filtering tables independently - remove genes with any 0 counts")
  
  Counts <- Counts[!apply(Counts,1,function(x){
    any(x==0)
  }),]
} else {
  cat ("Not filtering here - genes with any 0 counts stay")
}

######################## Create the metadata table
meta <- data.frame(do.call("rbind",strsplit(colnames(Counts),"_")),
                        row.names = colnames(Counts),stringsAsFactors = F)

# Manually set the column names
colnames(meta) <- c("Genotype","Time","Soil","Treatment","R")
head(meta)

####################### Sanity check for counts and meta data compatibility

checkNames <- function(){
  if ((all(colnames(Counts) %in% rownames(meta)))){
    cat ("all is good, # of samples:", nrow(meta))
  } else {
    warning("Check sample names")
  }}

checkNames()

####################### Filter out samples identified as outliers/low quality libraries
doFilter <- TRUE
if (doFilter){
  cat ("Filtering low counts and outlier samples \n")
  toFilter <- c("SQR_3_S_CTR_REP4") # Outlier 
  Outliers <- which(rownames(meta) %in% toFilter)
  Counts <- Counts[,-Outliers]
  meta <- meta[-Outliers,]
} else { cat ("No filter \n") }

checkNames()

####################### Removing lowly-expressed genes
dge <- DGEList(counts=Counts,remove.zeros = T)
minSampleMin <- 3 
minCPM <- 1
isexpr <- rowSums(cpm(dge) > minCPM) >= minSample
dge <- dge[isexpr,,keep.lib.size = FALSE]

#######################  Calculate the normalization factors with TMM
dge <- calcNormFactors(dge)
dge.norm$samples$norm.factors

par(mfrow=c(1,2))
lcpm <- cpm(dge, log=TRUE)
boxplot(lcpm)
title(main="A. Example: Unnormalised data",ylab="Log-cpm")

lcpm2 <- cpm(dge.norm, log=TRUE)
boxplot(lcpm2)
title(main="B. Example: Normalised data",ylab="Log-cpm")

##################### MDS plots
meta$GT <- paste(meta$Genotype,meta$Time,sep="_")
meta$GT <- as.factor(meta$GT)

par(mfrow=c(1,3))
col.g <- meta$Genotype
levels(col.g) <-  brewer.pal(nlevels(col.g), "Set1")
col.gl <- as.character(col.g)
col.t <- meta$Time
levels(col.t) <-  brewer.pal(nlevels(col.t), "Set2")
col.t <- as.character(col.t)
col.sampl <- meta$GT
levels(col.sampl) <-  brewer.pal(nlevels(col.sampl), "Set3")
col.sampl <- as.character(col.sampl)

plotMDS(lcpm2, labels=meta$Genotype,col=col.g)
title(main="A. Genotype")
plotMDS(lcpm2, labels=meta$Time, col=col.t, dim=c(3,4))
title(main="B. Time")
plotMDS(lcpm2, labels=meta$GT,col=col.sampl)
title(main="C. Sample")

####################### Differential expression - design matrix
meta$Genotype <- as.factor(meta$Genotype)
meta$Genotype <- relevel(meta$Soil,ref = "SQR")

meta$Time <- as.factor(meta$Time)
meta$Time <- relevel(meta$Time,ref = "2")

design <- model.matrix(~Genotype*Time, data=meta)
colnames(design) <- gsub("Soil|Day|Name|Genotype|Treatment|:|-|/|Groups","",colnames(design))
head(design)

####################### Differential expression - transform with voom (and renormalize with quantile).
v <- voom(dge, design, plot = TRUE,normalize.method = "quantile")

vfit <- lmFit(v, design)
vfit <- contrasts.fit(vfit, contrasts=contr.matrix)
efit <- eBayes(vfit)
plotSA(efit, main="Final model: Mean-variance trend")

####################### Save Expression Data
cpmExpression <- cpm(dge)
voomExpression <- v$E
colnames(cpmExpression) <- rownames(meta)
colnames(voomExpression) <-  rownames(meta)
write.csv(cpmExpression, "DE/cpmExpression_Genotypes.csv")
save(cpmExpression,file = "DE/cpmExpression_Genotypes.RData")
save(voomExpression,file = "DE/voomExpression_Genotypes.RData")

##########################  DEG
v2 <- lmFit(v, design)
fit <- eBayes(v2)
colSums(abs(decideTests(fit)))

coefficientsTested <- colnames(fit$coefficients)[-1]

DEList <- lapply(coefficientsTested,function(x){
  DEresults <- topTable(fit, coef=x, adjust="BH",number = Inf,sort.by = "none")
  colnames(DEresults) <- paste(x,colnames(DEresults),sep = ".")
  return(DEresults)
})
names(DEList) <- coefficientsTested
names(DEList)

########################## DEG filter FDR < 0.05
significantDE <- lapply(DEList,function(x){ x[x[,grep("adj.P.Val",colnames(x))] < 0.05,] })
sapply(significantDE,nrow)
shortName <- "anova_genotype*time"
saveDElist <- paste0("DEList_",shortName,".RData")
save(file = saveDElist,DEList)
sapply(DEList,colnames)

names(significantDE)

write.table( significantDE[["SRN39"]],file = "Genotype.csv", sep = ",",quote = F,col.names = T)
write.table( significantDE[["Time3"]],file = "Time.csv", sep = ",",quote = F,col.names = T)
write.table( significantDE[["SRN39Time3"]],file = "GenotypexTime.csv", sep = ",",quote = F,col.names = T)


