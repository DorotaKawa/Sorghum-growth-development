
library("goseq")
library(GenomicRanges)
library(rtracklayer)
library(Rsamtools)

setwd("/Users/Dorota/Dropbox/UC_DAVIS/PROMISE/SRN39 paper/RNAseq/")

############## Create lists with genes from each cluster
DEList1 <- read.csv("Clustering/All-Clusters-top75var.csv", header = T,stringsAsFactors = F)
DEList6 <- DEList5 <- DEList4 <- DEList3 <- DEList2 <- DEList1
DEList <- list(DEList1, DEList2, DEList3, DEList4, DEList5, DEList6)
names(DEList) <- c("C1", "C2","C3", "C4", "C5", "C6")
setwd("/Users/Dorota/Dropbox/UC_DAVIS/PROMISE/SRN39 paper/RNAseq/Clustering/")
filesUniq <- list.files(".",pattern = "cluster")
uniqueList <- lapply(filesUniq,function(x){read.csv(x, header = TRUE,stringsAsFactors = F)})
names(uniqueList) <- c("C1", "C2","C3", "C4", "C5", "C6")
sapply(uniqueList,nrow)
head(DEList[[1]])


############### Create seqlengths from gff file
setwd("/Users/Dorota/Dropbox/UC_DAVIS/PROMISE/POP#2/RNA/DEG/4.SQR2/GOseq")
GFFfile = "Sbicolor_454_v3.1.1.gff.txt"
GFF <- import.gff(GFFfile,version="3",feature.type="gene")
grl <- reduce(split(GFF, mcols(GFF)$Name))
reducedGTF <- unlist(grl, use.names=T)
mcols(reducedGTF)$Name <- rep(names(grl), elementNROWS(grl))
reducedGTF
AllLengths <- width(reducedGTF)
names(AllLengths) <- mcols(reducedGTF)$Name
head(reducedGTF)
head(AllLengths)
go.goseq <- read.delim("Sbicolor_454_v3.1.1.annotation_go.txt",header = F)
############### Enrichments
GOList <- list()
for (each in names(DEList)){
  print (each) 
  tmp <- DEList[[each]]
  rownames(tmp) <- tmp$ID
  assayed.genes <- rownames(tmp)
  de.genes <- uniqueList[[each]]$ID
  length(assayed.genes)
  length(de.genes)
  gene.vector=as.integer(assayed.genes%in%de.genes)
  names(gene.vector)=assayed.genes
  head(gene.vector)
  table(gene.vector)
  gene.vector2 <- gene.vector
  GeneLengths <- AllLengths[names(gene.vector2)]
  cbind(GeneLengths,gene.vector)
  
  all.genes <- names(gene.vector)
  pwf=nullp(gene.vector, "AGI", id=all.genes, bias.data=GeneLengths)
  head(pwf)
  
  # 
  go.all <- goseq(pwf, "ITAG3.1", gene2cat=go.goseq)
  head(go.all)
  dim(go.all)
  table(go.all$over_represented_pvalue < 0.05)
  
  go.sign <- go.all[go.all$over_represented_pvalue < 0.05,]
  
  #View(go.sign)
  GOList[[each]] <- go.all
  
}

sapply(GOList, dim)

signGOList <- lapply(GOList, function(x){  x[x$over_represented_pvalue < 0.05,] })

lapply(signGOList, head)

source("metaFunctions_forNetworkAnalysis.R")

list2env(signGOList,envir=.GlobalEnv)

setwd("/Users/Dorota/Dropbox/UC_DAVIS/PROMISE/SRN39 paper/RNAseq/Clustering/")

write.csv(C1, "cluster1_GOs.csv")
write.csv(C2, "cluster2_GOs.csv")
write.csv(C3, "cluster3_GOs.csv")
write.csv(C4, "cluster4_GOs.csv")
write.csv(C5, "cluster5_GOs.csv")
write.csv(C6, "cluster6_GOs.csv")

