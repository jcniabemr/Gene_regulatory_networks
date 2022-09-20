#Part1 
#Run DESeq2 with mycoprotein count files reps uncollapsed 
library(DESeq2)
library(tximport)
setwd("C:/Users/john.connell/Documents/Bioinformatics_2022/projects/Mycoprotein/RNA_seq")

#Create tx2gene table 
gene.t <- paste0("g", 1:13216, ".t", 1)
gene.t <- as.data.frame(gene.t)
gene <- paste0("g", 1:13216)
gene <- as.data.frame(gene)
tx2gene <- cbind(gene.t, gene)
colnames(tx2gene) <- c("Name", "Name")

#Import data
txi.reps <- tximport(paste(list.dirs("./", full.names=T,recursive=F),"/quant.sf",sep=""),type="salmon",tx2gene=tx2gene,txOut=T)
mysamples <- list.dirs("./",full.names=F,recursive=F)
txi.genes <- summarizeToGene(txi.reps,tx2gene)
invisible(sapply(seq(1,3), function(i) {colnames(txi.genes[[i]])<<-mysamples}))

#Create experiment design table 
designtable <- mysamples
designtable <- as.data.frame(designtable)
designtable$isolate_time <- substr(designtable$designtable,1,3)
designtable$isolate <- substr(designtable$designtable,1,2)

condition <- function(isolate){
  if(isolate == "WT")
    return ("control")
  else 
    return ("experimental")
}

designtable$condition <- lapply(X = designtable$isolate, FUN = condition)
colnames(designtable) <- c("sample", "isolate_time", "isolate", "condition")

#Set up experiment 
colData <- designtable[order(designtable$sample),]
colData$Group <- paste0(colData$isolate_time, colData$condition)

#Define the GLM parameters 
design <- ~ Group
dds <- DESeqDataSetFromTximport(txi.genes,colData,design)

#Set rowsums threshold 
keep <- rowSums(counts(dds)) >= 50
dds <- dds[keep,]

#Library normalisation
dds <- estimateSizeFactors(dds)

#Run Deseq
dds <- DESeq(dds)
resultsNames(dds)

#Adjusted p value cut off 
alpha <- 0.05 

#Create contrast of results 
res= results(dds, alpha=alpha,contrast=c("Group","011experimental","WT1control"))
sig.res <- subset(res,padj<=alpha)
sig.res <- sig.res[order(sig.res$padj),]
write.table(sig.res,"C1_vs_WT_TP1_ALL.txt", sep="\t", na="", quote=FALSE)

res= results(dds, alpha=alpha,contrast=c("Group","021experimental","WT1control"))
sig.res <- subset(res,padj<=alpha)
sig.res <- sig.res[order(sig.res$padj),]
write.table(sig.res,"C2_vs_WT_TP1_ALL.txt", sep="\t", na="", quote=FALSE)

res= results(dds, alpha=alpha,contrast=c("Group","031experimental","WT1control"))
sig.res <- subset(res,padj<=alpha)
sig.res <- sig.res[order(sig.res$padj),]
write.table(sig.res,"C3_vs_WT_TP1_ALL.txt", sep="\t", na="", quote=FALSE)

res= results(dds, alpha=alpha,contrast=c("Group","041experimental","WT1control"))
sig.res <- subset(res,padj<=alpha)
sig.res <- sig.res[order(sig.res$padj),]
write.table(sig.res,"C4_vs_WT_TP1_ALL.txt", sep="\t", na="", quote=FALSE)

res= results(dds, alpha=alpha,contrast=c("Group","051experimental","WT1control"))
sig.res <- subset(res,padj<=alpha)
sig.res <- sig.res[order(sig.res$padj),]
write.table(sig.res,"C5_vs_WT_TP1_ALL.txt", sep="\t", na="", quote=FALSE)

res= results(dds, alpha=alpha,contrast=c("Group","061experimental","WT1control"))
sig.res <- subset(res,padj<=alpha)
sig.res <- sig.res[order(sig.res$padj),]
write.table(sig.res,"C6_vs_WT_TP1_ALL.txt", sep="\t", na="", quote=FALSE)

res= results(dds, alpha=alpha,contrast=c("Group","071experimental","WT1control"))
sig.res <- subset(res,padj<=alpha)
sig.res <- sig.res[order(sig.res$padj),]
write.table(sig.res,"C7_vs_WT_TP1_ALL.txt", sep="\t", na="", quote=FALSE)

res= results(dds, alpha=alpha,contrast=c("Group","081experimental","WT1control"))
sig.res <- subset(res,padj<=alpha)
sig.res <- sig.res[order(sig.res$padj),]
write.table(sig.res,"C8_vs_WT_TP1_ALL.txt", sep="\t", na="", quote=FALSE)

res= results(dds, alpha=alpha,contrast=c("Group","091experimental","WT1control"))
sig.res <- subset(res,padj<=alpha)
sig.res <- sig.res[order(sig.res$padj),]
write.table(sig.res,"C9_vs_WT_TP1_ALL.txt", sep="\t", na="", quote=FALSE)

res= results(dds, alpha=alpha,contrast=c("Group","101experimental","WT1control"))
sig.res <- subset(res,padj<=alpha)
sig.res <- sig.res[order(sig.res$padj),]
write.table(sig.res,"C10_vs_WT_TP1_ALL.txt", sep="\t", na="", quote=FALSE)

res= results(dds, alpha=alpha,contrast=c("Group","111experimental","WT1control"))
sig.res <- subset(res,padj<=alpha)
sig.res <- sig.res[order(sig.res$padj),]
write.table(sig.res,"C11_vs_WT_TP1_ALL.txt", sep="\t", na="", quote=FALSE)

res= results(dds, alpha=alpha,contrast=c("Group","121experimental","WT1control"))
sig.res <- subset(res,padj<=alpha)
sig.res <- sig.res[order(sig.res$padj),]
write.table(sig.res,"C12_vs_WT_TP1_ALL.txt", sep="\t", na="", quote=FALSE)

res= results(dds, alpha=alpha,contrast=c("Group","131experimental","WT1control"))
sig.res <- subset(res,padj<=alpha)
sig.res <- sig.res[order(sig.res$padj),]
write.table(sig.res,"C13_vs_WT_TP1_ALL.txt", sep="\t", na="", quote=FALSE)

res= results(dds, alpha=alpha,contrast=c("Group","141experimental","WT1control"))
sig.res <- subset(res,padj<=alpha)
sig.res <- sig.res[order(sig.res$padj),]
write.table(sig.res,"C14_vs_WT_TP1_ALL.txt", sep="\t", na="", quote=FALSE)

res= results(dds, alpha=alpha,contrast=c("Group","151experimental","WT1control"))
sig.res <- subset(res,padj<=alpha)
sig.res <- sig.res[order(sig.res$padj),]
write.table(sig.res,"C15_vs_WT_TP1_ALL.txt", sep="\t", na="", quote=FALSE)

res= results(dds, alpha=alpha,contrast=c("Group","161experimental","WT1control"))
sig.res <- subset(res,padj<=alpha)
sig.res <- sig.res[order(sig.res$padj),]
write.table(sig.res,"C16_vs_WT_TP1_ALL.txt", sep="\t", na="", quote=FALSE)

res= results(dds, alpha=alpha,contrast=c("Group","171experimental","WT1control"))
sig.res <- subset(res,padj<=alpha)
sig.res <- sig.res[order(sig.res$padj),]
write.table(sig.res,"C17_vs_WT_TP1_ALL.txt", sep="\t", na="", quote=FALSE)

res= results(dds, alpha=alpha,contrast=c("Group","181experimental","WT1control"))
sig.res <- subset(res,padj<=alpha)
sig.res <- sig.res[order(sig.res$padj),]
write.table(sig.res,"C182_vs_WT_TP1_ALL.txt", sep="\t", na="", quote=FALSE)

res= results(dds, alpha=alpha,contrast=c("Group","191experimental","WT1control"))
sig.res <- subset(res,padj<=alpha)
sig.res <- sig.res[order(sig.res$padj),]
write.table(sig.res,"C19_vs_WT_TP1_ALL.txt", sep="\t", na="", quote=FALSE)


library(dplyr)
library(tidyr)
library(pheatmap)


b <- list ()

files <- list.files(pattern = "\\.txt$")

for (input_file in files){
  d <- read.table (paste0 (input_file), header = T)
  d$gene <- row.names (d)
  d$file <- input_file
  b [[length (b) + 1]] <- d [, c("log2FoldChange", "gene", "file")]
}

test<-c()
for (i in c(1:19)){
  test<-rbind(test,b[[i]])
  
}

test2<-test %>%
  pivot_wider(names_from = file , values_from = log2FoldChange)

test3<-as.data.frame(test2)

# plot
row.names (test3) <- test3[,1]
test3 <- test3 [,2:ncol (test3)]

test4<-data.matrix(test3)
test4[is.na(test4)] = 0

paletteLength <- 50
myBreaks <- c(seq(min(test4), 0, length.out=ceiling(paletteLength/2) + 1), 
              seq(max(test4)/paletteLength, max(test4), length.out=floor(paletteLength/2)))
testmap<-pheatmap(test4, breaks = myBreaks, color=colorRampPalette(c("darkblue", "white", "red"))(50), main = "TP1 all data") 

save_pheatmap_pdf <- function(x, filename, width=10, height=18) {
  stopifnot(!missing(x))
  stopifnot(!missing(filename))
  pdf(filename, width=width, height=height)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}
save_pheatmap_pdf(testmap, "tp1_all.pdf")

gene_of_interest <- c("g3465")
gene_data <- test4[rownames(test4) %in% gene_of_interest,]
gene_data <- as.data.frame(gene_data)
gene_data <- t(gene_data)

paletteLength <- 5
myBreaks <- c(seq(min(gene_data), 0, length.out=ceiling(paletteLength/2) + 1), 
              seq(max(gene_data)/paletteLength, max(gene_data), length.out=floor(paletteLength/2)))
testmap<-pheatmap(gene_data, cluster_cols = FALSE, cluster_rows = FALSE, breaks = myBreaks, color=colorRampPalette(c("darkblue", "white", "red"))(5), main = "TP1 all data") 

save_pheatmap_pdf <- function(x, filename, width=10, height=2.5) {
  stopifnot(!missing(x))
  stopifnot(!missing(filename))
  pdf(filename, width=width, height=height)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}
save_pheatmap_pdf(testmap, "3465.pdf")




