library(DESeq2)
library(tximport)
library(dplyr)
library(ggplot2)
library(magrittr)
library(Biobase)
library(pheatmap)
library(RColorBrewer)
library(IHW)
library(tidyr)
# library(apeglm)
samples <- read.csv("samples.tsv", header = TRUE, sep = "\t")
exprnames <- do.call(paste, c(samples[c("Name", "Rep")], sep = "."))
filenames <- files <- file.path("results","STAR",paste(exprnames,"ReadsPerGene.out.tab",sep=""))

counts.files <- lapply( filenames, read.table, skip = 4 )
counts <- as.data.frame( sapply( counts.files, function(x) x[ , number ] ) )
ff <- gsub( "[.]ReadsPerGene[.]out[.]tab", "", ff )
ff <- gsub( "results/STAR/", "", ff )
colnames(counts) <- ff
row.names(counts) <- counts.files[[1]]$V1


#txi.kallisto <- tximport(files, type = "kallisto", txOut = TRUE)
head(txi.kallisto$counts)
head(txi.kallisto$abundance)
colnames(txi.kallisto$counts) <- exprnames
colnames(txi.kallisto$abundance) <- exprnames
write.csv(txi.kallisto$abundance, "reports/kallisto.TPM.csv")
write.csv(txi.kallisto$counts, "reports/kallisto.counts.csv")

oldsamples = samples
newsamples <- samples %>% unite("Loc_Sucep", Location:ImidaclopridStatus, na.rm = TRUE,
  remove = FALSE)
samples = newsamples

# DEseq2 analyses
Location = factor(samples$Location)
ImidaclopridStatus = factor(samples$ImidaclopridStatus)
sampleTable <- data.frame(suceptibility = ImidaclopridStatus, location = Location,
  Name = Name)
rownames(sampleTable) = exprnames

dds <- DESeqDataSetFromTximport(txi.kallisto, sampleTable, design = ~location + suceptibility)

# nrow(dds)
dds <- dds[rowSums(counts(dds)) > 1, ]
# nrow(dds)

dds <- estimateSizeFactors(dds)
vsd <- vst(dds, blind = FALSE)
rld <- rlog(dds, blind = FALSE)
# head(assay(vsd), 3)

dds <- DESeq(dds, fitType = "mean")
resultsNames(dds)
# dds <-
# lfcShrink(dds,coef='condition_Temecula_Susceptible_vs_GeneralBeale_Resistant',type='apeglm')
# res <- results(dds, contrast=c('condition','Mycelia','Spherule48H')) res <-
# results(dds, name='condition_Spherule48h_vs_Mycelia')
res <- results(dds, filterFun = ihw)
summary(res)
res <- subset(res, res$padj < 0.05)
res <- res[order(res$pvalue), ]
summary(res)


df <- bind_rows(as_tibble(log2(counts(dds, normalized = TRUE)[, 1:2] + 1)) %>% mutate(transformation = "log2(x + 1)"),
  as_tibble(assay(rld)[, 1:2]) %>% mutate(transformation = "rlog"), as_tibble(assay(vsd)[,
    1:2]) %>% mutate(transformation = "vst"))

colnames(df)[1:2] <- c("x", "y")

pdf("plots/RNASeq_kallisto.pdf")
ggplot(df, aes(x = x, y = y)) + geom_hex(bins = 80) + coord_fixed() + facet_grid(. ~
  transformation)

select <- order(rowMeans(counts(dds, normalized = TRUE)), decreasing = TRUE)[1:75]
df <- as.data.frame(colData(dds)[, c("location", "suceptibility")])
rownames(df) = exprnames
colnames(df) = c("Location", "Suceptibility")

pheatmap(assay(vsd)[select, ], cluster_rows = FALSE, show_rownames = TRUE, fontsize_row = 7,
  fontsize_col = 7, cluster_cols = FALSE, annotation_col = df, main = "VSD ordered")

topVar <- head(order(rowVars(assay(vsd)), decreasing = TRUE), 60)
mat <- assay(vsd)[topVar, ]
mat <- mat - rowMeans(mat)
pheatmap(mat, show_rownames = TRUE, fontsize_row = 7, fontsize_col = 7, cluster_cols = FALSE,
  annotation_col = df, main = "VSD - most different")

pheatmap(assay(rld)[select, ], cluster_rows = FALSE, show_rownames = TRUE, fontsize_row = 7,
  fontsize_col = 7, cluster_cols = FALSE, annotation_col = df, main = "RLD")

sampleDists <- dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix(sampleDists)
# rownames(sampleDistMatrix) <- paste(vsd$location, vsd$suceptibility,sep='-')
rownames(sampleDistMatrix) <- paste(vsd$location, vsd$suceptibility, sep = "-")
colors <- colorRampPalette(rev(brewer.pal(9, "Blues")))(255)
colnames(sampleDistMatrix) <- NULL

pheatmap(sampleDistMatrix, clustering_distance_rows = sampleDists, clustering_distance_cols = sampleDists,
  col = colors)

pcaData <- plotPCA(vsd, intgroup = c("suceptibility"), returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

ggplot(pcaData, aes(PC1, PC2, color = Location, shape = ImidaclopridStatus)) + geom_point(size = 3) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) + ylab(paste0("PC2: ", percentVar[2],
  "% variance")) + coord_fixed()

## res <- results(dds, contrast=c('condition','Mycelia','Spherule48H'))
## summary(res) plotMA(res,main='Mycelium vs Spherule48H') res <-
## subset(res,res$padj<0.05) res <- res[order(res$pvalue ),] summary(res)
## write.csv(res,'reports/kallisto.DESeq_Mycelium_Spherule48H.csv')


## res <- results(dds, contrast=c('condition','Mycelia','Spherule8D'))
## summary(res) plotMA(res,main='Mycelium vs Spherule8D') res <-
## subset(res,res$padj<0.05) res <- res[order(res$pvalue ),] summary(res)
## write.csv(res,'reports/kallisto.DESeq_Mycelium_Spherule8D.csv')

## res <- results(dds, contrast=c('condition','Spherule48H','Spherule8D'))
## summary(res) plotMA(res,main='Spherule48H vs Spherule8D') res <-
## subset(res,res$padj<0.05) res <- res[order(res$pvalue ),] summary(res)
## write.csv(res,'reports/kallisto.DESeq_Spherule48H_Spherule8D.csv')
