BiocManager::install("DESeq", version = "3.8")
BiocManager::install("gplots")

library("DESeq")
library(plyr)

##################### import counts table ############################
##
## I imported my unfiltered count table after using tximport. Your table can come from HTSeq-count or Rsubread. I think DESeq has built in functions for HTSeq-count data.
##

#txi.kallisto.counts = read.table(file='txi.kallisto.counts.txt', sep="\t", header=TRUE, row.names = 1)
count_217_219 <- read.table(file = "217_219_count_table.txt",header = TRUE, sep = '\t') 
names(count_217_219) <- c("ENSEMBL_ID", "count_217", "count_219")
count_217_219 <- count_217_219[ ,-c(4) ]
count_221 <- read.table(file = "221_count_table.txt",header = FALSE, sep = '\t') 
names(count_221) <- c("ENSEMBL_ID", "count_221")
count_222 <- read.table(file = "222_count_table.txt",header = FALSE, sep = '\t')
names(count_222) <- c("ENSEMBL_ID", "count_222")

## import file with ENSEMBL IDs annotated with gene names

withIDs <- read.table(file = 'hisat_geneSymbolsIDs.csv', header = TRUE, sep = ',')
names(withIDs) <- c("ENSEMBL_ID", "gene_name")

counts_joined <- join(count_217_219, count_221, by = 'ENSEMBL_ID', type = "inner")
counts_joined <- join(counts_joined, count_222, by = 'ENSEMBL_ID', type = "inner")
#counts_joined <- join(counts_joined, withIDs, by = 'ENSEMBL_ID', type = "inner")

hisat.countData <- sapply(counts_joined, as.integer, USE.NAMES = TRUE)
rownames(hisat.countData) <- counts_joined$ENSEMBL_ID
colnames(counts_joined)<- c("ENSEMBL_ID", "AgedSyn","AdultSyn","AdultCntrl", "AgedCntrl")
colnames(hisat.countData)<- c("ENSEMBL_ID", "AgedSyn","AdultSyn","AdultCntrl", "AgedCntrl")
hisat.countData <- hisat.countData[,-1]
hisat.countData <- hisat.countData[,-5]
hisat.countData <- hisat.countData[,-5]
head(hisat.countData)
hisat.countData[is.na(hisat.countData)] <- 0 

## in the end, my count data table looks like this.

#                       X217_Aged_Syn X219_Adult_Syn X221_Adult_Ctrl X222_Aged_Ctrl
# ENSMUSG00000000001.4           1757           1107            2342           2006
# ENSMUSG00000000003.13             0              0               0              0
# ENSMUSG00000000028.12            77             80             122             47
# ENSMUSG00000000031.13         47150          50924           18127            621
# ENSMUSG00000000037.14            41             35              61             30
# ENSMUSG00000000049.9            108             11              12              5

################## make the samples table as before ######################
##
## add another column with the contrast we did before (AgedSyn vs AdultSyn)
##

samples = data.frame(row.names = c("AgedSyn","AdultSyn","AdultCntrl", "AgedCntrl"),
                     "Condition" = c("Aged","Adult","Adult","Aged"),
                     "Sample" = c("Syn","Syn","Ctrl","Ctrl"))
samples$group = factor(paste0(samples$Condition, samples$Sample))
samples

########################## create the count dataset ############################
condition = samples$group
cds = newCountDataSet(hisat.countData, condition )


########################## pre-filtering #########################
keep <- rowSums(counts(cds)) >= 1
sum(keep)
length(rownames(dds))
fract_kept = sum(keep)/length(rownames(cds))
fract_kept
cds <- cds[keep,]

###################### perform differential expression ###############################
cds = estimateSizeFactors( cds )
sizeFactors(cds)
head( counts( cds, normalized=TRUE ) )

cds = estimateDispersions(  cds, method="blind", sharingMode="fit-only" )
str( fitInfo(cds) )
plotDispEsts( cds )

res = nbinomTest( cds, "AgedSyn", "AdultSyn" )
write.table(res,'resHisat.csv',sep = ',')

######################## how many DE genes did we get? ############################
head(res)
head(res[order(res$padj),])
addmargins( table( res_sig = res$padj < .05) )


############################ MA plots, heatmaps, and PCA plots ##############################
plotMA(res)

vsd = varianceStabilizingTransformation( cds )
library("RColorBrewer")
library("gplots")
select = order(rowMeans(counts(cds)), decreasing=TRUE)[1:30]
hmcol = colorRampPalette(brewer.pal(9, "GnBu"))(100)
heatmap.2(exprs(vsd)[select,], col = hmcol, trace="none", margins = c(8,13))

print(plotPCA(vsd, intgroup=c("condition")))

# Annotate ENSEMBL IDs to gene IDs
library("AnnotationDbi")
library("org.Mm.eg.db")
geneid = rownames(counts_joined)
geneid = unlist(lapply(geneid, function(x) {
  x <- strsplit(x,split='.',fixed=TRUE)[[1]][1]
  return(x)
}))
length(unique(geneid))==length(rownames(counts_joined)) # making sure that we're not losing any genes
rownames(counts_joined) = geneid
head(rownames(counts_joined))

#rownames(counts_joined) <- counts_joined[,1]

columns(org.Mm.eg.db)
counts_joined$symbol = mapIds(org.Mm.eg.db,
                        keys=rownames(counts_joined),
                        column="SYMBOL",
                        keytype="ENSEMBL",
                        multiVals="first")
counts_joined$entrez <- mapIds(org.Mm.eg.db,
                         keys=row.names(counts_joined),
                         column="ENTREZID",
                         keytype="ENSEMBL",
                         multiVals="first")
#resOrdered <- counts_joined[order(counts_joined$pvalue),]
#top100genes <- rownames(resOrdered[1:100,])

############# Immitating their heatmap
ache = match('Ache', counts_joined$symbol)
musk = match('Musk', counts_joined$symbol)
hdac4 = match('Hdac4', counts_joined$symbol)
lrp4 = match('Lrp4', counts_joined$symbol)
chrna1 = match('Chrna1', counts_joined$symbol)
myog = match('Myog', counts_joined$symbol)
myod = match('Myod1', counts_joined$symbol)
ppargc1a = match('Ppargc1a', counts_joined$symbol)
esrra = match('Esrra', counts_joined$symbol)
mfn2 = match('Mfn2', counts_joined$symbol)
opa1 = match('Opa1', counts_joined$symbol)
pgp = match('Pgp', counts_joined$symbol)
itgb1 = match('Itgb1', counts_joined$symbol)
gadd45a = match('Gadd45a', counts_joined$symbol)
mrpl15 = match('Mrpl15', counts_joined$symbol)
myl4 = match('Myl4', counts_joined$symbol)
fst = match('Fst', counts_joined$symbol)

keep20 <- c(ache,musk,hdac4,lrp4,chrna1,myog,myod,ppargc1a,esrra,mfn2,opa1,pgp,itgb1,gadd45a,mrpl15,myl4,fst)
keep20
cds <- cds[keep20,]

cds = estimateSizeFactors( cds )
sizeFactors(cds)
head( counts( cds, normalized=TRUE ) )

cds = estimateDispersions(  cds, method="blind", sharingMode="fit-only" )
str( fitInfo(cds) )
plotDispEsts( cds )

res = nbinomTest( cds, "AgedSyn", "AdultSyn" )
head(res)
head(res[order(res$padj),])
addmargins( table( res_sig = res$padj < .05) )

vsd = varianceStabilizingTransformation( cds )
library("RColorBrewer")
library("gplots")
hmcol = colorRampPalette(brewer.pal(9, "GnBu"))(100)
heatmap.2(exprs(vsd), col = hmcol, trace="none", margins = c(8,13))


##### Pathway Analysis
#if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
#BiocManager::install("gage", version = "3.8")
