setwd("/PATH/TO/WD")

library("tximeta")
library("GenomicAlignments")
library("DESeq2")
library("dplyr")
library("ggplot2")
library("pheatmap")
library("RColorBrewer")
library("PoiClaClu")
library("genefilter")
library("AnnotationDbi")
library("org.Hs.eg.db")
library("GO.db")
library("edgeR")
library("ComplexHeatmap")
library("EnhancedVolcano")
library("pathview")
library("gage")
library("gageData")
library("tibble")

dir <- "/PATH/TO/dir"
setwd(dir)
list.files(dir)
list.files(file.path(dir, "quants"))
csvfile <- file.path(dir, "Samples.csv")
coldata <- read.csv(csvfile, stringsAsFactors=FALSE)

Var = "group1vs2"
coldata <- coldata[c(1,3:7,10:21),] # Adjust for number of samples samples 2, 8 and 9 technical failures
coldata$names <- coldata$Sample
coldata$files <- file.path(dir, "quants", coldata$names, "quant.sf")
file.exists(coldata$files)

se <- tximeta(coldata)
dim(se)
head(rownames(se))
gse <- summarizeToGene(se)
dim(gse)
head(rownames(gse))

# Summarized experiment

assayNames(gse)
head(assay(gse), 3)
colSums(assay(gse))
rowRanges(gse)
seqinfo(rowRanges(gse))
colData(gse)

# We quickly check the millions of mapped reads round to number of decimal points
round(colSums(assay(gse))/1e6, 1) # 

dds <- DESeqDataSet(gse, design = ~Group) # for SALMON data

# Explore analysis and visualisation

dds$Group <- relevel(dds$Group, ref = "1") # to specify what the reference is, otherwise alphabetical according to condition
nrow(dds)

#Look at the variance

lambda <- 10^seq(from = -1, to = 2, length = 1000)
cts <- matrix(rpois(1000*100, lambda), ncol = 100)
library("vsn")
meanSdPlot(cts, ranks = FALSE) #Plots row standard deviations versus row means

# and for logarithm-transformed counts

log.cts.one <- log2(cts + 1)
meanSdPlot(log.cts.one, ranks = FALSE)  

vsd <- vst(dds, blind = FALSE)
head(assay(vsd), 3)
colData(vsd)
rld <- vst(dds, blind = FALSE)
head(assay(rld), 3)


dds <- estimateSizeFactors(dds)

df <- bind_rows(
  as_tibble(log2(counts(dds, normalized=TRUE)[, 1:2]+1)) %>%
    mutate(transformation = "log2(x + 1)"),
  as_tibble(assay(vsd)[, 1:2]) %>% mutate(transformation = "vst"),
  as_tibble(assay(rld)[, 1:2]) %>% mutate(transformation = "rlog"))

colnames(df)[1:2] <- c("x", "y")  

ggplot(df, aes(x = x, y = y)) + geom_hex(bins = 80) +
  coord_fixed() + facet_grid( . ~ transformation)  

#Sample distances

sampleDists <- dist(t(assay(vsd)))
sampleDists

sampleDistMatrix <- as.matrix( sampleDists )
rownames(sampleDistMatrix) <- paste(vsd$Name)
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette(c("grey85","burlywood4"))(8)
pheatmap(sampleDistMatrix,
         clustering_distance_rows = sampleDists,
         clustering_distance_cols = sampleDists,
         # breaks = c(0,25,50,75,100,125,150,175),
         col = colors)


poisd <- PoissonDistance(t(counts(dds)))
samplePoisDistMatrix <- as.matrix( poisd$dd )
rownames(samplePoisDistMatrix) <- paste( dds$Time, dds$Concentration, sep=" - " )
colnames(samplePoisDistMatrix) <- NULL
pheatmap(samplePoisDistMatrix,
         clustering_distance_rows = poisd$dd,
         clustering_distance_cols = poisd$dd,
         col = colors)

# PCA plot

plotPCA(vsd, intgroup = c("Concentration","Time")) + stat_ellipse()

# Differential expression analysis

dds <- DESeq(dds)

# Building the results table

res <- results(dds)
res

mcols(res, use.names = TRUE)
summary(res)

# Add gene names

res$ensembl <- sapply( strsplit( rownames(res), split="\\+" ), "[", 1 )
library( "biomaRt" )
ensembl = useMart( "ensembl", dataset = "hsapiens_gene_ensembl")
genemap <- getBM(attributes = c("ensembl_gene_id", "entrezgene_id", "hgnc_symbol"),
                 filters = "ensembl_gene_id",
                 values = res$ensembl,
                 mart = ensembl)
idx <- match( res$ensembl, genemap$ensembl_gene_id )
res$entrez <- genemap$entrezgene[ idx ]
res$hgnc_symbol <- genemap$hgnc_symbol[ idx ]
head(res,4)
write.csv(res, "name.csv")

# Add label mitocarta

csvfile <- file.path(dir, "MitoCartaGenes3.csv")
CartaData <- read.csv(csvfile, row.names=1, stringsAsFactors=FALSE)
idx2 <- match(res$ensembl, CartaData$EnsemblGeneID_mapping_version_20200130)
res$mito_localization <- CartaData$MitoCarta3.0_SubMitoLocalization[idx2]
res$mito_pathways <- CartaData$MitoCarta3.0_MitoPathways[idx2]
res_mitoOnly <- res[!is.na(res$mito_localization),]
head(res,4)
write.csv(res,"filename.csv")
write.csv(res_mitoOnly,"filename.csv")

# Heatmap most differentially expressed genes

res.df <- as.data.frame(res)
res.df <- res.df[(abs(res.df$padj)< 0.001),]
mat <- counts(dds, normalized = T)
write.csv(mat,"filename.csv")
mat.z <- t(apply(mat, 1, scale))
colnames(mat.z) <- coldata$Name
mat.z
heatmap(mat.z, distfun = dist, hclustfun = hclust, labRow = res$hgnc_symbol)

sigGenes <- subset(res, padj < 0.001)
sigGenes <- head(order(sigGenes$padj), increasing = TRUE,50)
mat <- assay(rld) [sigGenes,]

# Heatmap
mat <- mat - rowMeans(mat)
anno <- as.data.frame(colData(rld)[, c("Group")])
rownames(anno) <- colnames(mat)
pheatmap(mat, annotation_col = anno, labels_row = res$hgnc_symbol)

# Heatmap fold2change Mito Only

pheatmap(mat, annotation_col = anno, labels_row = res_mitoOnly$hgnc_symbol)

# Annotating and exporting results

columns(org.Hs.eg.db)

ens.str <- substr(rownames(res), 1, 15)
res$symbol <- mapIds(org.Hs.eg.db,
                     keys=ens.str,
                     column="SYMBOL",
                     keytype="ENSEMBL",
                     multiVals="first")
res$entrez <- mapIds(org.Hs.eg.db,
                     keys=ens.str,
                     column="ENTREZID",
                     keytype="ENSEMBL",
                     multiVals="first")

# GO Pathway analysis

# topVarGenes <- head(order(rowVars(assay(vsd)), decreasing = TRUE), 100)
SigGenes <- head(res[order(res$padj), ],increasing = TRUE, 50)
# number is significant from function above 
ListSigGenes <- SigGenes$entrez

go <- goana(ListSigGenes, species="Hs")
topGO<-topGO(go,n=50)
topGO
write.csv(topGO, "Filename.csv")

topBP<- topGO(go, ontology="BP", n=20)
topMF<- topGO(go, ontology="MF", n=20)
topCC<- topGO(go, ontology="CC", n=20)
topCC <- head(topCC[order(topCC$P.DE), ],increasing = TRUE, 15)
topMF <- head(topMF[order(topCC$P.DE), ],increasing = TRUE, 15)
topBP <- head(topBP[order(topCC$P.DE), ],increasing = TRUE, 15)

write.csv(topBP, "filename.csv")
write.csv(topMF, "filename.csv")
write.csv(topCC, "filename.csv")

# KEGG pathway analysis

keg <- kegga(ListSigGenes, species="Hs")
topKEGG<- topKEGG(keg, n=10, truncate=60)
topKEGG
write.csv(topKEGG, "filename.csv")

df<-topKEGG
df$Pathway <- factor(df$Pathway, levels = df$Pathway [order(rev(df$P.DE))]) #rev for increase
p <-ggplot(df, aes(x=Pathway, y= DE, fill=N)) +
  scale_fill_gradient(low="darkseagreen1", high="darkseagreen4") +
  geom_bar(stat = "identity")
p + coord_flip()

# Barchart Term versus Genenumber with different colors for BP, CC and MF

go <- goana(ListSigGenes, species="Hs")
topGO<-topGO(go,n=50)
topGO

df<-topGO
df$Term <- factor(df$Term, levels = df$Term [order(rev(df$P.DE))]) #rev for increase
p <-ggplot(df, aes(x=Term, y= DE, fill=Ont)) +
  geom_bar(stat = "identity")
p + coord_flip()

# Diverging barchart Term versus Genenumber with different colors for BP, CC and MF

df<-topMF
df$Term <- factor(df$Term, levels = df$Term [order(rev(df$P.DE))]) #rev for increase
p <-ggplot(df, aes(x=Term, y= DE, fill=Ont)) +
  geom_bar(stat = "identity")
p + coord_flip()

#GO:0031424 keratinization
#GO:0008544	epidermis development	
#GO:0030216	keratinocyte differentiation
#GO:0009913	epidermal cell differentiation	
#GO:0043588	skin development
#GO:0030855	epithelial cell differentiation
#GO:0070062	extracellular exosome
#GO:0022626 cytosolic ribosome
#GO:0005761 mitochondrial ribosome
#GO:0016266 O-glycan processing
#GO:0006954 inflammatory response
#GO:0033559	unsaturated fatty acid metabolic process
#GO:0001816	cytokine production
#GO:0005125	cytokine activity
#GO:0042593	glucose homeostasis
#GO:0007049 cell cycle

GoNumber <- "GO:1902686"

ensembl = useMart("ensembl",dataset="hsapiens_gene_ensembl")
gene.data <- getBM(attributes=c('hgnc_symbol', 'ensembl_gene_id', 'go_id'),
                   filters = 'go', values = GoNumber, mart = ensembl)
mart<- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
ensembl_genes <- gene.data$ensembl_gene_id
AllGenesGo <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","external_gene_name", "entrezgene_id","chromosome_name"),
  values= ensembl_genes, mart= mart)

# Make own list with genes
csvfile2 <- file.path(dir, "filename.csv")
coldata2 <- read.csv(csvfile2, stringsAsFactors=FALSE)
coldata2
ensembl_genes <- coldata2$Gene.stable.ID

# Heatmap with counts from a certain GOnumber

res.df <- as.data.frame(res)
AllGenesGo.df <- as.data.frame(AllGenesGo)
idx1 <- match(AllGenesGo$external_gene_name, res$hgnc_symbol)
AllGenesGo.df$`baseMean` <- res$baseMean [idx1]
AllGenesGo.df$`log2FoldChange` <- res$log2FoldChange [idx1]
AllGenesGo.df$`padj` <- res$padj [idx1]
AllGenesGo.df <- subset(AllGenesGo.df, baseMean !="") #remove missing genes (although this might be unnecessary for the next step)

# Optional, in case there are too many genes:
AllGenesGo.df <- subset(AllGenesGo.df, padj <0.01)

listGo <- AllGenesGo.df$ensembl_gene_id

mat <- counts(dds, normalized = T)
mat2<- subset(mat, rownames(mat) %in% listGo)
#  mat2<- subset(mat, rownames(mat) %in% res.df) dit werkt niet, want ik gebruik verder geen res.df
mat.z <- t(apply(mat2, 1, scale))
colnames(mat.z) <- coldata$Name
mat.z
h<- Heatmap(mat.z, cluster_rows = T, cluster_columns = T, column_labels = colnames(mat.z), name = "Z-score", row_labels = res.df[rownames(mat.z),]$symbol)
h

#Volcano plot with certain go number as label

 EnhancedVolcano(res,
    lab = res$hgnc_symbol,
    x = 'log2FoldChange',
    y = 'padj',
    selectLab = c(AllGenesGo.df$external_gene_name),
    xlim = c(-8,8),
    labSize = 4.0,
    legendLabels = c("NS", "Log2 FC", "Adjusted p-value", "Adjusted p-value & Log2 FC"),
    pCutoff = 0.01,
    FCcutoff = 0.5,
    drawConnectors = TRUE,
    widthConnectors = 0.2,
    maxoverlapsConnectors = Inf,
    colConnectors = 'grey30',
    title = "name")
 
 #Test pathway

foldchanges = res$log2FoldChange
names(foldchanges) = res$entrez
head(foldchanges)

# set up kegg database
kg.hsa=kegg.gsets(species="human")
kegg.sigmet.gs<-kg.hsa$kg.sets[kg.hsa$sigmet.idx]
kegg.dise.gs <- kg.hsa$kg.sets[kg.hsa$dise.idx]

# set up go database
data(go.sets.hs)
data(go.subs.hs)

go.hs=go.gsets(species="human") # If there is error pointer, restart R and empty global environment
go.bp.gs <- go.hs$go.sets[go.hs$go.subs$BP]
go.mf.gs <- go.hs$go.sets[go.hs$go.subs$MF]
go.cc.gs <- go.hs$go.sets[go.hs$go.subs$CC]

# Run enrichment analysis on all log fc
fc.kegg.sigmet.p <- gage(foldchanges, gsets = kegg.sigmet.gs)
fc.kegg.dise.p <- gage(foldchanges, gsets = kegg.dise.gs)
fc.go.bp.p <- gage(foldchanges, gsets = go.bp.gs)
fc.go.mf.p <- gage(foldchanges, gsets = go.mf.gs)
fc.go.cc.p <- gage(foldchanges, gsets = go.cc.gs)

# convert the kegg results to data frames
fc.kegg.sigmet.p.up <- as.data.frame(fc.kegg.sigmet.p$greater)
fc.kegg.dise.p.up <- as.data.frame(fc.kegg.dise.p$greater)

fc.kegg.sigmet.p.down <- as.data.frame(fc.kegg.sigmet.p$less)
fc.kegg.dise.p.down <- as.data.frame(fc.kegg.dise.p$less)

# convert the go results to data frames
fc.go.bp.p.up <- as.data.frame(fc.go.bp.p$greater)
fc.go.mf.p.up <- as.data.frame(fc.go.mf.p$greater)
fc.go.cc.p.up <- as.data.frame(fc.go.cc.p$greater)

fc.go.bp.p.down <- as.data.frame(fc.go.bp.p$less)
fc.go.mf.p.down <- as.data.frame(fc.go.mf.p$less)
fc.go.cc.p.down <- as.data.frame(fc.go.cc.p$less)

# Check most significant up
head(fc.kegg.sigmet.p.up, 10) 
head(fc.kegg.dise.p.up, 10) 

# View the most significant pathway from the pathway analysis
fc.kegg.sigmet.p.up[grepl("hsa00030", rownames(fc.kegg.sigmet.p.up), fixed=TRUE),]

# Overlay the expression data onto this pathway
pathview(gene.data=foldchanges, species="hsa", pathway.id="hsa00030")

# View as diagrams
fc.kegg.sigmet.p.up[grepl("hsa00030", rownames(fc.kegg.sigmet.p.up), fixed=TRUE),]
pathview(gene.data=foldchanges, species="hsa", pathway.id="hsa00030", kegg.native=FALSE)

# Check most significant down and up
head(fc.kegg.sigmet.p.down, 10)

# Similar with Gene Ontology

# Check most significant up
head(fc.go.mf.p.up, 10)
head(fc.go.cc.p.up, 10)
head(fc.go.bp.p.up, 10)

# Check most significant down
head(fc.go.mf.p.down, 10)
head(fc.go.cc.p.down, 10)
head(fc.go.bp.p.down, 10)

# Make a barchart from one of these:

df1<-head(fc.go.bp.p.down,11)
df1 <- df1 %>% rownames_to_column(var="Pathway")
df1$Pathway <- factor(df1$Pathway, levels = df1$Pathway [order(rev(df1$p.val))]) #rev for increase
p <-ggplot(df1, aes(x=Pathway, y= set.size, fill=set.size)) +
scale_fill_gradient(low="lightyellow2", high="darkgreen") +
geom_bar(stat = "identity")
p + coord_flip()
