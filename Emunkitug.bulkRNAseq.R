# Load libraries
library(DESeq2)
library(RSQLite)
library(data.table)
library(biomaRt)
library("factoextra")
library(fgsea)
library(ggsci)
library(ggpubr)
library("edgeR")
library("BiocParallel")
library(fgsea)
library(readxl)
library(plyr)
library(stringr)
library(dplyr)
library(ComplexHeatmap);library(circlize)


# Differential Expression ----

##### A) Data Munge ----

# Load & clean gene info from biomaRt 
human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
humBM <- getBM(attributes = c("ensembl_gene_id","hgnc_symbol","description"), mart = human)
humBM$description = gsub("\\s*\\[.*?\\]", "",humBM$description)

# Read meta-data
# meta = read.csv('~/HFB2003/Data/bulkRNAseq/2025/metadata.csv')
# meta$Sample = gsub('HFB200301_','',meta$Sample)
# 
# meta <- meta %>%
#   dplyr::mutate(
#     Patient = str_extract(Sample, "\\d{6}"),
#     Timepoint = case_when(
#       str_detect(Sample, "Screening|SCR") ~ "SCR",
#       str_detect(Sample, "C2D8") ~ "C2D8",
#       TRUE ~ NA_character_
#     )
#   )
# 
# rownames(meta) = meta$Sample
# 
# tracker = read.csv('~/HFB2003/Data/Tracker/Tracker_2.4.25.csv')
# tracker = tracker[,c('Subject','Tumor','Sched', 'Combo',   'Dose','BOR','Preceding.Tx.MoA','ToT.Month', 'TNFa', 'sTNFR2', 'TNFa.sTNFR2.hi', 'Prior.IO')]
# tracker$BOR = gsub(c('Death|Clinical PD'),'PD',tracker$BOR)
# tracker$BOR = gsub(c('SD>18wks'),'SD',tracker$BOR)
# tracker$BOR = gsub(c('PD->SD'),'SD',tracker$BOR)
# tracker$BOR = gsub(c('SD/PR'),'PR',tracker$BOR,fixed = T)
# tracker$Condition = ifelse(tracker$Subject %in%  unique(subset(tracker, BOR %in% c('PR','SD'))$Subject), 'Responder','Non-Responder')
# tracker$Patient = gsub('-','', tracker$Subject)
# 
# meta = merge(meta, tracker, by = 'Patient')

# Read in data matrix (read counts)
# dat = as.data.frame(fread('~/HFB2003/Data/bulkRNAseq/5.16/normalized_mtx.tsv'))
# dat = as.data.frame(fread('~/HFB2003/Data/bulkRNAseq/2025/normalized_mtx.tsv'))

meta = read.csv('~/HFB2003/Data/bulkRNAseq/2025/Metadata.2.26.25.csv')
dat = as.data.frame(fread('~/HFB2003/Data/bulkRNAseq/2025/raw_mtx.tsv'))

colnames(dat) = gsub('HFB200301_','',colnames(dat))


# specify gene names
gene_id = dat$V1
dat = dat[,-which(colnames(dat) == 'V1')]


##### Alexs File -----
df.log2 = log2(dat+1)
df.log2 = as.data.frame(df.log2)
df.log2$ensembl_gene_id =gene_id
df.log2 = merge(humBM, df.log2, by =  'ensembl_gene_id')
df.log2 = subset(df.log2, hgnc_symbol %in% c('CD8A','TNFRSF1B','KLRC1','PDCD1','CD274'))
df.log2 = melt(df.log2, id.var = c( 'ensembl_gene_id', 'hgnc_symbol','description'))
df.log2 = merge(meta, df.log2, by.x = 'Sample',by.y = 'variable')
df.df.log2 = df.log2[,c('Subject','Timepoint', 'hgnc_symbol','value')]
df.df.log2 = reshape(df.df.log2, idvar = c('Subject','Timepoint'), timevar = 'hgnc_symbol',direction = 'wide')
colnames(df.df.log2) = gsub('value.','',colnames(df.df.log2))
# write.csv(df.df.log2, '~/HFB2003/Data/bulkRNAseq/2025/Emunkitug.bulkRNAseq.log2.counts.2.21.25.csv',row.names = F)


# Subset only data at baseline
meta = subset(meta, Timepoint == 'SCR')
meta.scr = subset(meta, Timepoint == 'SCR')
dat = dat[,which(colnames(dat) %in% meta$Sample)]
meta = meta[match( colnames(dat), meta$Sample),]

# round data to nearest integer and make sure it's in numerical format
# dat = apply(dat,2, function(x) round(x))
# dat = apply(dat,2, function(x) as.numeric(as.character(x)))
# add gene ids 
rownames(dat) = gene_id
# put data into DESeq2 object
dds <- DESeqDataSetFromMatrix(dat, meta, design = ~ Condition)


##### PCA -----

# normalize data
normalized.dat = vst(dds, blind = FALSE)

pca.dat = as.data.frame(assay(normalized.dat))
pca.dat <- prcomp(t(pca.dat))
rownames(pca.dat$x) = meta$Patient 
fviz_pca_ind(pca.dat,  habillage= paste0(meta$Condition, meta$Combo), # this specifies the trait to draw a circle around 
             repel = T,
             addEllipses=TRUE, ellipse.level=0.5) + 
  theme(text = element_text(size = 23.5),
        axis.title = element_text(size = 17.5),
        axis.text = element_text(size = 17.5), legend.position = 'bottom') +
  guides(fill = guide_legend( override.aes = aes(label = ""))) 


##### TNFR2 & CD8 -----
scr = assay(dds)
scr = log2(scr+1)
scr = as.data.frame(scr)
scr$ensembl_gene_id = rownames(scr)
scr = merge(humBM, scr, by =  'ensembl_gene_id')
scr = subset(scr, hgnc_symbol %in% c('CD8A','TNFRSF1B'))
scr = melt(scr, id.var = c( 'ensembl_gene_id', 'hgnc_symbol','description'))
scr = merge(meta, scr, by.x = 'Sample',by.y = 'variable')
df.scr = scr[,c('Subject','Combo','BOR','hgnc_symbol','value')]
df.scr = reshape(df.scr, idvar = c('Subject','Combo','BOR'), timevar = 'hgnc_symbol',direction = 'wide')
colnames(df.scr)[4:5] = c('TNFR2','CD8A')

ggplot(df.scr, aes(x = TNFR2, y = CD8A, color = BOR, shape = Combo)) +
  geom_point(size = 6) +
  geom_hline(yintercept = 11.384786, linetype = 'dashed',color = 'red',linewidth = 1.2) +
  geom_vline(xintercept = 11.724590, linetype = 'dashed',color = 'red',linewidth = 1.2) +
  scale_color_lancet() +
  labs(title = 'Baseline', color = '', shape = '',subtitle ='bulk RNA-seq', y = 'CD8A (log2 counts)', x = 'TNFR2 (log2 counts)' ) +
  theme_minimal(base_size = 25) + theme(legend.position = 'top')


##### DEGS-----
dds <- DESeq(dds)
# specify comparisons. In this case comparing R vs NR
results <- results(dds, contrast = c('Condition', 'Responder', 'Non-Responder'))
##### D) Storing DESeq2 output ----
res = as.data.frame(results)
res$ensembl_gene_id = rownames(res)
res = merge(humBM,res, by = 'ensembl_gene_id')
# saveRDS(res, '~/HFB2003/Data/bulkRNAseq/2025/DEGs_PR.and.SD_vs_PD.all.scr.2.20.25.rds')

library(ggfx)


res$significant <- ifelse(res$padj < 0.05 & res$log2FoldChange > 1, 'S','NS')
res$significant <- ifelse(res$padj < 0.05 & res$log2FoldChange < -1, 'D',res$significant)

ggplot(res, aes(x = log2FoldChange, y = -log10(padj))) +
  with_shadow(
    geom_point(aes(color = significant), size = ifelse(res$significant %in% c('S','D'),2.5,.2), alpha = 0.8),
    sigma = 5, x_offset = 3, y_offset = 3
  ) +
  geom_hline(yintercept = -log10(0.05), color = 'red', linetype = 'dashed', linewidth = 1.2) +
  geom_vline(xintercept = c(-1, 1), color = 'black', linetype = 'dashed', linewidth = 1.2) +
  scale_color_manual(values = rev(c('red', 'grey85', 'cornflowerblue'))) +
  labs(x = expression(Log[2]~Fold~Change), y = expression(-Log[10]~P[adj]),
       title = 'PR & SD vs PD', subtitle = 'Baseline') +
  theme_classic(base_size = 25) +
  theme(legend.position = 'none') +
  ylim(0, 7.5) +
  xlim(-7.5, 7.5)







##### GSEA-----


hfb.genes = readRDS('~/GeneSets/HFB.Genesets.7.24.24.rds')

makeGeneRankStat <- function(statistic, names) {
  generank <- statistic
  names(generank) <- names
  generank <- generank[which(!is.na(generank))]
  generank <- generank[which(!is.nan(generank))]
  generank <- generank[which(names(generank) != "")]
  generank[which(generank == "-Inf")] <- min(generank[which(generank != "-Inf")], na.rm=T) -
    (0.01*min(generank[which(generank != "-Inf")], na.rm=T))
  generank[which(generank == "Inf")] <- max(generank[which(generank != "Inf")], na.rm=T) +
    (0.01*max(generank[which(generank != "Inf")], na.rm=T))
  generank <- generank[order(abs(generank), decreasing=T)]
  generank <- generank[!duplicated(names(generank))]
  generank <- generank[order(generank, decreasing=T)]
  return(generank) }


GSEA.enrich = function(x) {
  generank <- makeGeneRankStat(x$stat, toupper(x$hgnc_symbol))
  gsea.res <- fgseaMultilevel(pathways=hfb.genes, stats=generank, minSize=15,maxSize=500, eps = 0, nproc = 15)
  return(gsea.res)
}


res.gsea = GSEA.enrich(res)
# saveRDS(res.gsea, '~/HFB2003/Data/bulkRNAseq/6.4.24/GSEA.Responder.vs.NonResponder.DESeq2.rds')
res.gsea.sig = subset(res.gsea, padj < .05 & abs(NES) > 2)
res.gsea.sig$pathway = gsub('_',' ',res.gsea.sig$pathway)

# saveRDS(res.gsea, '~/HFB2003/Data/bulkRNAseq/2025/GSEA_PR.and.SD_vs_PD.all.scr.2.20.25.rds')

# res.gsea.sig = res.gsea.sig[-grep(c('ZHANG|PHONG|TIAN')),]

res.gsea.sig$pathway = gsub('HFB2003.PD','Emunkitug Tumor RNA-seq Response Signature v1',res.gsea.sig$pathway, fixed = T)

ggplot(res.gsea.sig, aes(x = -log10(padj), y = reorder(pathway, -padj))) +
  geom_bar(stat = 'identity', fill = 'navy', position = position_dodge2(.9)) +
  geom_vline(xintercept = -log10(.05), linetype = 'dashed', linewidth = 1.2, color = 'red') +
  # scale_fill_manual( values  = 'springgreen') +
  # scale_fill_gradient2(low = 'navy', mid = 'white',high = 'red',midpoint = 1) +
  theme_classic(base_size = 25) +
  labs(x = expression(-Log[10]~P[adj]), y = '', title = 'GSEA PR & SD vs PD', subtitle = 'Baseline') +
  theme(axis.text.y = element_text(size = 11))




##### Heatmap -----
base.degs = readRDS('~/HFB2003/Data/bulkRNAseq/2025/DEGs_PR.and.SD_vs_PD.all.scr.2.20.25.rds')
res.gsea = readRDS('~/HFB2003/Data/bulkRNAseq/2025/GSEA_PR.and.SD_vs_PD.all.scr.2.20.25.rds')

res.gsea.sig = subset(res.gsea, padj < .05 & abs(NES) >2)
gsea.genes = unlist(res.gsea.sig$leadingEdge)
gsea.overalp = intersect(baseline.genes, response.genes)

deg.genes = subset(base.degs, pvalue < .01 & abs(log2FoldChange) > 1)$hgnc_symbol
genes.of.interest = sort(intersect(deg.genes, gsea.genes))

Responders = subset(meta, BOR != 'PD')$Patient
  
  
dat = as.data.frame(fread('~/HFB2003/Data/bulkRNAseq/2025/normalized_mtx.tsv'))
colnames(dat) = gsub('HFB200301_','',colnames(dat))


meta.scr$Condition = ifelse(meta.scr$Patient %in% Responders, 'Responder','Non-Responder')

gene_id = dat$V1
dat = dat[,-which(colnames(dat) == 'V1')]
#colnames(dat) = substr(colnames(dat), 1, nchar(colnames(dat)) - 5)
#colnames(dat) = gsub("([^_]+)_\\d{2,3}M_(.*)", "\\1_\\2", colnames(dat))

dat = dat[,which(colnames(dat) %in% meta.scr$Sample)]
meta.scr = meta.scr[match( colnames(dat), meta.scr$Sample),]

rownames(dat) = gene_id
norm.mat = as.data.frame(dat)
norm.mat$ensembl_gene_id = rownames(norm.mat)
norm.mat = merge(humBM, norm.mat, by = 'ensembl_gene_id')
norm.mat = subset(norm.mat, hgnc_symbol %in% genes.of.interest)

rownames(norm.mat) = norm.mat$hgnc_symbol
norm.mat = norm.mat[,-which(colnames(norm.mat) %in% colnames(humBM))]
colnames(norm.mat) = substr(colnames(norm.mat),1,6)

pr =  meta.scr$Patient[which(meta.scr$Patient %in% meta.scr$Patient[meta.scr$BOR == 'PR'])]
sd =  meta.scr$Patient[which(meta.scr$Patient %in% meta.scr$Patient[meta.scr$BOR == 'SD'])]
pd =  meta.scr$Patient[which(meta.scr$Patient %in% meta.scr$Patient[meta.scr$BOR == 'PD'])]

norm.mat = norm.mat[,match(c(pd,sd,pr), colnames(norm.mat))]
norm.mat.scale.screen = t(apply(norm.mat, 1, function(x) scale(x)))
colnames(norm.mat.scale.screen) = colnames(norm.mat)


col_fun <- colorRamp2(c(-1.5, 0, 1.5), c("blue", "white", "red"))

# Cluster rows on the full matrix so that the order is the same across panels.
row_clust <- hclust(dist(norm.mat.scale.screen))
row_order <- rownames(norm.mat.scale.screen)[row_clust$order]

# Reorder the full matrix based on the clustering.
norm.mat.scale.screen <- norm.mat.scale.screen[row_order, ]

# Split the matrix into three parts:
mat1 <- norm.mat.scale.screen[, 1:22]
mat2 <- norm.mat.scale.screen[, 23:28]
mat3 <- norm.mat.scale.screen[, 29:31]

# Define the genes for which row labels should be shown.
selected_genes <- c('CCL2', 'CCL8', 'CXCL5', 'CXCL8', 'CXCL9',
                    'IFIT1', 'IFIT2', 'IFITM1', 'IFITM2',
                    'IL6R', 'IL1RAP', 'CSF1R', 'IL6', 'IL1RN',
                    'NFKBIZ','CD44', 'CD59','CD163', 'CD14', 'CD58',
                    'CD82', "TNFRSF11B")

# Create custom row labels: only show the gene name if it is in selected_genes.
row_labels <- ifelse(rownames(norm.mat.scale.screen) %in% selected_genes,
                     rownames(norm.mat.scale.screen), "")

# Create the three heatmap objects.
# For mat1 and mat2, we don't show row names (to avoid repetition).
ht1 <- Heatmap(mat1, 
               cluster_rows = FALSE,  # already clustered
               cluster_columns = FALSE,
               show_row_names = FALSE,
               col = col_fun,
               column_names_rot = 75,
               column_names_gp = gpar(fontsize = 18), show_heatmap_legend = F)

ht2 <- Heatmap(mat2,
               cluster_rows = FALSE,
               cluster_columns = FALSE,
               show_row_names = FALSE,
               col = col_fun,
               column_names_rot = 75,
               column_names_gp = gpar(fontsize = 18), show_heatmap_legend = F)

# For the third panel, we display row labels.
ht3 <- Heatmap(mat3,
               cluster_rows = FALSE,
               cluster_columns = FALSE,
               show_row_names = TRUE,
               row_labels = row_labels,
               col = col_fun,
               column_names_side = "bottom",  # place column names at the bottom
               column_names_rot = 75,
               column_names_gp = gpar(fontsize = 18),
               row_names_gp = gpar(fontsize = 10),
               heatmap_legend_param = list(
                 title = "Z(VST)",
                 at = c(-2, -1, 0, 1, 2),
                 labels = c("-2", "-1", "0", "1", "2")
               ),
               name = "Panel3")

# Combine the three heatmaps horizontally.
ht_list <- ht1 + ht2 + ht3

# Draw the combined heatmaps.
draw(ht_list, heatmap_legend_side = "right" ,gap = unit(10, "mm"))





##### TCGA Enrich -----
base.degs = readRDS('~/HFB2003/Data/bulkRNAseq/2025/DEGs_PR.and.SD_vs_PD.all.scr.2.20.25.rds')
res.gsea =  readRDS('~/HFB2003/Data/bulkRNAseq/2025/GSEA_PR.and.SD_vs_PD.all.scr.2.20.25.rds')

res.gsea.sig = subset(res.gsea, padj < .05 & abs(NES) >2)
gsea.genes = unlist(res.gsea.sig$leadingEdge)

deg.genes = subset(base.degs, pvalue < .01 & abs(log2FoldChange) > 1)$hgnc_symbol
genes.of.interest = sort(intersect(deg.genes, gsea.genes))

hfb.genes = list(Emunkitug = genes.of.interest )


GSEA.enrich = function(x) {
  generank <- makeGeneRankStat(x$Mean.z, toupper(x$hgnc_symbol))
  gsea.res <- fgseaMultilevel(pathways=hfb.genes, stats=generank, minSize=15,maxSize=500, eps = 0, nproc = 15)
  return(gsea.res)
}


makeGeneRankStat <- function(statistic, names) {
  generank <- statistic
  names(generank) <- names
  generank <- generank[which(!is.na(generank))]
  generank <- generank[which(!is.nan(generank))]
  generank <- generank[which(names(generank) != "")]
  generank[which(generank == "-Inf")] <- min(generank[which(generank != "-Inf")], na.rm=T) -
    (0.01*min(generank[which(generank != "-Inf")], na.rm=T))
  generank[which(generank == "Inf")] <- max(generank[which(generank != "Inf")], na.rm=T) +
    (0.01*max(generank[which(generank != "Inf")], na.rm=T))
  generank <- generank[order(abs(generank), decreasing=T)]
  generank <- generank[!duplicated(names(generank))]
  generank <- generank[order(generank, decreasing=T)]
  return(generank) }

# 
# 
# tcga.z = readRDS('~/TCGA/Data/TCGA.Subtype.Zstats.rds')
# tcga.gsea.res <- plyr::ddply(tcga.z, .(Cancer), GSEA.enrich)
# saveRDS(tcga.gsea.res, '~/HFB2003/Data/bulkRNAseq/2025/TCGA.Emunkitug.Signature.Enrich.All.rds')
# 
# 
# tcga.z = readRDS('~/TCGA/Data/TCGA.Subtype.Zstats.PDL1.rds')
# tcga.gsea.res <- plyr::ddply(tcga.z, .(variable), GSEA.enrich)
# saveRDS(tcga.gsea.res, '~/HFB2003/Data/bulkRNAseq/2025/TCGA.Emunkitug.Signature.Enrich.PDL1.rds')
# 



tcga.all = readRDS('~/HFB2003/Data/bulkRNAseq/2025/TCGA.Emunkitug.Signature.Enrich.All.rds') 
tcga.pd1 = readRDS('~/HFB2003/Data/bulkRNAseq/2025/TCGA.Emunkitug.Signature.Enrich.PDL1.rds') 

tcga.pd1 = tcga.pd1[grep('PDL1',tcga.pd1$variable),]
tcga.pd1 = tcga.pd1[grep(c('MESO|HNSC|CESC|LUAD'),tcga.pd1$variable),]
colnames(tcga.pd1)[1] = 'Cancer'
tcga.pd1 = subset(tcga.pd1,  pathway == 'Emunkitug')

tcga.all = as.data.frame(rbind(tcga.all, tcga.pd1))
tcga.all = tcga.all[order(tcga.all$padj),]
tcga.all$rank = 1:64

tcga.pre = read.csv('~/HFB2003/Data/Indication.Rank.Preclinical.csv')

#tcga.pre$perc = (1 - (tcga.pre$Rank - 1) / max(tcga.pre$Rank))*100
#tcga.all$perc = (1 - (tcga.all$rank - 1) / max(tcga.all$rank))*100

tcga.all <- tcga.all %>%
  mutate(Cancer.broad = sub("_.*", "", Cancer))

tcga.both = merge(tcga.pre, tcga.all, by.x = 'Cancer', by.y = 'Cancer.broad')

pdl1.drop1 = intersect( grep('PDL1_high',tcga.both$Cancer.specific), grep('PDL1_low',tcga.both$Cancer.y))
pdl1.drop2 = intersect( grep('PDL1_low',tcga.both$Cancer.specific), grep('PDL1_high',tcga.both$Cancer.y))
stad.drop =  intersect(grep('STAD_EBV+',tcga.both$Cancer.specific, fixed = T),  grep(c('STAD_MSS|STAD_MSI'),tcga.both$Cancer.y))
stad.drop2 = which(tcga.both$Cancer.specific == 'STAD_EBV-' & tcga.both$Cancer.y == 'STAD_EBV+')

tcga.plot = tcga.both[-c(pdl1.drop1,pdl1.drop2,stad.drop,stad.drop2),]

tcga.plot$perc.x = (1 - (tcga.plot$Rank - 1) / max(tcga.plot$Rank))*100
tcga.plot$perc.y = (1 - (tcga.plot$rank - 1) / max(tcga.plot$rank))*100

library(ggrepel)
ggplot(tcga.plot, aes(x = perc.x, y =  perc.y)) +
  geom_point(aes(colour = Cancer), size = 2, alpha = .75)+
  geom_text_repel(aes(label = Cancer.y), size = 4, box.padding = unit(0.35, "lines"), point.padding = unit(0.5, "lines"), max.overlaps = Inf) +
  geom_hline(yintercept =  75, linetype = 'dashed', color = 'red', linewidth = 1.5) + geom_vline(xintercept = 75, color = 'red',linetype = 'dashed', linewidth = 1.5) +
  labs(title = '', color = '', x = 'Indication Rank Preclinical (%tile)' , y = 'Indication Rank Phase 1 (%tile)') +
  theme_classic(base_size = 30) + 
  theme( legend.position = 'none', strip.text = element_text(size=25),strip.background = element_blank())

tcga.plot2 = tcga.plot[-which(duplicated(tcga.plot$Cancer.y)),]
tcga.plot2$signif =  ifelse(tcga.plot2$padj < .01, 'Sig','NS')

ggplot(tcga.plot2, aes(x = reorder(Cancer.y, -padj), y = -log10(padj), fill = signif)) +
  geom_bar(stat = 'identity', position = position_dodge(.9)) +
  theme_classic(base_size = 21.5) + labs(x = '', y = expression(-Log[10]~P[adj])) + theme(legend.position = 'none') +
  scale_fill_manual(values = c(Sig = "springgreen", NS= "navyblue" )) +  # Set gradient colors from blue to red
  #scale_fill_gradient(low = "blue", high = "red" ) +  # Set gradient colors from blue to red
  geom_hline(yintercept = -log10(.01), linetype = 'dashed', color = 'black', linewidth = 1.5) +
  theme(aes(axis.text.x = element_text(colour = lab.red))) +
  coord_flip()




#### Indication ranking

tcga.pd1 = readRDS('~/HFB2003/Data/bulkRNAseq/6.4.24/TCGA_GSEA_PDL1.HFB2003.pval.025.Non.Rescue.rds') 
tcga.all = readRDS('~/HFB2003/Data/bulkRNAseq/6.4.24/TCGA_GSEA_HFB2003.pval.025.rds') 

tcga.pd1 = tcga.pd1[grep('PDL1',tcga.pd1$variable),]
tcga.pd1 = tcga.pd1[grep(c('MESO|HNSC|CESC|LUAD'),tcga.pd1$variable),]
colnames(tcga.pd1)[1] = 'Cancer'
tcga.pd1 = subset(tcga.pd1,  pathway == 'HFB2003_Tumor_RNAseq_PD_Biomarkers')

tcga.all = as.data.frame(rbind(tcga.all, tcga.pd1))
tcga.all = tcga.all[order(tcga.all$padj),]
tcga.all$rank = 1:64

tcga.pre = read.csv('~/HFB2003/Data/Indication.Rank.Preclinical.csv')

#tcga.pre$perc = (1 - (tcga.pre$Rank - 1) / max(tcga.pre$Rank))*100
#tcga.all$perc = (1 - (tcga.all$rank - 1) / max(tcga.all$rank))*100

tcga.all <- tcga.all %>%
  mutate(Cancer.broad = sub("_.*", "", Cancer))

tcga.both = merge(tcga.pre, tcga.all, by.x = 'Cancer', by.y = 'Cancer.broad')

pdl1.drop1 = intersect( grep('PDL1_high',tcga.both$Cancer.specific), grep('PDL1_low',tcga.both$Cancer.y))
pdl1.drop2 = intersect( grep('PDL1_low',tcga.both$Cancer.specific), grep('PDL1_high',tcga.both$Cancer.y))
stad.drop =  intersect(grep('STAD_EBV+',tcga.both$Cancer.specific, fixed = T),  grep(c('STAD_MSS|STAD_MSI'),tcga.both$Cancer.y))
stad.drop2 = which(tcga.both$Cancer.specific == 'STAD_EBV-' & tcga.both$Cancer.y == 'STAD_EBV+')

tcga.plot = tcga.both[-c(pdl1.drop1,pdl1.drop2,stad.drop,stad.drop2),]

tcga.plot$perc.x = (1 - (tcga.plot$Rank - 1) / max(tcga.plot$Rank))*100
tcga.plot$perc.y = (1 - (tcga.plot$rank - 1) / max(tcga.plot$rank))*100

library(ggrepel)
ggplot(tcga.plot, aes(x = perc.x, y =  perc.y)) +
  geom_point(aes(colour = Cancer), size = 2, alpha = .75)+
  geom_text_repel(aes(label = Cancer.y), size = 4, box.padding = unit(0.35, "lines"), point.padding = unit(0.5, "lines"), max.overlaps = Inf) +
  geom_hline(yintercept =  75, linetype = 'dashed', color = 'red', linewidth = 1.5) + geom_vline(xintercept = 75, color = 'red',linetype = 'dashed', linewidth = 1.5) +
  labs(title = '', color = '', x = 'Indication Rank Preclinical (%tile)' , y = 'Indication Rank Ph1 A (%tile)') +
  theme_classic(base_size = 30) + 
  theme( legend.position = 'none', strip.text = element_text(size=25),strip.background = element_blank())

tcga.plot2 = tcga.plot[-which(duplicated(tcga.plot$Cancer.y)),]
tcga.plot2$signif =  ifelse(tcga.plot2$padj < .01, 'Sig','NS')

ggplot(tcga.plot2, aes(x = reorder(Cancer.y, -padj), y = -log10(padj), fill = signif)) +
  geom_bar(stat = 'identity', position = position_dodge(.9)) +
  theme_classic(base_size = 21.5) + labs(x = '', y = expression(-Log[10]~P[adj])) + theme(legend.position = 'none') +
  scale_fill_manual(values = c(Sig = "springgreen", NS= "navyblue" )) +  # Set gradient colors from blue to red
  #scale_fill_gradient(low = "blue", high = "red" ) +  # Set gradient colors from blue to red
  geom_hline(yintercept = -log10(.01), linetype = 'dashed', color = 'black', linewidth = 1.5) +
  theme(aes(axis.text.x = element_text(colour = lab.red))) +
  coord_flip()






# PD EFFECTS -----
library("variancePartition")
library("edgeR")
library("BiocParallel")

meta = read.csv('~/HFB2003/Data/bulkRNAseq/2025/metadata.csv')
meta$Sample = gsub('HFB200301_','',meta$Sample)

meta <- meta %>%
  dplyr::mutate(
    Patient = str_extract(Sample, "\\d{6}"),
    Timepoint = case_when(
      str_detect(Sample, "Screening|SCR") ~ "SCR",
      str_detect(Sample, "C2D8") ~ "C2D8",
      TRUE ~ NA_character_
    )
  )

rownames(meta) = meta$Sample

tracker = read.csv('~/HFB2003/Data/Tracker/Tracker_2.4.25.csv')
tracker = tracker[,c('Subject','Tumor','Sched', 'Combo',   'Dose','BOR','Preceding.Tx.MoA','ToT.Month', 'TNFa', 'sTNFR2', 'TNFa.sTNFR2.hi', 'Prior.IO')]
tracker$BOR = gsub(c('Death|Clinical PD'),'PD',tracker$BOR)
tracker$BOR = gsub(c('SD>18wks'),'SD',tracker$BOR)
tracker$BOR = gsub(c('PD->SD'),'SD',tracker$BOR)
tracker$BOR = gsub(c('SD/PR'),'PR',tracker$BOR,fixed = T)
tracker$Condition = ifelse(tracker$Subject %in%  unique(subset(tracker, BOR %in% c('PR','SD'))$Subject), 'Responder','Non-Responder')
tracker$Patient = gsub('-','', tracker$Subject)

meta = merge(meta, tracker, by = 'Patient')

dat = as.data.frame(fread('~/HFB2003/Data/bulkRNAseq/2025/raw_mtx.tsv'))
colnames(dat) = gsub('HFB200301_','',colnames(dat))
gene_id = dat$V1
dat = dat[,-which(colnames(dat) == 'V1')]

paired = subset(as.data.frame(table(meta$Subject)), Freq>1)$Var

# Subset only data at baseline
# meta.pair = subset(meta, Subject %in%  paired)

meta.pair = subset(meta, Subject %in%  paired & BOR != 'PD')

# meta.pair = subset(meta, Subject %in%  paired & Combo == 'Combo')
# meta.pair = subset(meta, Subject %in%  paired & Combo == 'Mono')

dat = dat[,which(colnames(dat) %in% meta.pair$Sample)]
meta.pair = meta.pair[match( colnames(dat), meta.pair$Sample),]

rownames(dat) = gene_id
dds <- DESeqDataSetFromMatrix(dat, meta.pair, design = ~ Timepoint)

# normalized.dat = vst(dds, blind = FALSE)
# 
# pca.dat = as.data.frame(assay(normalized.dat))
# pca.dat <- prcomp(t(pca.dat))
# rownames(pca.dat$x) = meta$Patient 
# fviz_pca_ind(pca.dat,  habillage= paste0(meta.pair$Condition, meta.pair$Timepoint), # this specifies the trait to draw a circle around 
#              repel = T,
#              addEllipses=TRUE, ellipse.level=0.5) + 
#   theme(text = element_text(size = 23.5),
#         axis.title = element_text(size = 17.5),
#         axis.text = element_text(size = 17.5), legend.position = 'bottom') +
#   guides(fill = guide_legend( override.aes = aes(label = ""))) 


isexpr <- rowSums(cpm(assay(dds)) > 1) >= 4

rownames(meta.pair) = meta.pair$Sample

# Standard usage of limma/voom
dge <- dds[isexpr, ]
dge <- calcNormFactors(dge)
all(colnames(dge) %in% rownames(meta.pair))
design = model.matrix(~Timepoint, meta.pair)
# param <- SnowParam(4, "SOCK", progressbar = TRUE)
# clusterEvalQ(param$cluster, { library(stats) })
form <- ~ Timepoint + (1 | Subject)
vobjDream <- voomWithDreamWeights(dge, form, meta.pair)
fitmm <- dream(vobjDream, form, meta.pair)
fitmm <- eBayes(fitmm)


head(fitmm$design, 3)
res.pre.post = topTable(fitmm, coef = "TimepointSCR", number = dim(fitmm)[1])
res.pre.post$ensembl_gene_id = rownames(res.pre.post)
res.pre.post = merge(humBM, res.pre.post, by = 'ensembl_gene_id')


resp.all = res.pre.post
resp.all$Condition = 'Emunkitug.Post.vs.Pre.Responders'

mono.all = res.pre.post
mono.all$Condition = 'Emunkitug.Post.vs.Pre.Mono'

combo.all = res.pre.post
combo.all$Condition = 'Emunkitug.Post.vs.Pre.Combo'

pd.all = res.pre.post
pd.all$Condition = 'Emunkitug.Post.vs.Pre.All'

pd.rna = as.data.frame(rbind(mono.all, combo.all,pd.all,resp.all))

# saveRDS(pd.rna, '~/HFB2003/Data/bulkRNAseq/2025/DEG.DREAM.PD.Mono.Combo.All.Responder.rds')

pd.rna = readRDS('~/HFB2003/Data/bulkRNAseq/2025/DEG.DREAM.PD.Mono.Combo.All.Responder.rds')
pd.rna$logFC = pd.rna$logFC*-1
pd.rna$z.std = pd.rna$z.std*-1
pd.rna$t = pd.rna$t*-1

pd.rna$Condition = gsub('Emunkitug.Post.vs.Pre.','',pd.rna$Condition, fixed = T)
pd.rna$Condition = gsub('Responders','PR & SD',pd.rna$Condition, fixed = T)
pd.rna$se = abs(  (pd.rna$logFC / pd.rna$z.std)   )

ggplot(subset(pd.rna, hgnc_symbol %in% c('TNFRSF1B','CD8A','PDCD1','CD274') ) , aes(x = hgnc_symbol, y = logFC, fill = hgnc_symbol)) +
  geom_bar(position = position_dodge(.9), stat = "identity") + ylim(-1.25,2) +
  geom_errorbar(aes(ymin = logFC-se, ymax = logFC+se),position = position_dodge(.9), width = 0.3, color = "black") +
  geom_hline(yintercept = c(0), color = 'black', linewidth = 1) +
  scale_fill_simpsons() +
  geom_text(position = position_dodge(.9), color = 'black', aes(  y = logFC +sign(logFC) * se * 1.25 ,  label = ifelse(P.Value < 0.0001, '****', ifelse(P.Value < .001, '***', ifelse( P.Value < .01, '**',       ifelse(P.Value < .05, '*', ifelse(P.Value < .1 & P.Value > .05, '', '')))))) , size = 7.5) +
  facet_wrap(~Condition,scales = 'free',nrow = 1) +
  labs(x = "", y = expression(Log[2]~FC), fill = "",title = 'Emunkitug Post- vs Pre-treatment',subtitle = 'Tumor bulkRNA-seq') +
  theme_minimal(base_size = 30) + rotate_x_text(angle = 65) +
  theme( legend.position = 'none') 


###### 2) PD Boxplot----
meta = read.csv('~/HFB2003/Data/bulkRNAseq/2025/Metadata.2.26.25.csv')
 dat = as.data.frame(fread('~/HFB2003/Data/bulkRNAseq/2025/raw_mtx.tsv'))
 colnames(dat) = gsub('HFB200301_','',colnames(dat))
 gene_id = dat$V1
 dat = dat[,-which(colnames(dat) == 'V1')]
 
 paired = subset(as.data.frame(table(meta$Subject)), Freq>1)$Var
 
 # Subset only data at baseline
 meta.pair = subset(meta, Subject %in%  paired)
 # 
 # meta.pair = subset(meta, Subject %in%  paired & BOR != 'PD')
 # # meta.pair = subset(meta, Subject %in%  paired & Combo == 'Combo')
 # meta.pair = subset(meta, Subject %in%  paired & Combo == 'Mono')
 
 dat = dat[,which(colnames(dat) %in% meta.pair$Sample)]
 meta.pair = meta.pair[match( colnames(dat), meta.pair$Sample),]
 
 rownames(dat) = gene_id
 dds <- DESeqDataSetFromMatrix(dat, meta.pair, design = ~ Timepoint)
 
 # isexpr <- rowSums(cpm(assay(dds)) > 1) >= 0
 
 rownames(meta.pair) = meta.pair$Sample
 
 # Standard usage of limma/voom
 # dge <- dds[isexpr, ]
 dge <- calcNormFactors(dds)
 
 df = as.data.frame(dge$counts)
 df$gene = rownames(df)
 df = merge(humBM, df, by.x = 'ensembl_gene_id',by.y = 'gene')
 df = melt(df, id.vars  = c('ensembl_gene_id', 'hgnc_symbol', 'description'))
 df = merge( meta, df,  by.x = 'Sample',by.y = 'variable')
 df$log2.counts = log2(df$value+1)
 
 df$Timepoint = factor(df$Timepoint, levels = c('SCR','C2D8'))
 ggplot(subset(df, hgnc_symbol %in% c('TNFRSF1B','PDCD1','CD8A','CD274' )), aes(x = Timepoint , y = log2.counts, color = Timepoint)) +
   geom_boxplot() +  # Add boxplots
   geom_line(aes(group = Patient), color = 'black', linetype = 'dashed') +  # Add lines connecting pre and post
   geom_point(aes(color = Timepoint), size = 3, position = position_dodge(width = 0.3)) +
   labs(x = "", y = "Normalized Counts (Log2+1)", title = 'Emunkitug PD Effects') + scale_color_jama() +
   theme_classic(base_size = 28) + 
   theme(legend.position = 'none', strip.background = element_blank()) +  rotate_x_text(angle = 65) + 
   facet_grid(Combo~hgnc_symbol, scales = 'free')

 ggplot(subset(df, hgnc_symbol %in% c('TNFRSF1B','PDCD1','CD8A','CD274' )), aes(x = Timepoint , y = log2.counts, color = Timepoint)) +
   geom_boxplot() +  # Add boxplots
   geom_line(aes(group = Patient), color = 'black', linetype = 'dashed') +  # Add lines connecting pre and post
   geom_point(aes(color = Timepoint), size = 3, position = position_dodge(width = 0.3)) +
   labs(x = "", y = "Normalized Counts (Log2+1)", title = 'Emunkitug PD Effects') + scale_color_jama() +
   theme_classic(base_size = 28) + 
   theme(legend.position = 'none', strip.background = element_blank()) +  rotate_x_text(angle = 65) + 
   facet_grid(BOR~hgnc_symbol, scales = 'free')
 
 
 
 ggplot(subset(df, Dose == '50 mg'& hgnc_symbol %in% c('TNFRSF1B','PDCD1','CD8A','CD274' )), aes(x = Timepoint , y = log2.counts, color = Timepoint)) +
   geom_boxplot() +  # Add boxplots
   geom_line(aes(group = Patient), color = 'black', linetype = 'dashed') +  # Add lines connecting pre and post
   geom_point(aes(color = Timepoint), size = 3, position = position_dodge(width = 0.3)) +
   labs(x = "", y = "Normalized Counts (Log2+1)", title = 'Emunkitug 50mg PD Effects') + scale_color_jama() +
   theme_classic(base_size = 28) + 
   theme(legend.position = 'none', strip.background = element_blank()) +  rotate_x_text(angle = 65) + 
   facet_grid(paste0(Combo,Sched)~hgnc_symbol, scales = 'free')


hfb.genes = readRDS('~/GeneSets/HFB.Genesets.2.24.25.rds')


makeGeneRankStat <- function(statistic, names) {
  generank <- statistic
  names(generank) <- names
  generank <- generank[which(!is.na(generank))]
  generank <- generank[which(!is.nan(generank))]
  generank <- generank[which(names(generank) != "")]
  generank[which(generank == "-Inf")] <- min(generank[which(generank != "-Inf")], na.rm=T) -
    (0.01*min(generank[which(generank != "-Inf")], na.rm=T))
  generank[which(generank == "Inf")] <- max(generank[which(generank != "Inf")], na.rm=T) +
    (0.01*max(generank[which(generank != "Inf")], na.rm=T))
  generank <- generank[order(abs(generank), decreasing=T)]
  generank <- generank[!duplicated(names(generank))]
  generank <- generank[order(generank, decreasing=T)]
  return(generank) }


GSEA.enrich = function(x) {
  generank <- makeGeneRankStat(x$t, toupper(x$hgnc_symbol))
  gsea.res <- fgseaMultilevel(pathways=hfb.genes, stats=generank, minSize=15,maxSize=500, eps = 0, nproc = 15)
  return(gsea.res)
}

pd.rna = readRDS('~/HFB2003/Data/bulkRNAseq/2025/DEG.DREAM.PD.Mono.Combo.All.Responder.rds')
pd.rna$logFC = pd.rna$logFC*-1
pd.rna$z.std = pd.rna$z.std*-1
pd.rna$t = pd.rna$t*-1

pd.rna$Condition = gsub('Emunkitug.Post.vs.Pre.','',pd.rna$Condition, fixed = T)
pd.rna$Condition = gsub('Responders','PR & SD',pd.rna$Condition, fixed = T)
pd.rna$se = abs(  (pd.rna$logFC / pd.rna$z.std)   )

pd.gsea.res <- plyr::ddply(pd.rna, .(Condition), GSEA.enrich)

pathways = c('Cytotoxicity','Cytokine_Signaling', 'Chemokine_Signaling','Interferon_Signaling','Lymphocyte_Activation','Emunkitug','NF_kappaB_Signaling','T_Cell_Functions','NK_Cell_Functions','TNF_Family_Signaling','Adaptive_Immune_System','Innate_Immune_System')

hfb.gsea = subset(pd.gsea.res, pathway %in% pathways)

hfb.gsea$pval_label <- case_when(
  hfb.gsea$pval < 0.0001 ~ "****",
  hfb.gsea$pval < 0.001 ~ "***",
  hfb.gsea$pval < 0.01 ~ "**",
  hfb.gsea$pval < 0.05 ~ "*",
  TRUE ~ ""
)

color_fun <- scales::col_numeric(
  palette = c("blue", "white", "red"),
  domain = c(min(hfb.gsea$NES), 0, max(hfb.gsea$NES)))

ggplot(hfb.gsea , aes(x = Condition, y = pathway)) +
  geom_tile(aes(fill = NES), color = "white") +
  scale_fill_gradientn(colours = c("blue", "white", "red"),
                       breaks = c(min(hfb.gsea$NES), 0, max(hfb.gsea$NES)),
                       labels = scales::label_number(),
                       limits = c(min(hfb.gsea$NES), max(hfb.gsea$NES))) +
  geom_text(aes(label = pval_label), color = "black", size = 5) +
  theme_minimal(base_size = 20) +
  labs(x = "", y = "", title = 'GSEA: Emunkitug Post- vs Pre-Treatment') +
  rotate_x_text( angle = 45)











# 
# OX40 ----
base.degs = readRDS('~/HFB2003/Data/bulkRNAseq/2025/DEGs_PR.and.SD_vs_PD.all.scr.2.20.25.rds')
res.gsea =  readRDS('~/HFB2003/Data/bulkRNAseq/2025/GSEA_PR.and.SD_vs_PD.all.scr.2.20.25.rds')

res.gsea.sig = subset(res.gsea, padj < .05 & abs(NES) >2)
gsea.genes = unlist(res.gsea.sig$leadingEdge)

deg.genes = subset(base.degs, pvalue < .01 & abs(log2FoldChange) > 1)$hgnc_symbol
genes.of.interest = sort(intersect(deg.genes, gsea.genes))

Emunkitug = list(Emunkitug = genes.of.interest )


meta = read.csv('~/HFB3010/Data/bulkRNAseq/metadata_final.csv')

trt = read_xlsx('~/HFB3010/Data/Clinical/OX40 Clinical Data_16AUG2024.xlsx')
trt = trt[,c('Subject', 'Preceding.Therapy')]
trt$Preceding.Therapy[which(trt$Preceding.Therapy %in% c('C','T','P+T'))] = 'NonIO'
trt$Preceding.Therapy[which(trt$Preceding.Therapy %in% c('P'))] = 'IO'
colnames(trt) = c('Subject','PriorTx')

# meta = merge(meta, trt, by = 'Subject')
meta = merge(meta, trt, by = 'Subject')

# What is normalized ?? Batch normalized? Should I round to nearest integer?
dat = as.data.frame(fread('~/HFB3010/Data/bulkRNAseq/normalized_mtx.tsv'))
# Remove outlier
# meta = subset(meta, Patient != '401008')

gene_id = dat$gene_id
dat = dat[,-which(colnames(dat) == 'gene_id')]

meta = subset(meta, TimePoint == 'Screening' )
# meta = subset(meta, TimePoint == 'Screening' & Batch == 'Batch3' & !(Subject %in% c('301-004','303-001')))
dat = dat[,which(colnames(dat) %in% meta$Sample)]
meta = meta[match( colnames(dat), meta$Sample),]

dat = apply(dat,2, function(x) round(x))
dat = apply(dat,2, function(x) as.numeric(as.character(x)))
rownames(dat) = gene_id

dds <- DESeqDataSetFromMatrix(dat, meta, design = ~ PriorTx)


# normalize data
normalized.dat = vst(dds, blind = FALSE)

pca.dat = as.data.frame(assay(normalized.dat))
pca.dat <- prcomp(t(pca.dat))
rownames(pca.dat$x) = meta$Patient 
fviz_pca_ind(pca.dat,  habillage= paste0(meta$Batch), # this specifies the trait to draw a circle around 
             repel = T,
             addEllipses=TRUE, ellipse.level=0.5) + 
  theme(text = element_text(size = 23.5),
        axis.title = element_text(size = 17.5),
        axis.text = element_text(size = 17.5), legend.position = 'bottom') +
  guides(fill = guide_legend( override.aes = aes(label = ""))) 






dds <- DESeq(dds)
# specify comparisons. In this case comparing R vs NR
results <- results(dds, contrast = c('PriorTx', 'IO', 'NonIO'))
res = as.data.frame(results)
res$ensembl_gene_id = rownames(res)
res = merge(humBM,res, by = 'ensembl_gene_id')


hfb.genes = readRDS('~/GeneSets/HFB.Genesets.7.24.24.rds')

makeGeneRankStat <- function(statistic, names) {
  generank <- statistic
  names(generank) <- names
  generank <- generank[which(!is.na(generank))]
  generank <- generank[which(!is.nan(generank))]
  generank <- generank[which(names(generank) != "")]
  generank[which(generank == "-Inf")] <- min(generank[which(generank != "-Inf")], na.rm=T) -
    (0.01*min(generank[which(generank != "-Inf")], na.rm=T))
  generank[which(generank == "Inf")] <- max(generank[which(generank != "Inf")], na.rm=T) +
    (0.01*max(generank[which(generank != "Inf")], na.rm=T))
  generank <- generank[order(abs(generank), decreasing=T)]
  generank <- generank[!duplicated(names(generank))]
  generank <- generank[order(generank, decreasing=T)]
  return(generank) }


GSEA.enrich = function(x) {
  generank <- makeGeneRankStat(x$stat, toupper(x$hgnc_symbol))
  gsea.res <- fgseaMultilevel(pathways=hfb.genes, stats=generank, minSize=15,maxSize=500, eps = 0, nproc = 15)
  return(gsea.res)
}


hfb.genes = c(hfb.genes,Emunkitug)

res.gsea = GSEA.enrich(res)
# saveRDS(res.gsea, '~/HFB3010/Data/bulkRNAseq/GSEA.Preceding.IO.vs.NonIO.rds')



# BTLA ----
meta = read.csv('~/HFB200603/Data/bulkRNA/metadata_final.csv')

trt = read.csv('~/HFB200603/Data/Clinical/BTLA.metadat.csv')
trt = trt[,c('Subject', 'Preceding.Ther.MOA')]
trt$Preceding.Ther.MOA[which(trt$Preceding.Ther.MOA %in% c( 'C', 'C+T','T'))] = 'NonIO'
trt$Preceding.Ther.MOA[which(trt$Preceding.Ther.MOA != 'NonIO')] = 'IO'
colnames(trt) = c('Subject','PriorTx')

# meta = merge(meta, trt, by = 'Subject')
meta = merge(meta, trt, by = 'Subject')

# What is normalized ?? Batch normalized? Should I round to nearest integer?
dat = as.data.frame(fread('~/HFB200603/Data/bulkRNA/raw_mtx.tsv'))
# Remove outlier
# meta = subset(meta, Patient != '401008')

gene_id = dat$V1
dat = dat[,-which(colnames(dat) == 'V1')]

meta = subset(meta, TimePoint == 'Screening')
# meta = subset(meta, TimePoint == 'Screening' & Batch == 'Batch2')
# meta = subset(meta, !(Patient %in% c('607005','602008')))

dat = dat[,which(colnames(dat) %in% meta$Sample)]
meta = meta[match( colnames(dat), meta$Sample),]

dat = apply(dat,2, function(x) round(x))
dat = apply(dat,2, function(x) as.numeric(as.character(x)))
rownames(dat) = gene_id

dds <- DESeqDataSetFromMatrix(dat, meta, design = ~ PriorTx)


# normalize data
normalized.dat = vst(dds, blind = FALSE)

pca.dat = as.data.frame(assay(normalized.dat))
pca.dat <- prcomp(t(pca.dat))
rownames(pca.dat$x) = meta$Patient 
fviz_pca_ind(pca.dat,  habillage= paste0(meta$Batch), # this specifies the trait to draw a circle around 
             repel = T,
             addEllipses=TRUE, ellipse.level=0.5) + 
  theme(text = element_text(size = 23.5),
        axis.title = element_text(size = 17.5),
        axis.text = element_text(size = 17.5), legend.position = 'bottom') +
  guides(fill = guide_legend( override.aes = aes(label = ""))) 







dds <- DESeq(dds)
# specify comparisons. In this case comparing R vs NR
results <- results(dds, contrast = c('PriorTx', 'IO', 'NonIO'))
res = as.data.frame(results)
res$ensembl_gene_id = rownames(res)
res = merge(humBM,res, by = 'ensembl_gene_id')


hfb.genes = readRDS('~/GeneSets/HFB.Genesets.7.24.24.rds')

makeGeneRankStat <- function(statistic, names) {
  generank <- statistic
  names(generank) <- names
  generank <- generank[which(!is.na(generank))]
  generank <- generank[which(!is.nan(generank))]
  generank <- generank[which(names(generank) != "")]
  generank[which(generank == "-Inf")] <- min(generank[which(generank != "-Inf")], na.rm=T) -
    (0.01*min(generank[which(generank != "-Inf")], na.rm=T))
  generank[which(generank == "Inf")] <- max(generank[which(generank != "Inf")], na.rm=T) +
    (0.01*max(generank[which(generank != "Inf")], na.rm=T))
  generank <- generank[order(abs(generank), decreasing=T)]
  generank <- generank[!duplicated(names(generank))]
  generank <- generank[order(generank, decreasing=T)]
  return(generank) }


GSEA.enrich = function(x) {
  generank <- makeGeneRankStat(x$stat, toupper(x$hgnc_symbol))
  gsea.res <- fgseaMultilevel(pathways=hfb.genes, stats=generank, minSize=15,maxSize=500, eps = 0, nproc = 15)
  return(gsea.res)
}


hfb.genes = c(hfb.genes,Emunkitug)

res.gsea = GSEA.enrich(res)
# saveRDS(res.gsea, '~/HFB200603/Data/bulkRNA/GSEA.Preceding.IO.vs.NonIO.rds')




# TNFR2 -----
meta = read.csv('~/HFB2003/Data/bulkRNAseq/2025/metadata.csv')
meta$Sample = gsub('HFB200301_','',meta$Sample)

meta <- meta %>%
  dplyr::mutate(
    Patient = str_extract(Sample, "\\d{6}"),
    Timepoint = case_when(
      str_detect(Sample, "Screening|SCR") ~ "SCR",
      str_detect(Sample, "C2D8") ~ "C2D8",
      TRUE ~ NA_character_
    )
  )

rownames(meta) = meta$Sample

batch3 = c('202005','202005','402015','402015','403010','403010','403011','204002','401005','401005','402013','402013','403016','403016')

batch12 = c('401001','401001','402006','402006','402007','402012','402012','403008','403008','202001','202001','202002','202003','402004','402007')


tracker = read.csv('~/HFB2003/Data/Tracker/Tracker_2.4.25.csv')
tracker = tracker[,c('Subject','Tumor','Sched', 'Combo',   'Dose','BOR','Preceding.Tx.MoA','ToT.Month', 'TNFa', 'sTNFR2', 'TNFa.sTNFR2.hi', 'Prior.IO')]
tracker$BOR = gsub(c('Death|Clinical PD'),'PD',tracker$BOR)
tracker$BOR = gsub(c('SD>18wks'),'SD',tracker$BOR)
tracker$BOR = gsub(c('PD->SD'),'SD',tracker$BOR)
tracker$BOR = gsub(c('SD/PR'),'PR',tracker$BOR,fixed = T)
tracker$Condition = ifelse(tracker$Subject %in%  unique(subset(tracker, BOR %in% c('PR','SD'))$Subject), 'Responder','Non-Responder')
tracker$Patient = gsub('-','', tracker$Subject)

meta = merge(meta, tracker, by = 'Patient')

dat = as.data.frame(fread('~/HFB2003/Data/bulkRNAseq/2025/raw_mtx.tsv'))
colnames(dat) = gsub('HFB200301_','',colnames(dat))
# specify gene names
gene_id = dat$V1
dat = dat[,-which(colnames(dat) == 'V1')]


# Subset only data at baseline
meta = subset(meta, Timepoint == 'SCR' & Patient %in% batch3)
meta.scr = subset(meta, Timepoint == 'SCR')
dat = dat[,which(colnames(dat) %in% meta$Sample)]
meta = meta[match( colnames(dat), meta$Sample),]

meta$Preceding.Tx.MoA[which(meta$Preceding.Tx.MoA %in% c( 'C', 'C+T','T'))] = 'NonIO'
meta$Preceding.Tx.MoA[which(meta$Preceding.Tx.MoA != 'NonIO')] = 'IO'

# round data to nearest integer and make sure it's in numerical format
# dat = apply(dat,2, function(x) round(x))
# dat = apply(dat,2, function(x) as.numeric(as.character(x)))
# add gene ids 
rownames(dat) = gene_id
# put data into DESeq2 object
dds <- DESeqDataSetFromMatrix(dat, meta, design = ~ Preceding.Tx.MoA)


normalized.dat = vst(dds, blind = FALSE)

pca.dat = as.data.frame(assay(normalized.dat))
pca.dat <- prcomp(t(pca.dat))
rownames(pca.dat$x) = meta$Patient 
fviz_pca_ind(pca.dat,  habillage= paste0(meta$Preceding.Tx.MoA), # this specifies the trait to draw a circle around 
             repel = T,
             addEllipses=TRUE, ellipse.level=0.5) + 
  theme(text = element_text(size = 23.5),
        axis.title = element_text(size = 17.5),
        axis.text = element_text(size = 17.5), legend.position = 'bottom') +
  guides(fill = guide_legend( override.aes = aes(label = ""))) 



dds <- DESeq(dds)
# specify comparisons. In this case comparing R vs NR
results <- results(dds, contrast = c('Preceding.Tx.MoA', 'IO', 'NonIO'))
res = as.data.frame(results)
res$ensembl_gene_id = rownames(res)
res = merge(humBM,res, by = 'ensembl_gene_id')

hfb.genes = readRDS('~/GeneSets/HFB.Genesets.7.24.24.rds')

makeGeneRankStat <- function(statistic, names) {
  generank <- statistic
  names(generank) <- names
  generank <- generank[which(!is.na(generank))]
  generank <- generank[which(!is.nan(generank))]
  generank <- generank[which(names(generank) != "")]
  generank[which(generank == "-Inf")] <- min(generank[which(generank != "-Inf")], na.rm=T) -
    (0.01*min(generank[which(generank != "-Inf")], na.rm=T))
  generank[which(generank == "Inf")] <- max(generank[which(generank != "Inf")], na.rm=T) +
    (0.01*max(generank[which(generank != "Inf")], na.rm=T))
  generank <- generank[order(abs(generank), decreasing=T)]
  generank <- generank[!duplicated(names(generank))]
  generank <- generank[order(generank, decreasing=T)]
  return(generank) }


GSEA.enrich = function(x) {
  generank <- makeGeneRankStat(x$stat, toupper(x$hgnc_symbol))
  gsea.res <- fgseaMultilevel(pathways=hfb.genes, stats=generank, minSize=15,maxSize=500, eps = 0, nproc = 15)
  return(gsea.res)
}


hfb.genes = c(hfb.genes,Emunkitug)

res.gsea = GSEA.enrich(res)
# saveRDS(res.gsea, '~/HFB2003/Data/bulkRNAseq/2025/GSEA.Preceding.IO.vs.NonIO.Batch3.rds')


# IO GSEA -----
io.ox40 = readRDS('~/HFB3010/Data/bulkRNAseq/GSEA.Preceding.IO.vs.NonIO.rds')
io.btla = readRDS('~/HFB200603/Data/bulkRNA/GSEA.Preceding.IO.vs.NonIO.rds')
io.tnf2 = readRDS('~/HFB2003/Data/bulkRNAseq/2025/GSEA.Preceding.IO.vs.NonIO.Batch3.rds')

io.tnf2$Condition = 'Emunkitug Ph1'
io.btla$Condition = 'HFB200603 Ph1'
io.ox40$Condition = 'HFB301001 Ph1'
io.all = as.data.frame(rbind(io.tnf2,io.btla,io.ox40))

sort(table(subset(io.all, pval < .05 & NES < 0)$pathway))

pathways =c('PHONG_TNF_TARGETS_UP','NK_Cell_Functions','T_Cell_Functions','NF_kappaB_Signaling', 'TNF_Family_Signaling','Lymphocyte_Activation','Cytotoxicity', 'T_Cell_Functions','Adaptive_Immune_System','Emunkitug','Cytokines','Costimulatory_Signaling','Antigen Presentation')

hfb.gsea = subset(io.all, pathway %in% pathways)

hfb.gsea$pathway = gsub('PHONG_','',hfb.gsea$pathway )
hfb.gsea$pathway = gsub('Emunkitug','Emunkitug_Response_Signature_v2',hfb.gsea$pathway )
hfb.gsea$pathway = str_to_title(hfb.gsea$pathway)
hfb.gsea$pathway = gsub('Tnf','TNF',hfb.gsea$pathway )
hfb.gsea$pathway = gsub('Nk','NK',hfb.gsea$pathway )
hfb.gsea$pathway = gsub('Nf','NF',hfb.gsea$pathway )
hfb.gsea$pathway = gsub('kappab','kappaB',hfb.gsea$pathway )


hfb.gsea$pathway = factor(hfb.gsea$pathway, levels = rev(c("Emunkitug_response_signature_v2",
                                                       "TNF_targets_up", "TNF_family_signaling" , 'NF_kappaB_signaling',  "Adaptive_immune_system",
                                                       "Cytokines","Cytotoxicity" ,
                                                       "Antigen Presentation","Costimulatory_signaling",
                                                       "Lymphocyte_activation","NK_cell_functions", "T_cell_functions")))

hfb.gsea$Condition = factor(hfb.gsea$Condition, levels = c('Emunkitug Ph1','HFB200603 Ph1', 'HFB301001 Ph1'))



# hfb.gsea$NES = ifelse(hfb.gsea$NES > 0 & hfb.gsea$padj > .05, 0, hfb.gsea$NES)

hfb.gsea$pval_label <- case_when(
  hfb.gsea$pval < 0.0001 ~ "****",
  hfb.gsea$pval < 0.001 ~ "***",
  hfb.gsea$pval < 0.01 ~ "**",
  hfb.gsea$pval < 0.05 ~ "*",
  TRUE ~ ""
)

color_fun <- scales::col_numeric(
  palette = c("blue", "white", "red"),
  domain = c(min(hfb.gsea$NES), 0, max(hfb.gsea$NES)))

ggplot(hfb.gsea , aes(x = Condition, y = pathway)) +
  geom_tile(aes(fill = NES), color = "white") +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", 
                       midpoint = 0,
                       limits = c(-3, 3),
                       breaks = c(-2, 0, 2),
                       labels = scales::label_number()) +
  geom_text(aes(label = pval_label), color = "black", size = 5) +
  theme_minimal(base_size = 20) +
  labs(x = "", y = "", title = 'GSEA: Preceding IO vs Non-IO at Baseline') +
  rotate_x_text( angle = 45)





















# Preced Tx Boxplot -----


###### 1) Emunkitug ----
genes = readRDS('~/GeneSets/HFB.Genesets.2.24.25.rds')
genes.of.interest = genes$Emunkitug

meta = read.csv('~/HFB2003/Data/bulkRNAseq/2025/Metadata.2.26.25.csv')
dat = as.data.frame(fread('~/HFB2003/Data/bulkRNAseq/2025/normalized_mtx.tsv'))
colnames(dat) = gsub('HFB200301_','',colnames(dat))

meta.scr = subset(meta, Timepoint == 'SCR')

gene_id = dat$V1
dat = dat[,-which(colnames(dat) == 'V1')]
#colnames(dat) = substr(colnames(dat), 1, nchar(colnames(dat)) - 5)
#colnames(dat) = gsub("([^_]+)_\\d{2,3}M_(.*)", "\\1_\\2", colnames(dat))

dat = dat[,which(colnames(dat) %in% meta.scr$Sample)]
meta.scr = meta.scr[match( colnames(dat), meta.scr$Sample),]

rownames(dat) = gene_id
norm.mat = as.data.frame(dat)
norm.mat$ensembl_gene_id = rownames(norm.mat)
norm.mat = merge(humBM, norm.mat, by = 'ensembl_gene_id')
norm.mat = subset(norm.mat, hgnc_symbol %in% genes.of.interest)

rownames(norm.mat) = norm.mat$hgnc_symbol
norm.mat = norm.mat[,-which(colnames(norm.mat) %in% colnames(humBM))]
colnames(norm.mat) = substr(colnames(norm.mat),1,6)

norm.mat.scale.screen = t(apply(norm.mat, 1, function(x) scale(x)))
colnames(norm.mat.scale.screen) = colnames(norm.mat)


emun.scr.z = norm.mat.scale.screen[which(rownames(norm.mat.scale.screen) %in% genes.of.interest),]
emun.scr.z = data.frame( SampleID = names(colMeans(emun.scr.z)), PD.Signature = colMeans(emun.scr.z))
emun.scr.z = merge(meta.scr, emun.scr.z, by.x = 'Patient', by.y = 'SampleID')


emun.scr.z$Preceding.Tx.MoA = gsub('P+T','IO',emun.scr.z$Preceding.Tx.MoA, fixed = T)
emun.scr.z$Preceding.Tx.MoA = gsub('T / P','IO',emun.scr.z$Preceding.Tx.MoA, fixed = T)
emun.scr.z$Preceding.Tx.MoA = gsub('P','IO',emun.scr.z$Preceding.Tx.MoA, fixed = T)


emun.scr.z$Preceding.Tx.MoA = factor(emun.scr.z$Preceding.Tx.MoA, levels = c('C','T','C+T','IO'))
ggplot(emun.scr.z, aes(x = Preceding.Tx.MoA, y = PD.Signature, color = Preceding.Tx.MoA)) +
  geom_boxplot(outliers = F) + geom_point(size =4, color = 'black',position = position_dodge2(.25) )+
  theme_classic(base_size = 30) + theme(legend.position = 'none') +
  labs(y = 'Emunkitug Response Signature (Z)', x= 'Preceding Tx', title = 'Emunkitug Ph1 Baseline') +
  geom_label_repel(aes(label = Tumor),position = position_dodge2(.25), size =4, color = 'grey25', alpha = .5  ) +
  geom_hline(yintercept = 0, linetype = 'dashed', linewidth = 1.2) +
  rotate_x_text(angle = 65) 



###### 2) BTLA -----
meta = read.csv('~/HFB200603/Data/bulkRNA/metadata_final.csv')

trt = read.csv('~/HFB200603/Data/Clinical/BTLA.metadat.csv')
trt = trt[,c('Subject', 'Preceding.Ther.MOA')]
# trt$Preceding.Ther.MOA[which(trt$Preceding.Ther.MOA %in% c( 'C', 'C+T','T'))] = 'NonIO'
# trt$Preceding.Ther.MOA[which(trt$Preceding.Ther.MOA != 'NonIO')] = 'IO'
colnames(trt) = c('Subject','PriorTx')

# meta = merge(meta, trt, by = 'Subject')
meta = merge(meta, trt, by = 'Subject')

# What is normalized ?? Batch normalized? Should I round to nearest integer?
dat = as.data.frame(fread('~/HFB200603/Data/bulkRNA/raw_mtx.tsv'))
# Remove outlier
# meta = subset(meta, Patient != '401008')

gene_id = dat$V1
dat = dat[,-which(colnames(dat) == 'V1')]

meta = subset(meta, TimePoint == 'Screening' & Batch == 'Batch2')
# meta = subset(meta, TimePoint == 'Screening')
meta = subset(meta, !(Patient %in% c('607005','602008')))

dat = dat[,which(colnames(dat) %in% meta$Sample)]
meta = meta[match( colnames(dat), meta$Sample),]

dat = apply(dat,2, function(x) round(x))
dat = apply(dat,2, function(x) as.numeric(as.character(x)))
rownames(dat) = gene_id

dds <- DESeqDataSetFromMatrix(dat, meta, design = ~ PriorTx)
# normalize data
normalized.dat = vst(dds, blind = FALSE)
df = as.data.frame(assay(normalized.dat))
df$ensembl_gene_id = rownames(df)
df = merge(humBM, df, by = 'ensembl_gene_id')
df = subset(df, hgnc_symbol %in% genes.of.interest)
rownames(df) = df$hgnc_symbol
df = df[,-which(colnames(df) %in% colnames(humBM))]


norm.mat.scale.screen = t(apply(df, 1, function(x) scale(x)))
colnames(norm.mat.scale.screen) = colnames(df)


btla.scr.z = norm.mat.scale.screen[which(rownames(norm.mat.scale.screen) %in% genes.of.interest),]
btla.scr.z = data.frame( SampleID = names(colMeans(btla.scr.z)), PD.Signature = colMeans(btla.scr.z))
btla.scr.z = merge(meta, btla.scr.z, by.x = 'Sample', by.y = 'SampleID')

btla.scr.z$PriorTx = gsub('P+T','IO',btla.scr.z$PriorTx, fixed = T)
# btla.scr.z$PriorTx = gsub('T / P','IO',btla.scr.z$PriorTx, fixed = T)
btla.scr.z$PriorTx = gsub('P','IO',btla.scr.z$PriorTx, fixed = T)

btla.scr.z$PriorTx = factor(btla.scr.z$PriorTx, levels = c('C','T','C+T','IO'))
ggplot(btla.scr.z, aes(x = PriorTx, y = PD.Signature, color = PriorTx)) +
  geom_boxplot(outliers = F) + geom_point(size =4, color = 'black',position = position_dodge2(.25) )+
  theme_classic(base_size = 30) + theme(legend.position = 'none') +
  labs(y = 'Emunkitug Response Signature (Z)', x= 'Preceding Tx', title = 'HFB200603 Ph1 Baseline',subtitle = 'Batch2') +
  # geom_label_repel(aes(label = Tumor),position = position_dodge2(.25), size =4, color = 'grey25', alpha = .5  ) +
  geom_hline(yintercept = 0, linetype = 'dashed', linewidth = 1.2) +
  rotate_x_text(angle = 65) 



##### 3) OX40 ----
meta = read.csv('~/HFB3010/Data/bulkRNAseq/metadata_final.csv')

trt = read_xlsx('~/HFB3010/Data/Clinical/OX40 Clinical Data_16AUG2024.xlsx')
trt = trt[,c('Subject', 'Preceding.Therapy')]
# trt$Preceding.Therapy[which(trt$Preceding.Therapy %in% c('C','T','P+T'))] = 'NonIO'
# trt$Preceding.Therapy[which(trt$Preceding.Therapy %in% c('P'))] = 'IO'
colnames(trt) = c('Subject','PriorTx')

# meta = merge(meta, trt, by = 'Subject')
meta = merge(meta, trt, by = 'Subject')

# What is normalized ?? Batch normalized? Should I round to nearest integer?
dat = as.data.frame(fread('~/HFB3010/Data/bulkRNAseq/normalized_mtx.tsv'))
# Remove outlier
# meta = subset(meta, Patient != '401008')

gene_id = dat$gene_id
dat = dat[,-which(colnames(dat) == 'gene_id')]

meta = subset(meta, TimePoint == 'Screening' & Batch == 'Batch3' & !(Subject %in% c('301-004','303-001')))
dat = dat[,which(colnames(dat) %in% meta$Sample)]
meta = meta[match( colnames(dat), meta$Sample),]

dat = apply(dat,2, function(x) round(x))
dat = apply(dat,2, function(x) as.numeric(as.character(x)))
rownames(dat) = gene_id

dds <- DESeqDataSetFromMatrix(dat, meta, design = ~ PriorTx)



normalized.dat = vst(dds, blind = FALSE)
df = as.data.frame(assay(normalized.dat))
df$ensembl_gene_id = rownames(df)
df = merge(humBM, df, by = 'ensembl_gene_id')
df = subset(df, hgnc_symbol %in% genes.of.interest)
rownames(df) = df$hgnc_symbol
df = df[,-which(colnames(df) %in% colnames(humBM))]


norm.mat.scale.screen = t(apply(df, 1, function(x) scale(x)))
colnames(norm.mat.scale.screen) = colnames(df)

ox40.scr.z = norm.mat.scale.screen[which(rownames(norm.mat.scale.screen) %in% genes.of.interest),]
ox40.scr.z = data.frame( SampleID = names(colMeans(ox40.scr.z)), PD.Signature = colMeans(ox40.scr.z))
ox40.scr.z = merge(meta, ox40.scr.z, by.x = 'Sample', by.y = 'SampleID')

ox40.scr.z$PriorTx = gsub('P+T','IO',ox40.scr.z$PriorTx, fixed = T)
# ox40.scr.z$PriorTx = gsub('T / P','IO',ox40.scr.z$PriorTx, fixed = T)
ox40.scr.z$PriorTx = gsub('P','IO',ox40.scr.z$PriorTx, fixed = T)

ox40.scr.z$PriorTx = factor(ox40.scr.z$PriorTx, levels = c('C','T','C+T','IO'))
ggplot(ox40.scr.z, aes(x = PriorTx, y = PD.Signature, color = PriorTx)) +
  geom_boxplot(outliers = F) + geom_point(size =4, color = 'black',position = position_dodge2(.25) )+
  theme_classic(base_size = 30) + theme(legend.position = 'none') +
  labs(y = 'Emunkitug Response Signature (Z)', x= 'Preceding Tx', title = 'Nuvustotug Ph1 Baseline',subtitle = 'Batch3') +
  # geom_label_repel(aes(label = Tumor),position = position_dodge2(.25), size =4, color = 'grey25', alpha = .5  ) +
  geom_hline(yintercept = 0, linetype = 'dashed', linewidth = 1.2) +
  rotate_x_text(angle = 65) 


# RNA/mIF Cor -----


##### 1) TNFR2 -----
df.log2 = read.csv('~/HFB2003/Data/bulkRNAseq/2025/Emunkitug.bulkRNAseq.log2.counts.2.21.25.csv')
df.meta = read.csv('~/HFB2003/Data/bulkRNAseq/2025/Metadata.2.26.25.csv')
df.meta = subset(df.meta, Timepoint == 'SCR')

df= merge(df.meta, df.log2, by = 'Subject')
df = subset(df, Timepoint.y == 'SCR' )

ggplot(df, aes(x = TNFRSF1B, y = sTNFR2) ) +
  geom_point(aes(color = BOR), size = 4) +
  geom_smooth(method = 'lm', color = 'black', se = F) +
  stat_cor(size = 8) +
  labs(x = 'bulk RNA TNFR2 (log2)', y = 'sTNFR2 (pg/mL)', title = 'Baseline', color = '') +
  scale_color_lancet() + 
  theme_classic(base_size = 30) + theme(legend.position = 'top') 



mif = read.csv('~/HFB2003/Data/mIF/2025/Alex/mIF.CellType.TNFR2.expression.1.10.25.csv')
mif = subset(mif, Visit == 'SCR' & CellType %in% c('CD8T','Macrophage','CD4T','Treg','Tumor'))

df.rna.mif = merge(subset(df.log2, Timepoint == 'SCR'), mif, by.x = 'Subject', by.y = 'SubjectID')
df.rna.mif = df.rna.mif[-which(df.rna.mif$Log2.Density == '-Inf'),]
df.rna.mif = merge(df.rna.mif, df.meta, by = 'Subject')

ggplot(df.rna.mif, aes(x = TNFRSF1B, y = Log2.Density) ) +
  geom_point(aes(color = BOR), size = 4) +
  geom_smooth(method = 'lm', color = 'black', se = F) +
  stat_cor(size = 6) +
  labs(x = 'bulk RNA TNFR2 (log2)', y = 'mIF TNFR2 density (log2)', title = 'Baseline', color = '') +
  scale_color_lancet() + 
  facet_wrap(~CellType,scales = 'free_y', nrow = 2) +
  theme_minimal(base_size = 25) + theme(legend.position = 'top') 



##### 2) CD8A -----
df.log2 = read.csv('~/HFB2003/Data/bulkRNAseq/2025/Emunkitug.bulkRNAseq.log2.counts.2.21.25.csv')
df.meta = read.csv('~/HFB2003/Data/bulkRNAseq/2025/Metadata.2.26.25.csv')
df.meta = subset(df.meta, Timepoint == 'SCR')

mif = read.csv('~/HFB2003/Data/mIF/2025/Alex/mIF.CellTypeDensity.1.10.25.csv')
mif = subset(mif, Visit == 'SCR' & CellType %in% c('CD8T'))

df.rna.mif = merge(subset(df.log2, Timepoint == 'SCR'), mif, by.x = 'Subject', by.y = 'SubjectID')
df.rna.mif = df.rna.mif[-which(df.rna.mif$Log2.Density == '-Inf'),]
df.rna.mif = merge(df.rna.mif, df.meta, by = 'Subject')

ggplot(df.rna.mif, aes(x = CD8A, y = Log2.Density) ) +
  geom_point(aes(color = BOR), size = 4) +
  geom_smooth(method = 'lm', color = 'black', se = F) +
  stat_cor(size = 6) +
  labs(x = 'bulk RNA CD8A (log2)', y = 'mIF CD8 T-cell density (log2)', title = 'Baseline', color = '') +
  scale_color_lancet() + 
  # facet_wrap(~CellType,scales = 'free_y', nrow = 2) +
  theme_minimal(base_size = 25) + theme(legend.position = 'top') 



##### 3) PD -----


###### a) TNFR2 / PD1 expression ------

# mif = read.csv('C:/Users/SpencerHuggett/OneDrive - HifiBiO/Documents/HFB2003/Data/mIF/2025/HFB200301_mif_summary.1.8.25.csv')
mif = read.csv('~/HFB2003/Data/mIF/2025/HFB200301_mif_summary.1.8.25.csv')

# resp.meta = read.csv('C:/Users/SpencerHuggett/OneDrive - HifiBiO/Documents/HFB2003/Data/Dose.Selection/TNFR2.Efficacy.csv')
resp.meta = read.csv('~/HFB2003/Data/Dose.Selection/TNFR2.Efficacy.csv')
resp.meta = resp.meta[c('Subject','TumorSTD','Dose','Sched','Combo')]
resp.meta$Dose = gsub(' ','', resp.meta$Dose )
resp.meta$Dose = factor(resp.meta$Dose, levels = c('5mg','15mg','50mg','150mg','300mg'))
resp.meta$TumorSTD = gsub('Soft Tissue Sarcoma','STS',resp.meta$TumorSTD)
resp.meta$TumorSTD = gsub(c('Testicular Germ Cell |Testicular Germ Cell'),'Testicular',resp.meta$TumorSTD)
# resp.meta$TumorSTD = factor(resp.meta$TumorSTD, levels = c('Cervical','HNSCC','NSCLC','Mesothelioma','RCC','Gastric', 'STS','Melanoma','Testicular'))


mif$TNFR2 = 'Other'
mif$TNFR2[grep('TNFR2',mif$phenotypes)] = 'TNFR2 expression'

mif$PD1 = 'Other'
mif$PD1[grep('PD1',mif$phenotypes)] = 'PD1 expression'

mif$PDL1 = 'Other'
mif$PDL1[grep('PDL1',mif$phenotypes)] = 'PDL1 expression'

mif.tnfr2 = subset(mif, value != 'Inf') %>%
  dplyr::group_by(SubjectID, Visit,CellType, TNFR2) %>%
  dplyr::summarise(TNFR2.Log2.Density = log2(sum(value , na.rm =T )))

mif.pd1 = subset(mif, value != 'Inf') %>%
  dplyr::group_by(SubjectID, Visit,CellType, PD1) %>%
  dplyr::summarise(PD1.Log2.Density = log2(sum(value , na.rm =T )))

mif.pdl1 = subset(mif, value != 'Inf') %>%
  dplyr::group_by(SubjectID, Visit,CellType, PDL1) %>%
  dplyr::summarise(PDL1.Log2.Density = log2(sum(value , na.rm =T )))

mif.pd1 = subset(mif.pd1, PD1 == 'PD1 expression')
mif.pdl1 = subset(mif.pdl1, PDL1 == 'PDL1 expression')
mif.tnfr2 = subset(mif.tnfr2, TNFR2 == 'TNFR2 expression')

mif.tnfr2$ID =paste0(mif.tnfr2$Visit, mif.tnfr2$SubjectID,mif.tnfr2$CellType)
mif.pd1$ID =  paste0(mif.pd1$Visit, mif.pd1$SubjectID,mif.pd1$CellType)
mif.pdl1$ID = paste0(mif.pdl1$Visit, mif.pdl1$SubjectID,mif.pdl1$CellType)

mif.pd1 = mif.pd1[,c('ID','PD1.Log2.Density')]
mif.pdl1 = mif.pdl1[,c('ID','PDL1.Log2.Density')]

mif.marker = mif.tnfr2 %>%
  left_join(mif.pd1,  by='ID') %>%
  left_join(mif.pdl1, by='ID')

mif.marker.scr = subset(mif.marker, Visit == 'SCR')
mif.marker.c2d = subset(mif.marker, Visit == 'C2D8')
mif.marker.c2d = mif.marker.c2d[,c('ID','TNFR2.Log2.Density','PD1.Log2.Density','PDL1.Log2.Density')]
colnames(mif.marker.c2d) = c('ID','TNFR2.C2D8','PD1.C2D8','PDL1.C2D8')

mif.marker.scr$ID = gsub(c('SCR|C2D8'),'',mif.marker.scr$ID )
mif.marker.c2d$ID = gsub(c('SCR|C2D8'),'',mif.marker.c2d$ID )

mif.pre.post = merge(mif.marker.scr, mif.marker.c2d, by = 'ID')
mif.pre.post$TNFR2.Delta = mif.pre.post$TNFR2.C2D8 - mif.pre.post$TNFR2.Log2.Density
mif.pre.post$PD1.Delta = mif.pre.post$PD1.C2D8 - mif.pre.post$PD1.Log2.Density
mif.pre.post$PDL1.Delta = mif.pre.post$PDL1.C2D8 - mif.pre.post$PDL1.Log2.Density


df.log2 = read.csv('~/HFB2003/Data/bulkRNAseq/2025/Emunkitug.bulkRNAseq.log2.counts.2.21.25.csv')
df.scr = subset(df.log2, Timepoint == 'SCR')
df.c2d = subset(df.log2, Timepoint == 'C2D8')
df.pair = merge(df.scr, df.c2d, by = 'Subject')

df.pair$TNFR2.Delta.RNA = df.pair$TNFRSF1B.y - df.pair$TNFRSF1B.x
df.pair$PD1.Delta.RNA = df.pair$PDCD1.y - df.pair$PDCD1.x
df.pair$PDL1.Delta.RNA = df.pair$CD274.y  - df.pair$CD274.x
df.pair$CD8A.Delta.RNA = df.pair$CD8A.y  - df.pair$CD8A.x

df.pair = df.pair[,grep(c('Subject|Delta'), colnames(df.pair))]

df.meta = read.csv('~/HFB2003/Data/bulkRNAseq/2025/Metadata.2.26.25.csv')
df.meta = subset(df.meta, Timepoint == 'SCR')

df= merge(df.meta, df.pair, by = 'Subject')

df.rna.mif = merge(df, mif.pre.post, by.x = 'Subject',by.y = 'SubjectID')

df.rna.mif.cells = subset(df.rna.mif , CellType %in% c('CD8T','CD4T','Treg','Tumor'))

ggplot(subset(df.rna.mif.cells, !(TNFR2.Delta %in% c('Inf','-Inf'))), aes(x = TNFR2.Delta.RNA, y = TNFR2.Delta) ) +
  geom_point(aes(color = BOR), size = 4) +
  geom_smooth(method = 'lm', color = 'black', se = F) +
  stat_cor(size = 6) +
  geom_hline(yintercept = 0, linetype = 'dashed',linewidth = 1.3)+ 
  geom_vline(xintercept = 0, linetype = 'dashed',linewidth = 1.3)+ 
  labs(x = expression(Delta~RNA~TNFR2~Log[2]) , y = expression(Delta~mIF~TNFR2~Density~Log[2]), title = 'Emunkitug PD Effects on TNFR2 Expression', color = '') +
  scale_color_lancet() + 
  facet_wrap(~CellType,scales = 'free_y', nrow = 1) +
  theme_minimal(base_size = 25) + theme(legend.position = 'top') 


ggplot(subset(df.rna.mif.cells, CellType != 'Tumor' & !(PD1.Delta %in% c('Inf','-Inf'))), aes(x = PD1.Delta.RNA, y = PD1.Delta) ) +
  geom_point(aes(color = BOR), size = 4) +
  geom_smooth(method = 'lm', color = 'black', se = F) +
  stat_cor(size = 6) +
  geom_hline(yintercept = 0, linetype = 'dashed',linewidth = 1.3)+ 
  geom_vline(xintercept = 0, linetype = 'dashed',linewidth = 1.3)+ 
  labs(x = expression(Delta~RNA~PD1~Log[2]) , y = expression(Delta~mIF~PD1~Density~Log[2]), title = 'Emunkitug PD Effects on PD1 Expression', color = '') +
  scale_color_lancet() + 
  facet_wrap(~CellType,scales = 'free_y', nrow = 1) +
  theme_minimal(base_size = 25) + theme(legend.position = 'top') 




ggplot(subset(df.rna.mif.cells, !(PDL1.Delta %in% c('Inf','-Inf'))), aes(x = PDL1.Delta.RNA, y = PDL1.Delta) ) +
  geom_point(aes(color = BOR), size = 4) +
  geom_smooth(method = 'lm', color = 'black', se = F) +
  stat_cor(size = 6) +
  geom_hline(yintercept = 0, linetype = 'dashed',linewidth = 1.3)+ 
  geom_vline(xintercept = 0, linetype = 'dashed',linewidth = 1.3)+ 
  labs(x = expression(Delta~RNA~PDL1~Log[2]) , y = expression(Delta~mIF~PDL1~Density~Log[2]), title = 'Emunkitug PD Effects on PDL1 Expression', color = '') +
  scale_color_lancet() + 
  facet_wrap(~CellType,scales = 'free_y', nrow = 1) +
  theme_minimal(base_size = 25) + theme(legend.position = 'top') 


