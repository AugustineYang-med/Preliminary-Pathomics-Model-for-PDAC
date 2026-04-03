########################### Differential analysis using edgeR ##################
library(readr)
library(data.table)  
library(dplyr)       
library(ggplot2)     
library(pheatmap)    
library(DESeq2)    
library(limma)   
library(tinyarray)   
library(openxlsx)

library(edgeR)

identical(colnames(tcga_paad_exp), rownames(colData))
group <- colData$group
class(group)
d <- DGEList(counts = tcga_paad_exp, group = group)
keep <- rowSums(cpm(d) > 1) >= 2
keep <- filterByExpr(d)

table(keep)

d <- d[keep, , keep.lib.sizes = FALSE]
d$samples$lib.size <- colSums(d$counts)
d <- calcNormFactors(d)

head(d$samples)
#                     group lib.size norm.factors
# TCGA-2J-AAB1-01  Low_risk 44668882    1.0066464
# TCGA-2J-AAB4-01  Low_risk 50555052    1.0941152
# TCGA-2J-AAB6-01  Low_risk 39274565    0.8879193
# TCGA-2J-AAB8-01 High_risk 22335107    1.1227573
# TCGA-2J-AAB9-01 High_risk 21818885    0.9881490
# TCGA-2J-AABA-01 High_risk 33357374    1.1312226
# Assign normalized data to dge variable
dge = d

design <- model.matrix(~0 + factor(group))
rownames(design) <- colnames(dge)
colnames(design) <- levels(factor(group))

dge <- estimateGLMCommonDisp(dge, design)
dge <- estimateGLMTrendedDisp(dge, design)
dge <- estimateGLMTagwiseDisp(dge, design)

fit <- glmFit(dge, design)
lrt <- glmLRT(fit, contrast = c(-1, 1))

nrDEG <- topTags(lrt, n = nrow(dge))

DEG_edgeR <- as.data.frame(nrDEG)

head(DEG_edgeR)
#             logFC   logCPM       LR       PValue          FDR
# STAR    -4.531300 2.522249 99.55831 1.904718e-23 4.100477e-19
# CYP21A2 -3.972406 1.501827 79.52857 4.753011e-19 5.116141e-15
# NXPH4   -2.921580 2.516680 75.25021 4.146870e-18 2.975794e-14
# FGB     -4.027070 5.775333 67.59081 2.012035e-16 1.082877e-12
# ELANE   -3.807738 1.742607 62.50358 2.659612e-15 1.145122e-11
# PPBP    -4.268197 3.089323 62.02345 3.393910e-15 1.217735e-11

save(DEG_edgeR, file = './Result/DFS_DEGs/Data/DEG_edgeR.Rdata')

logFC = 1
P.Value = 0.05
k1 <- (DEG_edgeR$FDR < P.Value) & (DEG_edgeR$logFC < -logFC)
k2 <- (DEG_edgeR$FDR < P.Value) & (DEG_edgeR$logFC > logFC)
DEG_edgeR <- mutate(DEG_edgeR, change = ifelse(k1, "down", ifelse(k2, "up", "stable")))
table(DEG_edgeR$change)

############## Gradient color volcano plot ################
DEG_edgeR$genename <- rownames(DEG_edgeR)
DEG_edgeR$color <- -log10(DEG_edgeR$FDR)
pdf(file = "./Result/DFS_DEGs/Plot/DEG_edgeR_Volcano Plot.pdf", width = 6, height = 6)
gradual_volcano(DEG_edgeR, x = "logFC", y = "FDR", pointSizeRange = c(0.5, 4),
                label = "genename", 
                # custom_label = c("STAR", "CYP21A2", "NXPH4", "FGB", "ELANE", "CYP4F22", "RPL37P6", "SNORD17", "PLA2G2A", "PAPPA2"),
                label_number = 10,
                fills = c("#0D0D4C", "#4063AF","#39bbec", "#7DCA99", "#FAED31","#f38466","#E51D6C", "#b81f25"),
                colors = c("#17194e", "#1E4D8B", "#1E9A78", "#49B8B0", "#F48B7C","#a22f27", "#211f1f"),
                log2FC_cut = 1, FDR_cut = 0.05, 
                output = FALSE)+ scale_x_continuous(limits = c(-8, 8))+ scale_y_continuous(limits = c(0, 35))
dev.off()
