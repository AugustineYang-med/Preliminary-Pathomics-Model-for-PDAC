library(readr)
library(data.table)  # For efficient large dataset processing
library(dplyr)       # For data manipulation and transformation
library(ggplot2)     # For plotting
library(pheatmap)    # For heatmap visualization
library(DESeq2)      # Differential analysis method 1
library(edgeR)       # Differential analysis method 2
library(limma)       # Differential analysis method 3
library(tinyarray)   # For various plots, used here for Venn diagrams
library(openxlsx)

tcga_gene_expression <- read_tsv("./Data/TCGA Gene Expression/TCGA-PAAD.star_counts.tsv.gz",
                        col_types = cols(),  # Use cols() for automatic type guessing
                        progress = TRUE)     # Show read progress bar

group <- read.xlsx("./Result/Risk Score/tcga DFS risk group.xlsx")

probemap_data <- read.delim("./Data/TCGA Gene Expression/gencode.v36.annotation.gtf.gene.probemap",
                            header = TRUE,          # First row contains column names
                            sep = "\t",             # Tab-separated
                            stringsAsFactors = FALSE)

id <- colnames(tcga_gene_expression)
id <- id[id != "Ensembl_ID"]
id <- substr(id, 1, 15)
col <- c("Ensembl_ID", id)
colnames(tcga_gene_expression) <- col

cols <- c("gene", colnames(tcga_gene_expression))
genename <- subset(probemap_data, select = c("id","gene"))
tcga_gene_expression <- merge(genename, tcga_gene_expression, by.x = "id", by.y = "Ensembl_ID")
# Remove duplicates
tcga_gene_expression <- distinct(tcga_gene_expression, gene, .keep_all = T)
rownames(tcga_gene_expression) <- tcga_gene_expression$gene
tcga_gene_expression <- subset(tcga_gene_expression, select = -c(id, gene))
genename <- rownames(tcga_gene_expression)
tcga_gene_expression <- as.data.frame(apply(tcga_gene_expression, 2, as.numeric))
rownames(tcga_gene_expression) <- genename
tcga_paad_exp <- 2^tcga_gene_expression-1
tcga_paad_exp <- round(tcga_paad_exp)
tcga_paad_exp <- as.data.frame(t(tcga_paad_exp))
tcga_paad_exp$sampleid <- rownames(tcga_paad_exp)
group <- subset(group, select = c(id, group))
tcga_paad_exp <- merge(group, tcga_paad_exp, by.x = "id", by.y = "sampleid")

group <- subset(tcga_paad_exp, select = c(id, group))
rownames(group) <- group$id
group <- subset(group, select = -id)
group$group <- factor(group$group, levels = c("Low Risk", "High Risk"))

rownames(tcga_paad_exp) <- tcga_paad_exp$id
tcga_paad_exp <- subset(tcga_paad_exp, select = -c(id, group))
tcga_paad_exp <- as.data.frame(t(tcga_paad_exp))
identical(rownames(group), colnames(tcga_paad_exp))
colData <- group

### Differential analysis using DESeq2
dds <- DESeqDataSetFromMatrix(countData = tcga_paad_exp, # Expression matrix
                              colData = colData,           # Mapping between expression matrix columns and group info
                              design = ~group)           # group refers to the group column in the colData data frame
head(dds)
dds <- DESeq(dds) # Perform differential analysis
resultsNames(dds) # View result names

res <- results(dds, contrast = c("group", "High Risk", "Low Risk"))
res <- res[order(res$padj), ] # Sort by adjusted p-value
DEG <- as.data.frame(res)

DEG_deseq2 <- na.omit(DEG)

# Add change column to label up/down-regulated genes; adjust thresholds as needed
logFC = 1
P.Value = 0.05
k1 <- (DEG_deseq2$padj < P.Value) & (DEG_deseq2$log2FoldChange < -logFC)
k2 <- (DEG_deseq2$padj < P.Value) & (DEG_deseq2$log2FoldChange > logFC)
DEG_deseq2 <- mutate(DEG_deseq2, change = ifelse(k1, "down", ifelse(k2, "up", "stable")))
table(DEG_deseq2$change)

save(tcga_gene_expression, tcga_paad_exp, colData, DEG_deseq2, file = "./Result/DFS_DEGs/Data/DESeq2.RData")

DEG_deseq2$genename <- rownames(DEG_deseq2)
DEG_deseq2_plot <- DEG_deseq2
DEG_deseq2_plot <- DEG_deseq2_plot[rownames(DEG_deseq2_plot) != c("C6orf58", "GAST"),]

library(ggVolcano)

pdf(file = "./Result/DFS_DEGs/Plot/DEG_DESeq2_Volcano Plot.pdf", width = 6, height = 6)
gradual_volcano(DEG_deseq2_plot, x = "log2FoldChange", y = "padj", pointSizeRange = c(0.5, 4),
                label = "genename",
                # custom_label = c("ERICH3","RPL35P5", "STYXL2", "CST4", "PAX3", "REG4", "MEG3", "NXPH4", "CEACAM20", "DPY19L2P1"),
                label_number = 10,
                fills = c("#0D0D4C", "#4063AF","#39bbec", "#7DCA99", "#FAED31","#f38466","#E51D6C", "#b81f25"),
                colors = c("#17194e", "#1E4D8B", "#1E9A78", "#49B8B0", "#F48B7C","#a22f27", "#211f1f"),
                log2FC_cut = 1, FDR_cut = 0.05, 
                output = FALSE)+ scale_x_continuous(limits = c(-6, 6))+ scale_y_continuous(limits = c(0, 15))
dev.off()

