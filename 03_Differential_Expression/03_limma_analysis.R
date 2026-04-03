########################### Differential analysis using limma ##################

library(limma)
exprSet <- tcga_paad_exp
design <- model.matrix(~0 + factor(group))
colnames(design) = levels(factor(group))
rownames(design) = colnames(exprSet)

dge <- DGEList(counts = exprSet, group = group)

keep <- filterByExpr(dge)
dge <- dge[keep, , keep.lib.sizes = FALSE]

dge <- calcNormFactors(dge)

v <- voom(dge, design, plot = TRUE, normalize = "quantile")

fit <- lmFit(v, design)

con <- paste(rev(levels(group)), collapse = "-")
con
# [1] "High_risk-Low_risk"

cont.matrix <- makeContrasts(contrasts = c(con), levels = design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2)

# Get differential expression results
tempOutput <- topTable(fit2, coef = con, n = Inf)
DEG_limma_voom <- na.omit(tempOutput)
head(DEG_limma_voom)
#               logFC  AveExpr         t      P.Value adj.P.Val         B
# DHRS9     1.2527093 4.058359  3.378189 0.0008968156 0.9999236 -1.765191
# C22orf46 -0.2069325 4.263962 -3.287879 0.0012169993 0.9999236 -1.907436
# SESN2     0.2505554 4.546049  3.195374 0.0016530302 0.9999236 -1.988507
# BBS2     -0.1581822 5.064252 -3.039053 0.0027319733 0.9999236 -2.203175
# KRCC1    -0.1784422 5.158943 -2.946298 0.0036477766 0.9999236 -2.362523
# EEF1A1P5  0.4729866 3.863443  3.025150 0.0028541809 0.9999236 -2.449084

logFC = 1
P.Value = 0.05
k1 <- (DEG_limma_voom$P.Value < P.Value) & (DEG_limma_voom$logFC < -logFC)
k2 <- (DEG_limma_voom$P.Value < P.Value) & (DEG_limma_voom$logFC > logFC)
DEG_limma_voom <- mutate(DEG_limma_voom, change = ifelse(k1, "down", ifelse(k2, "up", "stable")))
table(DEG_limma_voom$change)

# Save differential expression results to RData file
save(DEG_limma_voom, file = './Result/DEG_limma_voom.Rdata')

