gene_for_survival <- c(intersect(DEG_edgeR_down_gene, DEG_deseq2_down_gene), intersect(DEG_edgeR_up_gene, DEG_deseq2_up_gene))
df_gene_for_survival <- tcga_paad_exp[gene_for_survival, ]
df_gene_for_survival <- as.data.frame(t(df_gene_for_survival))

library(openxlsx)
tcga_survival <- read.xlsx("./Analysis/Data/TCGA-PAAD_pathomics_survival_data.xlsx")
tcga_survival <- tcga_survival[,1:5]

df_gene_for_survival$id <- rownames(df_gene_for_survival)
df_gene_for_survival <- merge(tcga_survival, df_gene_for_survival, by.x = "ID", by.y = "id")
DFS_DEGs_for_survival <- df_gene_for_survival
save(DFS_DEGs_for_survival, file = "./Result/DFS_DEGs/Data/DFS_DEGs_gene_for_survival.RData")
#---------------------------------------------------------------
##-------------------#### Survival Analysis ####-----------------
##---------------- # Based on the best cutpoint -----------------

# ---------------- ### OS Plot ### -------------------

## Survival Analysis --
library(survminer)
library(ggplot2)
library(ggsurvfit)
library(survival)
library(openxlsx)

# data <- read.xlsx("./Data/Targetedgene_summary_patient.xlsx")
data <- df_gene_for_survival
colnames(data)

variables <- colnames(data)[6:18]
variables

#### -------------- OS Plot --------------- ####
for (var in variables){
  rt = data[, c("OS.time", "OS", var)]
  rt$OS.time <- as.numeric(rt$OS.time)
  res.cut <- surv_cutpoint(rt, time = "OS.time", event = "OS", variables = var) # Calculate optimal cutpoint
  res.cat <- surv_categorize(res.cut) # Categorize based on cutpoint
  res.cat$group <- factor(res.cat[[var]], labels = c(paste0(var,"_High"),paste0(var,"_Low")))
  fit <- survfit2(Surv(OS.time, OS) ~ group, data = res.cat)
  diff <- survdiff(Surv(OS.time, OS) ~ group, data = res.cat)
  pValue <- 1 - pchisq(diff$chisq, df = 1)
  if (pValue < 0.001) {
    pValue = "p<0.001"
  } else {
    pValue = paste0("p=", sprintf("%.3f", pValue))
  }
  
  KMplot_best_cutpoint <- ggsurvfit(fit) +  
    add_risktable(risktable_height = 0.25, 
                  risktable_stats = c("{n.risk} ({cum.censor})"),
                  stats_label = list(n.risk = "Number at risk", cum.censor = "number censored"),
                  size = 4.5,
                  hjust = 0.5,
                  theme =   # Increase font size for risk table title and y-axis labels
                    list(
                      theme_risktable_default(axis.text.y.size = 14,
                                              plot.title.size = 14),
                      theme(plot.title = element_text(face = "bold")))  
    ) +  
    add_risktable_strata_symbol(symbol = "\U25CF", size = 20) +
    add_censor_mark(size = 4, shape = 73) +
    add_quantile(y_value = 0.5, linetype = "dashed", color = "black", linewidth = 0.6) +
    labs(title = "",
         x = "Time since surgery (days)",
         y = "OS (%)") +
    scale_x_continuous(breaks = seq(0, 2000, 500), expand = c(0.04, 0)) + 
    scale_y_continuous(breaks = seq(0, 1, 0.25), labels = seq(0, 100, 25), expand = c(0, 0)) +  
    scale_color_manual(values = c('#E05133', '#4877B9')) +  
    scale_fill_manual(values = c('#E05133', '#4877B9')) +  
    guides(color = guide_legend(ncol = 1)) +  
    theme_classic() +  
    theme(axis.text = element_text(size = 14, color = "black", family = "Calibri"),
          axis.ticks.length = unit(2, "mm"),
          axis.title.y = element_text(size = 14, color = "black", family = "Calibri"),
          axis.title.x = element_text(size = 14, color = "black", family = "Calibri"),
          panel.grid = element_blank(),
          legend.text = element_text(size = 13, color = "black", family = "Calibri"),
          legend.background = element_blank(),
          legend.position.inside = c(0.12, 0.12))+
    annotate("text", x = 1600, y = 0.05, label = paste(pValue), size = 5, color = "black", , fontface = "bold", family = "Calibri")
  filename <- paste0("./Result/OS Plot/DEGs/OS_KMplot_best_cutpoint_",var,".png")
  # Save plot to file
  ggsave(filename, plot = KMplot_best_cutpoint, width = 10, height = 7, dpi = 600)
}

