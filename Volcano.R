########################火山图，热图
BiocManager::install('C')
install.packages("E:/NEWRNA/RNASEQ/testrna/EnhancedVolcano_1.18.0.zip",repos=NULL)
library('EnhancedVolcano')
library(dplyr)
EnhancedVolcano(de_result,
                lab = de_result$geneid,
                x ='logFC',
                y = 'PValue',
                title = 'vocano_plot',
                subtitle = 'gr_vs_gr',
                xlim = c(-10,10),
                pCutoff = 0.05,
                FCcutoff = 1)
library(dplyr)
library(tibble)
library(tidyverse)
top_de_exp <- dplyr::slice(de_result, 1:20) %>%
  dplyr::select(geneid, gr, re ) %>%
column_to_rownames(var = 'geneid')
library(pheatmap)
sample_info <- read.table("E:/NEWRNA/RNASEQ/testrna/sample.txt",header = T , row.names=1)
pheatmap(top_de_exp,
         scale = 'row',
         color = colorRampPalette(c("green","white","red"))(200),
         annotation_col = dplyr::select(sample_info,color),
         annotation_colors = list(
           color =c(green = '#4DBBD5FF',
                    red = '#E27B25FF')),
         cutree_rows = 2)
