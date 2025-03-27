# 安装 clusterProfiler
BiocManager::install("clusterProfiler")
# 加载 clusterProfiler 包
library(clusterProfiler)
# 读取 phase0.pep.fasta.tsv 文件
data <- read.delim("phase0.pep.fasta.tsv", header = FALSE, sep = "\t")
# 从数据中提取基因ID和GO术语
gene_go_df <- data %>%
  select(GeneID, GO)

# 提取GO术语列表
bg_terms <- unique(gene_go_df$GO)

head(bg_terms)

bg_split <- tibble(original = bg_terms) %>%
  separate_rows(original, sep = "\\|") %>%
  filter(original != "-" & original != "")

# 查看拆分后的结果
head(bg_terms)
head(bg_split)
library(GO.db)
go_data <- select(GO.db, keys(GO.db, keytype = "GOID"), columns = c("TERM", "DEFINITION"), keytype = "GOID")

# 查看数据
head(go_data)
# 合并GO编号与GO名称和描述
bg_terms_merged <- merge(bg_split, go_data, by.x = "original", by.y = "GOID", all.x = TRUE)

library(dplyr)
library(tidyr)

# 拆分GO term，并将数据展平
gene_go_df_cleaned <- gene_go_df %>%
  filter(GO != "-") %>% # 移除GO列为空或为"-"的行
  separate_rows(GO, sep = "\\|") %>% # 按"|"分隔GO term
  drop_na() # 移除任何剩余的NA值

# 检查数据
head(term2name)
library(clusterProfiler)
library(dplyr)

extracted_columns <- gene_go_df_cleaned %>%
  dplyr::select(GO, GeneID)


# 查看提取后的结果
print(extracted_columns)
head(gene)
head(extracted_columns)
head(term2name)
significant_genes_red$Geneid <- rownames(significant_genes_red)
significant_genes_yellow$Geneid <- rownames(significant_genes_yellow)

gene_YELLOW <- as.vector(significant_genes_yellow$Geneid)

gene_red <- as.vector(significant_genes_red$Geneid)
significant_genes_green$Geneid <- rownames(significant_genes_green)

gene_green <- as.vector(significant_genes_green$Geneid)


res_Red_vs_Green_sig_new$Geneid <- rownames(res_Red_vs_Green_sig_new)
Red_vs_Green  <- as.vector(significant_genes_red$Geneid)

res_Red_vs_Yellow_sig_new$Geneid <- rownames(res_Red_vs_Yellow_sig_new)
Red_vs_Yellow  <- as.vector(res_Red_vs_Yellow_sig_new$Geneid)

res_Green_vs_Yellow_sig_new$Geneid <- rownames(res_Green_vs_Yellow_sig_new)
Green_vs_Yellow  <- as.vector(res_Green_vs_Yellow_sig_new$Geneid)


# 进行富集分析
x <- enricher(gene = Green_vs_Yellow, 
              TERM2GENE = extracted_columns, 
              TERM2NAME = bg_terms_merged, 
              pvalueCutoff = 0.05, 
              #universe = gene,
              pAdjustMethod = "BH"
             )

# 检查结果
print(x)
write.csv(as.data.frame(x),"Green_vs_Yellow.GO.csv",row.names = F)
