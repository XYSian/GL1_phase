options(BioC_mirror = "https://bioconductor.org")
BiocManager::install(c("BiocParallel", "DESeq2"))
library(DESeq2)
install.packages("dplyr")
library(dplyr)
library(readxl)

read_and_extract_counts <- function(file_path, sample_name) {
  print(paste("File path is:", file_path))  # 输出当前文件路径以便调试
  df <- read.table(file_path, header=TRUE, sep="\t", stringsAsFactors=FALSE)
  counts <- df[, c(1, ncol(df))]
  colnames(counts) <- c("Geneid", sample_name) # 设置列名
  return(counts)
}
files <- list(
  "Red1" = "E:\\rstudio\\work\\ASE\\featurecount\\count0\\count74.txt",
  "Green3" = "E:\\rstudio\\work\\ASE\\featurecount\\count0\\count75.txt",
  "Green2" = "E:\\rstudio\\work\\ASE\\featurecount\\count0\\count76.txt",
  "Green1" = "E:\\rstudio\\work\\ASE\\featurecount\\count0\\count77.txt",
  "Yellow3" = "E:\\rstudio\\work\\ASE\\featurecount\\count0\\count78.txt",
  "Yellow2" = "E:\\rstudio\\work\\ASE\\featurecount\\count0\\count79.txt",
  "Yellow1" = "E:\\rstudio\\work\\ASE\\featurecount\\count0\\count80.txt",
  "Red3" = "E:\\rstudio\\work\\ASE\\featurecount\\count0\\count81.txt",
  "Red2" = "E:\\rstudio\\work\\ASE\\featurecount\\count0\\count82.txt"
)

count_data_list <- lapply(names(files), function(sample) {
  file_path <- files[[sample]] # 获取文件路径
  read_and_extract_counts(file_path, sample)
})
# 合并所有数据框，按基因ID对齐
count_data_combined <- Reduce(function(x, y) merge(x, y, by = "Geneid", all = TRUE), count_data_list)

# 保存合并后的计数数据
write.table(count_data_combined, file="count_data_combined.txt", sep="\t", row.names=FALSE, quote=FALSE)
# 定义一个函数来读取文件并提取基因ID和基因长度
extract_gene_lengths <- function(file_path) {
  df <- read.table(file_path, header=TRUE, sep="\t", stringsAsFactors=FALSE)
  # 提取基因ID和倒数第二列的基因长度
  gene_lengths <- df[, c(1, ncol(df)-1)]
  colnames(gene_lengths) <- c("Geneid", "Length")
  return(gene_lengths)
}

# 提取基因长度
gene_lengths_df <- extract_gene_lengths("E:\\rstudio\\work\\ASE\\featurecount\\count0\\count74.txt")

# 检查结果
head(gene_lengths_df)
# 将基因长度数据框转换为向量
gene_lengths <- gene_lengths_df$Length
names(gene_lengths) <- gene_lengths_df$Geneid
# 计算每个样本的总 Reads 数
total_counts <- colSums(count_data_combined[, -1])

# 计算每个基因的FPKM
calculate_fpkm <- function(counts, gene_length, total_counts) {
  fpkm <- (counts / (gene_length / 1000)) / (total_counts / 1e6)
  return(fpkm)
}

# 对每个样本应用FPKM计算
fpkm_matrix <- count_data_combined
for (i in 2:ncol(fpkm_matrix)) {
  sample_name <- colnames(fpkm_matrix)[i]
  fpkm_matrix[, i] <- calculate_fpkm(fpkm_matrix[, i], gene_lengths[fpkm_matrix$Geneid], total_counts[sample_name])
}

# 保存FPKM矩阵
write.table(fpkm_matrix, file="fpkm_matrix.txt", sep="\t", row.names=FALSE, quote=FALSE)

#########################################
# 更新文件路径
files <- list(
  "Red1" = "E:\\rstudio\\work\\ASE\\featurecount\\count1\\count74.txt",
  "Green3" = "E:\\rstudio\\work\\ASE\\featurecount\\count1\\count75.txt",
  "Green2" = "E:\\rstudio\\work\\ASE\\featurecount\\count1\\count76.txt",
  "Green1" = "E:\\rstudio\\work\\ASE\\featurecount\\count1\\count77.txt",
  "Yellow3" = "E:\\rstudio\\work\\ASE\\featurecount\\count1\\count78.txt",
  "Yellow2" = "E:\\rstudio\\work\\ASE\\featurecount\\count1\\count79.txt",
  "Yellow1" = "E:\\rstudio\\work\\ASE\\featurecount\\count1\\count80.txt",
  "Red3" = "E:\\rstudio\\work\\ASE\\featurecount\\count1\\count81.txt",
  "Red2" = "E:\\rstudio\\work\\ASE\\featurecount\\count1\\count82.txt"
)

# 读取并提取每个文件的计数数据
count_data_list1 <- lapply(names(files), function(sample) {
  file_path <- files[[sample]]
  read_and_extract_counts(file_path, sample)
})

# 合并所有数据框，按基因ID对齐，并将 NA 值替换为 0
count_data_combined1 <- Reduce(function(x, y) merge(x, y, by = "Geneid", all.x = TRUE, all.y = TRUE), count_data_list1)
count_data_combined1[is.na(count_data_combined1)] <- 0  # 替换 NA 值为 0

# 保存合并后的计数数据
write.table(count_data_combined1, file="count_data_combined1.txt", sep="\t", row.names=FALSE, quote=FALSE)

# 定义一个函数来读取文件并提取基因ID和基因长度
extract_gene_lengths1 <- function(file_path) {
  df <- read.table(file_path, header=TRUE, sep="\t", stringsAsFactors=FALSE)
  gene_lengths <- df[, c(1, ncol(df)-1)]
  colnames(gene_lengths) <- c("Geneid", "Length")
  return(gene_lengths)
}

# 提取基因长度
gene_lengths_df1 <- extract_gene_lengths1("E:\\rstudio\\work\\ASE\\featurecount\\count1\\count74.txt")

# 检查结果
head(gene_lengths_df1)

# 将基因长度数据框转换为向量
gene_lengths1 <- gene_lengths_df1$Length
names(gene_lengths1) <- gene_lengths_df1$Geneid

# 计算每个样本的总 Reads 数
total_counts1 <- colSums(count_data_combined1[, -1])

# 计算每个基因的FPKM
calculate_fpkm1 <- function(counts, gene_length, total_counts) {
  # 确保基因长度和总计数不为零
  fpkm <- ifelse(gene_length > 0 & total_counts > 0, 
                 (counts / (gene_length / 1000)) / (total_counts / 1e6), 
                 0)
  return(fpkm)
}

# 对每个样本应用FPKM计算
fpkm_matrix1 <- count_data_combined1
for (i in 2:ncol(fpkm_matrix1)) {
  sample_name <- colnames(fpkm_matrix1)[i]
  fpkm_matrix1[, i] <- calculate_fpkm1(fpkm_matrix1[, i], gene_lengths1[fpkm_matrix1$Geneid], total_counts1[sample_name])
}

# 保存FPKM矩阵
write.table(fpkm_matrix1, file="fpkm_matrix1.txt", sep="\t", row.names=FALSE, quote=FALSE)
############################################33
# 加载 readxl 包
library(readxl)

# 读取 Excel 文件
file_path <- "E:\\rstudio\\work\\ASE\\alle\\alle99-100.xlsx"
allele_data <- read_excel(file_path)

# 查看读取的数据
print(allele_data)

# 假设数据框有两列，分别是单倍型基因组1和单倍型基因组2
colnames(allele_data) <- c("Haplotype1", "Haplotype2")

# 检查数据框的前几行
head(allele_data)
# 假设 fpkm_matrix 是你的FPKM矩阵
fpkm_matrix$Geneid <- paste0(fpkm_matrix$Geneid, ".0")

# 检查结果
head(fpkm_matrix)
fpkm_matrix1$Geneid <- paste0(fpkm_matrix1$Geneid, ".1")

# 检查结果
head(fpkm_matrix1)

merged_fpkm <- merge(fpkm_matrix, fpkm_matrix1, by.x = "Geneid", by.y = "Geneid", all = TRUE)

# 对应关系表基因筛选并将基因ID重命名以进行匹配
hap1_exp <- merged_fpkm %>%
  filter(Geneid %in% allele_data$Haplotype1) %>%
  mutate(Haplotype2 = allele_data$Haplotype2[match(Geneid, allele_data$Haplotype1)])

hap2_exp <- merged_fpkm %>%
  filter(Geneid %in% allele_data$Haplotype2)

# 将单倍型1的表达数据与单倍型2的表达数据合并
final_data <- merge(hap1_exp, hap2_exp, by.x = "Haplotype2", by.y = "Geneid", suffixes = c("_Hap1", "_Hap2"))


# 处理合并后的数据，删除完全为空的列
final_data_clean <- final_data %>%
  select(where(~ !all(is.na(.))))

# 查看合并结果
print(final_data_clean)

# 重新命名 final_data_clean 的 Geneid 列为 Haplotype1
final_data_clean <- final_data_clean %>%
  rename(Haplotype1 = Geneid)
#######################################################################################
#################################################3333333333333333333333
############3333333333333333333333333333333333333
# 转换为长格式并添加颜色信息
library(tidyr)
library(dplyr)

# 将数据转换为长格式
long_data <- final_data_clean %>%
  pivot_longer(
    cols = starts_with(c("Red", "Green", "Yellow")),  # 根据实际列名调整
    names_to = "Sample_Haplotype",
    values_to = "FPKM"
  ) %>%
  # 提取颜色和单倍型信息
  mutate(
    Sample = sub("\\..*$", "", Sample_Haplotype),  # 提取颜色和样本编号
    Haplotype = ifelse(grepl("\\.x_", Sample_Haplotype), "hap1", "hap2"),
    Sample = sub("\\.x_|\\.y_", "", Sample),  # 去掉单倍型部分
    Sample = factor(Sample, levels = c("Red1", "Green1", "Yellow1", "Red2", "Green2", "Yellow2", "Red3", "Green3", "Yellow3"))
  ) %>%
  select(-Sample_Haplotype)  # 删除不需要的列

wide_data <- long_data %>%
  pivot_wider(
    names_from = c(Sample, Haplotype),
    values_from = FPKM
  )
head(wide_data)
library(tibble)

# 创建表达矩阵，使用 Haplotype1 作为行名
expression_matrix <- wide_data %>%
  select(-Haplotype2) %>%  # 删除 Haplotype2 列
  column_to_rownames(var = "Geneid") %>%  # 使用 Haplotype1 作为行名
  as.data.frame()

# 打印表达矩阵以确认其正确性
print(head(expression_matrix))
# 如果你确实有计数数据，确保数据为整数
expression_matrix <- round(expression_matrix)
# 保留每个基因在所有样本中 count 大于 6 的行
filtered_expression_matrix <- expression_matrix[rowSums(expression_matrix > 6) == ncol(expression_matrix), ]

# 查看过滤后的结果
print(head(filtered_expression_matrix))

# 创建样本条件数据框
sample_conditions <- data.frame(
  row.names = colnames(expression_matrix),  # 使用表达矩阵的列名
  Sample = gsub("_hap[12]$", "", colnames(expression_matrix)),
  Haplotype = gsub(".*_(hap[12])$", "\\1", colnames(expression_matrix))
)

# 创建 DESeq2 数据集
dds <- DESeqDataSetFromMatrix(
  countData = filtered_expression_matrix,
  colData = sample_conditions,
  design = ~ Sample + Haplotype
)

# 进行差异表达分析
dds <- DESeq(dds)
results <- results(dds)

# 查看结果
print(head(results))
######################################################################################
#########################################################################################
##############################################################################

# 为每种颜色提取对应的列，并分别构建表达矩阵
red_expression_matrix <- filtered_expression_matrix %>%
  select(starts_with("Red"))
green_expression_matrix <- filtered_expression_matrix %>%
  select(starts_with("Green"))

yellow_expression_matrix <- filtered_expression_matrix %>%
  select(starts_with("Yellow"))

# 提取对应颜色的样本信息
red_sample_conditions <- sample_conditions %>%
  filter(Sample %in% c("Red1", "Red2", "Red3"))

green_sample_conditions <- sample_conditions %>%
  filter(Sample %in% c("Green1", "Green2", "Green3"))

yellow_sample_conditions <- sample_conditions %>%
  filter(Sample %in% c("Yellow1", "Yellow2", "Yellow3"))

# 创建 DESeq2 数据集
dds_red <- DESeqDataSetFromMatrix(
  countData = round(red_expression_matrix),  # 将表达矩阵中的值四舍五入为整数
  colData = red_sample_conditions,
  design = ~ Haplotype
)

dds_green <- DESeqDataSetFromMatrix(
  countData = round(green_expression_matrix),
  colData = green_sample_conditions,
  design = ~ Haplotype
)

dds_yellow <- DESeqDataSetFromMatrix(
  countData = round(yellow_expression_matrix),
  colData = yellow_sample_conditions,
  design = ~ Haplotype
)

# 进行差异表达分析
dds_red <- DESeq(dds_red)
dds_green <- DESeq(dds_green)
dds_yellow <- DESeq(dds_yellow)

# 获取结果
results_red <- results(dds_red)
results_green <- results(dds_green)
results_yellow <- results(dds_yellow)
print(head(results_red))
library(dplyr)

# 设置筛选阈值
log2FC_threshold <- 1 # Log2 fold change threshold
padj_threshold <- 0.05  # Adjusted p-value (FDR) threshold

# 筛选差异表达基因
significant_genes_red <- results_red %>%
  as.data.frame() %>%
  filter(abs(log2FoldChange) >= log2FC_threshold & padj <= padj_threshold)

# 统计差异基因的数量
num_significant_genes_red <- nrow(significant_genes_red)
num_upregulated_red <- sum(significant_genes_red$log2FoldChange > 0)  # 上调基因数量
num_downregulated_red <- sum(significant_genes_red$log2FoldChange < 0)  # 下调基因数量

# 打印统计结果
cat("Total significant genes:", num_significant_genes_red, "\n")
cat("Upregulated genes:", num_upregulated_red, "\n")
cat("Downregulated genes:", num_downregulated_red, "\n")

# 查看差异基因列表
head(significant_genes_red)
# 统计差异基因的数量###########################grenn#########
significant_genes_green <- results_green %>%
  as.data.frame() %>%
  filter(abs(log2FoldChange) >= log2FC_threshold & padj <= padj_threshold)
num_significant_genes_green <- nrow(significant_genes_green)
num_upregulated_green <- sum(significant_genes_green$log2FoldChange > 0)  # 上调基因数量
num_downregulated_green <- sum(significant_genes_green$log2FoldChange < 0)  # 下调基因数量

# 打印统计结果
cat("Total significant genes:", num_significant_genes_green, "\n")
cat("Upregulated genes:", num_upregulated_green, "\n")
cat("Downregulated genes:", num_downregulated_green, "\n")

# 查看差异基因列表
head(significant_genes_green)

head(results_red)
# 统计差异基因的数量###########################YELLOW#########
significant_genes_yellow<- results_yellow %>%
  as.data.frame() %>%
  filter(abs(log2FoldChange) >= log2FC_threshold & padj <= padj_threshold)
num_significant_genes_yellow <- nrow(significant_genes_yellow)
num_upregulated_yellow <- sum(significant_genes_green$log2FoldChange > 0)  # 上调基因数量
num_downregulated_yellow <- sum(significant_genes_green$log2FoldChange < 0)  # 下调基因数量
cat("Total significant genes:", num_significant_genes_yellow, "\n")

####################匹配基因信息##############
# 读取 merged1.final 文件
merged_anno <- read.table("merged0_final", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
colnames(merged_anno) <- c("Geneid", "AnotherColumn", "GeneAnnotation")
#########打印
significant_genes_red$RowNames <- rownames(significant_genes_red)
# Step 2: 使用新创建的 RowNames 列进行合并
matched_results_red<- merge(significant_genes_red, merged_anno, by.x = "RowNames", by.y = "Geneid", all.x = TRUE)
# 从 wide_data 中提取 Geneid 对应的 Haplotype2 列
matched_results_red<- merge(matched_results_red, wide_data[, c("Geneid", "Haplotype2")], 
                                by.x = "RowNames", by.y = "Geneid", all.x = TRUE)

# 查看结果
head(matched_results_red)
write.csv(matched_results_red, file = "matched_results_red_fd2.csv", quote = FALSE, row.names = FALSE)
####grenn####
# Step 1: 将行名转换为一列
significant_genes_green$RowNames <- rownames(significant_genes_green)
# Step 2: 使用新创建的 RowNames 列进行合并
matched_results_green <- merge(significant_genes_green, merged_anno, by.x = "RowNames", by.y = "Geneid", all.x = TRUE)
# 从 wide_data 中提取 Geneid 对应的 Haplotype2 列
matched_results_green <- merge(matched_results_green, wide_data[, c("Geneid", "Haplotype2")], 
                             by.x = "RowNames", by.y = "Geneid", all.x = TRUE)

# 查看结果
head(matched_results_green)
write.csv(matched_results_green, file = "matched_results_greenfd2.csv", quote = FALSE, row.names = FALSE)
#######yellow###
significant_genes_yellow$RowNames <- rownames(significant_genes_yellow)
# Step 2: 使用新创建的 RowNames 列进行合并
matched_results_yellow <- merge(significant_genes_yellow, merged_anno, by.x = "RowNames", by.y = "Geneid", all.x = TRUE)
# 从 wide_data 中提取 Geneid 对应的 Haplotype2 列
matched_results_yellow <- merge(matched_results_yellow, wide_data[, c("Geneid", "Haplotype2")], 
                             by.x = "RowNames", by.y = "Geneid", all.x = TRUE)

# 查看结果
head(matched_results_yellow)
write.csv(matched_results_yellow, file = "matched_results_yellow.csv", quote = FALSE, row.names = FALSE)
write.csv(expression_matrix, file = "expression_matrix.csv", quote = FALSE, row.names = TRUE)


########################################火山图
i# 设置 Bioconductor 版本
# 首先安装 devtools 包（如果尚未安装）
install.packages("devtools")

# 从 GitHub 安装 EnhancedVolcano
devtools::install_github("kevinblighe/EnhancedVolcano")
BiocManager::install("EnhancedVolcano")
library(EnhancedVolcano)
library(EnhancedVolcano)
significant_genes_red$Geneid <- rownames(significant_genes_red)

# 提取分析结果
results_red <- results(dds_red)

# 将结果转换为数据框格式
results_red_df <- as.data.frame(results_red)

# 绘制火山图
EnhancedVolcano(results_red_df,
                lab = rownames(results_red_df),
                x = 'log2FoldChange',
                y = 'padj',
                pCutoff = 0.01,
                FCcutoff = 2,
                pointSize = 3.0,
                labSize = 3.0,
                colAlpha = 0.5,
                title = "Volcano Plot",
                subtitle = "Phase0 vs Phase1 DEA in Red Leaf")
# 提取分析结果
results_red <- results(dds_red)

# 将结果转换为数据框格式
results_red_df <- as.data.frame(results_red)

# 绘制火山图
EnhancedVolcano(results_red_df,
                lab = rownames(results_red_df),
                x = 'log2FoldChange',
                y = 'padj',
                #selectLab =c('Aco_HBLgroup10g016157.0'),
                pCutoff = 0.01,
                FCcutoff = 2,
                pointSize = 3.0,
                labSize = 3.0,
                colAlpha = 0.5,
                labCol ='black',
                labFace = 'bold',
                #boxedLabels = TRUE,
                legendPosition ='right',
                drawConnectors = TRUE,
                title = "Volcano Plot",
                subtitle = "Phase0 vs Phase1 DEA in Red Tissue")
# 火山图-------绿色
results_green <- results(dds_green)

# 将结果转换为数据框格式
results_green_df <- as.data.frame(results_green)

# 绘制火山图
EnhancedVolcano(results_green_df,
                lab = rownames(results_green_df),
                x = 'log2FoldChange',
                y = 'padj',
                pCutoff = 0.01,
                FCcutoff = 2,
                pointSize = 3.0,
                labSize = 3.0,
                colAlpha = 0.5,
                title = "Volcano Plot",
                subtitle = "Phase0 vs Phase1 DEA in Green Leaf")
# 火山图-------黄色
results_yellow<- results(dds_yellow)

# 将结果转换为数据框格式
results_yellow_df <- as.data.frame(results_yellow)

# 绘制火山图
EnhancedVolcano(results_yellow_df,
                lab = rownames(results_yellow_df),
                x = 'log2FoldChange',
                y = 'padj',
                pCutoff = 0.01,
                FCcutoff = 2,
                pointSize = 3.0,
                labSize = 3.0,
                colAlpha = 0.5,
                title = "Volcano Plot",
                subtitle = "Phase0 vs Phase1 DEA in Yellow Leaf")

