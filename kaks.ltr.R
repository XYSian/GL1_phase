
getwd()
setwd("E:\\rstudio\\work\\ASE\\ltr")
all_genes <- read.table("all.extracted_genes.txt", header = FALSE, stringsAsFactors = FALSE)
LTR_gene1_closest_within_2kb <- read.table("LTR_gene1_closest_within_2kb.bed", header = FALSE, stringsAsFactors = FALSE)
LTR_gene1_closest_within_2kb$V8 <- paste0(LTR_gene1_closest_within_2kb$V8, ".1")
print(LTR_gene1_closest_within_2kb)


head(merged_anno1)
print(LTR_gene0_overlap)
LTR_gene1_closest_within_2kb <- merge(LTR_gene1_closest_within_2kb, merged_anno1[, c("Haplotype1", "GeneAnnotation")], 
                           by.x = "V8", by.y = "Haplotype1", all.x = TRUE)
write.csv(as.data.frame(LTR_gene1_closest_within_2kb),"LTR_gene1_closest_within_2kb_anno.csv",row.names = F)

library(clusterProfiler)
library(dplyr)
library(ggplot2)
#install.packages("viridis")
library(viridis)
# 进行富集分析
x <- enricher(gene = LTR_exon1_closest$V1, 
              TERM2GENE = extracted_columns, 
              TERM2NAME = bg_terms_merged, 
              pvalueCutoff = 0.05, 
              #universe = gene,
              pAdjustMethod = "BH"
             )

# 检查结果
print(x)
write.csv(as.data.frame(x),"LTR_gene0_overlap.go.csv",row.names = F)
# 绘制 dotplot 并提取数据
p <- dotplot(x, showCategory=20)

# 修改图形对象以添加颜色渐变
p + scale_color_gradient(low = "blue", high = "red") +
  theme_minimal() + 
  ggtitle("GO Enrichment Analysis of Overlap Genes") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))



###############################
z <- enricher(
  gene = LTR_gene0_overlap$V8, 
  TERM2GENE = koid_to_gene, 
  TERM2NAME = koid2name,
  pvalueCutoff = 0.05, 
   pAdjustMethod = "BH"

)

# 查看结果
print(z)
barplot(z, showCategory = 20)+ ggtitle("KEGG Pathway Enrichment Analysis Of Overlap Genes")

# 绘制条形图
write.csv(as.data.frame(z),"LTR_gene0_overlap.kegg.csv",row.names = F)
# 绘制条形图
ggplot(top10 , aes(x = Description, y = Count, fill = p.adjust)) +
  geom_bar(stat = "identity") +
  coord_flip() +  # 让条形图水平显示
  scale_fill_gradient(low = "red", high = "skyblue") +  # 设置颜色梯度（红色为显著）
  labs(
    x = "Pathway", 
    y = "Count", 
    title = "KEGG Analysis : Red vs Green",
    fill = "pvalue"
  ) +
  theme_minimal()
#############################
setwd("E:\\rstudio\\work\\ASE\\varient")

indel.exonic_variant <- read.table("indel.exonic_variant_function", header = FALSE, sep = "\t", stringsAsFactors = FALSE)
head(indel.exonic_variant)



# 加载所需的包
library(ggplot2)
library(dplyr)
# 加载所需的包
library(ggplot2)

# 创建数据框，使用 V4 列作为染色体信息
mutation_data <- data.frame(
  mutation_types = indel.exonic_variant$V2,
  positions = indel.exonic_variant$V6,  # 基因组位置
  chromosome = indel.exonic_variant$V4  # 染色体信息
)

# 查看数据框头部确认
head(mutation_data)
# 加载所需的包
library(ggplot2)

# 创建数据框，使用 V4 列作为染色体信息
mutation_data <- data.frame(
  mutation_types = indel.exonic_variant$V2,
  positions = indel.exonic_variant$V6,  # 基因组位置
  chromosome = indel.exonic_variant$V4  # 染色体信息
)

# 查看数据框头部确认
head(mutation_data)
# 将基因组位置转换为Mb（百万碱基对）
mutation_data$positions_mb <- mutation_data$positions / 1e6
# 确保按1-25顺序排列chromosome
# 去掉chromosome列中的"group"前缀
mutation_data$chromosome <- gsub("group", "", mutation_data$chromosome)

# 确保按1-25顺序排列chromosome
mutation_data$chromosome <- factor(mutation_data$chromosome, levels = as.character(25:1))

# 绘制散点图，并美化图形，调整网格线
ggplot(mutation_data, aes(y = chromosome, x = positions_mb, color = mutation_types)) +
  geom_jitter(alpha = 0.6, height = 0.2) +  # 添加轻微的抖动，避免点重叠，调整为height
  scale_color_brewer(palette = "Set1") +  # 使用颜色调色板
  labs(
    title = "Genomic Position Distribution by Mutation Type",
    x = "Genomic Position (Mb)",  # 更改X轴标签为Mb单位
    y = "Chromosome"  # 更改Y轴标签为染色体
  ) +
  theme_minimal() +  # 简洁主题
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),  # 中央对齐标题，增大字体，粗体
    axis.title = element_text(size = 14, face = "bold"),  # 设置轴标题字体大小和粗体
    axis.text = element_text(size = 12)  # 设置轴标签字体大小
  ) +
  theme(panel.grid.major.x = element_blank(),  # 隐藏X轴的主网格线
        panel.grid.minor.x = element_blank(),  # 隐藏X轴的次网格线
        panel.grid.major.y = element_line(color = "gray", size = 0.5),  # 显示Y轴的主网格线
        panel.grid.minor.y = element_line(color = "lightgray", size = 0.25))  # 显示Y轴的次网格线



# 计算每个染色体和突变类型的数量
mutation_counts_by_type <- table(mutation_data$chromosome, mutation_data$mutation_types)

# 绘制堆积柱形图
ggplot(as.data.frame(mutation_counts_by_type), aes(x = Var1, y = Freq, fill = Var2)) +
  geom_bar(stat = "identity") +  # 堆积柱形图
  scale_fill_brewer(palette = "Set3") +  # 使用颜色调色板
  labs(title = "Mutation Type Distribution by Chromosome", 
       x = "Chromosome", y = "Mutation Count") +
  theme_minimal() +  # 简洁主题
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # 调整X轴标签角度
# 使用 ggplot2 绘制热图
ggplot(as.data.frame(mutation_counts_by_type), aes(x = Var1, y = Var2, fill = Freq)) +
  geom_tile() +  # 绘制热图
  scale_fill_gradient(low = "white", high = "red") +  # 设置颜色梯度
  labs(title = "Heatmap of Mutation Counts by Chromosome and Mutation Type", 
       x = "Chromosome", y = "Mutation Type") +
  theme_minimal() +  # 简洁主题
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # 调整X轴标签角度

##########################################
