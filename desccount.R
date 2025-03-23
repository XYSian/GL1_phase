# 加载必要的库
library(ggplot2)
options(repos='http://cran.rstudio.com/')
install.packages("ComplexHeatmap")
library(pheatmap)
library(dplyr)

# 加载表达矩阵和 ASE 基因对
expression_matrix <- read.csv("expression_matrix.csv", row.names = 1, check.names = FALSE)
ase_genes <- read.table("log1.kaks_genes.txt", header = FALSE, stringsAsFactors = FALSE)
ase_kaks <- read.table("high_kaks_genes.txt", header = FALSE, stringsAsFactors = FALSE)
gene_data <- read.table("merged0_final", header = FALSE, sep = " ", stringsAsFactors = FALSE, fill = TRUE)
print(ase_kaks)
library(dplyr)
library(tidyr)
# 选择V1和V5列
ase_kaks <- ase_kaks[, c("V1", "V5")]

# 去除V5列中包含NA的行
ase_kaks <- ase_kaks[!is.na(ase_kaks$V5), ]

# 查看清理后的数据框
head(ase_kaks)

# 合并第 3 到第 14 列
gene_data <- gene_data %>%
  mutate(Description = apply(select(., 3:13), 1, function(row) paste(na.omit(row), collapse = " "))) %>%
  select(V1, V2, Description) %>%
  rename(GeneID = V1, ProteinID = V2)

# 查看结果
head(gene_data)
# 给列命名
colnames(gene_data) <- c("GeneID", "ProteinID", "Description")
colnames(ase_genes) <- c("GeneID", "kaks")

# 检查前几行
head(gene_data)

head(expression_matrix)
# 假设 ase_kaks 数据框已加载
# 假设 ase_kaks 数据框已加载
ase_kaks$GeneID <- sub("-.*", ".0", ase_kaks$GeneID)

# 查看结果
head(ase_kaks)
ase_genes$GeneID <- sub("-.*", ".0", ase_genes$GeneID)

# 查看结果
head(ase_kaks)

head(ase_genes)
head(ase_kaks)

# 提取匹配基因的表达矩阵行
matched_rows <- expression_matrix[rownames(expression_matrix) %in% ase_kaks$GeneID, ]
# 提取包含 "hap1" 的列
hap1_columns <- grep("hap1", colnames(matched_rows), value = TRUE)
hap1_expression <- matched_rows[, hap1_columns, drop = FALSE]
head(hap1_expression)
write.csv(hap1_expression, "hap1_expression.csv", row.names = TRUE)

hap2_columns <- grep("hap2", colnames(matched_rows), value = TRUE)
hap2_expression <- matched_rows[, hap2_columns, drop = FALSE]
head(hap2_expression)
write.csv(hap2_expression, "hap2_expression.csv", row.names = TRUE)
ase_kaks$V1 <- sub("-.*", ".0", ase_kaks$V1)
# 查看修改后的数据

# 加载必要包
library(ggplot2)
library(ComplexHeatmap)
library(gridExtra)
BiocManager::install("ComplexHeatmap")
# 1. 准备数据
# 提取 Ka/Ks 比率数据
kaks_data <- ase_kaks[-1, ] # 移除第一行标题行
colnames(kaks_data) <- c("Gene", "Method", "Ka", "Ks", "Ka_Ks")
kaks_data$kaks <- as.numeric(kaks_data$kaks)

# 按 Gene 排序
kaks_data <- kaks_data[order(kaks_data$Gene), ]
print(hap1_expression)
# 提取红色样本的表达量列
hap1_red <- hap1_expression[, grepl("Red", colnames(hap1_expression))]
hap2_red <- hap2_expression[, grepl("Red", colnames(hap2_expression))]
head(kaks_data)
# 提取绿色样本的表达量列
hap1_green <- hap1_expression[, grepl("Green", colnames(hap1_expression))]
hap2_green <- hap2_expression[, grepl("Green", colnames(hap2_expression))]

# 提取黄色样本的表达量列
hap1_yellow <- hap1_expression[, grepl("Yellow", colnames(hap1_expression))]
hap2_yellow <- hap2_expression[, grepl("Yellow", colnames(hap2_expression))]

# 计算每种颜色的平均表达量
hap1_red_avg <- rowMeans(hap1_red)
hap1_green_avg <- rowMeans(hap1_green)
hap1_yellow_avg <- rowMeans(hap1_yellow)

hap2_red_avg <- rowMeans(hap2_red)
hap2_green_avg <- rowMeans(hap2_green)
hap2_yellow_avg <- rowMeans(hap2_yellow)
# 合并表达量数据和 Ka/Ks 数据
# 通过 V1 或 V2 列与表达量数据对齐
colnames(hap2_expression) <- c("GeneID", "kaks")

hap1_expression_ordered <- hap1_expression[rownames(hap1_expression) %in% ase_genes$GeneID, ]
hap2_expression_ordered <- hap2_expression[rownames(hap2_expression) %in% ase_genes$GeneID, ]

# 对 Ka/Ks 数据按基因排序
kaks_data_ordered <- kaks_data[order(kaks_data$Gene), ]

# 提取相关基因的表达量
expression_avg <- data.frame(
  Gene = rownames(hap1_expression_ordered),
  Hap1_Red = rowMeans(hap1_red[rownames(hap1_expression_ordered), ]),
  Hap1_Green = rowMeans(hap1_green[rownames(hap1_expression_ordered), ]),
  Hap1_Yellow = rowMeans(hap1_yellow[rownames(hap1_expression_ordered), ]),
  Hap2_Red = rowMeans(hap2_red[rownames(hap2_expression_ordered), ]),
  Hap2_Green = rowMeans(hap2_green[rownames(hap2_expression_ordered), ]),
  Hap2_Yellow = rowMeans(hap2_yellow[rownames(hap2_expression_ordered), ])
)

# 整合到一个数据框
plot_data <- merge(kaks_data_ordered, expression_avg, by.x = "GeneID", by.y = "Gene")
# 1. 过滤 Ka/Ks 小于 10 的基因
filtered_plot_data <- plot_data[plot_data$kaks <= 10, ]

# 2. 提取需要标准化的表达量数据列
expression_data <- filtered_plot_data[, c("Hap1_Red", "Hap1_Green", "Hap1_Yellow", "Hap2_Red", "Hap2_Green", "Hap2_Yellow")]

# 3. 对表达量数据进行 Z-score 标准化
standardized_expression_data <- t(scale(t(expression_data)))

# 4. 将标准化后的表达量数据与基因 ID 和 Ka/Ks 值结合
standardized_plot_data <- cbind(filtered_plot_data[, c("GeneID", "kaks")], standardized_expression_data)

# 查看标准化后的数据
head(standardized_plot_data)

library(ggplot2)
install.packages("circlize")
library(circlize)
install.packages("plotly")
library(gtable)

head(plot_data)
head(gene_data)

library(ComplexHeatmap)
###############################################
# 准备热图数据
heatmap_data <- as.matrix(standardized_plot_data[, c("Hap1_Red", "Hap1_Green", "Hap1_Yellow", 
                                        "Hap2_Red", "Hap2_Green", "Hap2_Yellow")])

# 确保基因描述和数据匹配
gene_names <- gene_data$Description[match(standardized_plot_data$Gene, gene_data$GeneID)]
head(gene_names)
# 用基因描述替换行名

# 准备 Ka/Ks 比率注释
annotation_data <- data.frame(Ka_Ks = standardized_plot_data$kaks)
rownames(annotation_data) <- rownames(heatmap_data)  # 确保行名与热图一致
annotation_data <- edit(annotation_data)
write.csv(heatmap_data, "heatmap_data.csv", row.names = TRUE)
heatmap_data <- read.csv("heatmap_data(1).csv", row.names = 1, check.names = FALSE)
head(heatmap_data)
head(annotation_data)

# 绘制热图
library(pheatmap)

# 假设 heatmap_data 和 annotation_data 已经定义好
library(pheatmap)
library(ggplot2)

# 将数据框转换为矩阵
heatmap_data_matrix <- as.matrix(heatmap_data)

# 生成热图对象
heatmap_obj <- pheatmap(
  heatmap_data_matrix,
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  display_numbers = FALSE,
  color = colorRampPalette(c("white", "red"))(100),
  main = "Expression Heatmap with Ka/Ks",
  annotation_row = annotation_data,  # 只显示 Ka/Ks 注释
  annotation_colors = list(
    Ka_Ks = colorRampPalette(c("pink", "red"))(100)
  ),  # 定义 Ka/Ks 注释颜色
  fontsize_row = 10,  # 行名称字体大小，减小字体
  fontsize_col = 10   # 列名称字体大小，减小字体
)


# 使用ggsave保存热图并调整尺寸
ggsave(
  "heatmap.png", 
  plot = heatmap_obj$gtable,  # 获取热图的gtable对象
  width = 6,  # 调整宽度（单位为英寸）
  height = 12  # 调整高度（单位为英寸）
)

###############################################
library(corrplot)

# 转置数据，使得行表示基因，列表示样本
cor_matrix <- cor(t(heatmap_data))

# 绘制基因相关性热图
corrplot(cor_matrix, 
         method = "color", 
         col = colorRampPalette(c("white", "red"))(100),  # 使用渐变色
         type = "full",  # 显示完整的相关性矩阵
      #   addCoef.col = "black",  # 在热图中显示相关性系数
         tl.col = "black",  # 设置标签颜色
         tl.srt = 45,  # 设置标签旋转角度
         title = "Gene Expression Correlation Heatmap",  # 设置标题
         mar = c(5, 5, 4, 2),  # 设置边距
         tl.cex = 0.6,  # 设置标签字体大小
         cl.lim = c(-1, 1),  # 设置色条的范围
         cl.cex = 0.6,  # 设置色条字体大小
         number.cex = 0.7,  # 设置相关系数数字字体大小
         
       diag = FALSE  # 去掉对角线的数值（自相关）
)

###################################################
all_kaks <- read.table("all.txt", header = FALSE, stringsAsFactors = FALSE)
head(all_kaks)
all_kaks <- all_kaks %>%
  filter(V1 != "Sequence") %>%  # 去除标题行
  mutate(Chromosome = gsub("Aco_HBLgroup(\\d+)g.*", "\\1", V1),  # 提取染色体编号
         Ka_Ks = as.numeric(V5))  # 确保 Ka/Ks 是数值类型
library(gtable)
library(dplyr)
library(ggplot2)
library(viridis)  # 提供美观配色方案
library(ggplot2)
library(viridis)  # 用于配色
library(devtools)
# 创建箱线图
# 对 Ka/Ks 比率进行 log10 转换
all_kaks$logKa_Ks <- log10(all_kaks$Ka_Ks)

# 过滤掉 logKa_Ks 范围之外的数据
all_kaks_filtered <- all_kaks[all_kaks$logKa_Ks >= -2 & all_kaks$logKa_Ks <= 2, ]
# 确保 Chromosome 列为因子并按顺序排列
all_kaks_filtered$Chromosome <- factor(all_kaks_filtered$Chromosome, 
                                       levels = as.character(1:25))  # 设置染色体的顺序为 1 到 25

# 创建箱线图
boxplot(logKa_Ks ~ factor(Chromosome), data = all_kaks_filtered, 
        col = rainbow(length(unique(all_kaks_filtered$Chromosome))),
        main = "Log10(Ka/Ks) Ratio Distribution by Chromosome",
        xlab = "Chromosome", ylab = "Log10(Ka/Ks) Ratio",
        las = 2)  # 旋转 x 轴标签
# 加载RColorBrewer包
library(RColorBrewer)

# 创建箱线图，使用马卡龙色调色板
boxplot(logKa_Ks ~ factor(Chromosome), data = all_kaks_filtered, 
        col = brewer.pal(length(unique(all_kaks_filtered$Chromosome)), "Set3"),  # 使用Set3色盘
        main = "Log10(Ka/Ks) Ratio Distribution by Chromosome",
        xlab = "Chromosome", ylab = "Log10(Ka/Ks) Ratio",
        las = 2)  # 旋转 x 轴标签
#################################################################3
#################################################################33
#######################ASE 染色体分布##################
alle <- read.table("matched_genes_location.txt", header = FALSE, stringsAsFactors = FALSE)
ase <- read.table("matched_location.txt", header = FALSE, stringsAsFactors = FALSE)
head(alle)
library(ggplot2)
library(dplyr)
colnames(ase) <- c("Chromosome", "Start", "End", "GeneID")

# 计算每个基因的大小（结束位置 - 起始位置）
genes <- ase %>%
  mutate(Gene_Length = End - Start + 1)

# 将起始位置转换为Mb (百万碱基对)
genes <- genes %>%
  mutate(Start_MB = Start / 1e6, End_MB = End / 1e6)

# 去掉染色体前缀 'group'，仅保留数字部分
genes$Chromosome <- gsub("group", "", genes$Chromosome)

# 将染色体列转换为有序因子，以确保按照1-25的顺序排列
genes$Chromosome <- factor(genes$Chromosome, levels = as.character(1:25))

# 查看数据
head(genes)

# 绘制基因分布图
p <- ggplot(genes, aes(x = Start_MB, y = Gene_Length)) +
  geom_point(color = "#123543", size = 1) +  # 基因的位置点
  facet_wrap(~ Chromosome, scales = "free_x", ncol = 25) +  # 每个染色体分面
  labs(title = "DEA Distribution Across Chromosomes", 
       x = "Position (Mb)", 
       y = "Gene Length (bp)") +
  theme_minimal(base_size = 14) +  # 使用简洁的主题
  theme(
    legend.position = "none",  # 隐藏图例
    axis.text.x = element_text(angle = 90, hjust = 1, size = 10),  # 设置x轴标签
    axis.text.y = element_text(size = 10),  # 设置y轴标签
    strip.text = element_text(size = 12),  # 调整分面标签字体大小
    strip.background = element_blank(),  # 去掉分面背景
    axis.title.x = element_text(size = 14),  # 设置x轴标题
    axis.title.y = element_text(size = 14),  # 设置y轴标题
    panel.spacing = unit(0.5, "lines"),  # 设置面板间距
    panel.grid.major = element_blank(),  # 去除主网格线
    panel.grid.minor = element_blank()   # 去除次网格线
  ) 
dev.off() 
p
# 保存图形
ggsave("log1_distribution_plot.png", plot = p, width = 20, height = 2, dpi = 300)  # 设置合适的宽度和高度
#####################################
###############################3韦恩图##############################33
library(data.table)

green <- fread("ase//matched_results_green.csv", check.names = FALSE, fill = TRUE)
red <- fread("ase//matched_results_red100.csv", check.names = FALSE, fill = TRUE)
yellow <- fread("ase//matched_results_yellow.csv", check.names = FALSE, fill = TRUE)




head(green)
head(red)
head(yellow)

# 提取三个数据框的第一列基因
green_genes <- green$RowNames
red_genes <- red$RowNames
yellow_genes <- yellow$RowNames
print(green_genes)
print(green_genes)

print(green_genes)

# 安装 VennDiagram 包（如果还没安装）
install.packages("VennDiagram")
# 加载 VennDiagram 包
library(VennDiagram)

# 绘制三组基因的韦恩图
venn.plot <- venn.diagram(
  x = list(
    Green = green_genes,
    Red = red_genes,
    Yellow = yellow_genes
  ),
  category.names = c("Green", "Red", "Yellow"),
  filename = NULL,  # 不保存到文件，只显示
  output = TRUE,
  scaled = FALSE,  # 是否按比例缩放区域大小
  fill = c("lightgreen", "lightcoral", "lightyellow"),
  alpha = 0.5,  # 设置透明度
  cex = 1.5,  # 标签字体大小
  cat.cex = 1.5,  # 类别标签字体大小
  cat.pos = c(0, 0, 180),  # 调整标签位置，180度使得yellow标签在顶部
  cat.dist = c(0.05, 0.05, 0.05)  # 类别标签距离圆圈的距离
)

# 显示图形
grid.draw(venn.plot)

dev.off() 

# 显示韦恩图
grid.draw(venn.plot)
##################################################3333
library(ggplot2)
library(reshape2)
head(ase_kaks)

# 重命名列
colnames(ase_kaks) <- c("Geneid", "KaKs")
#kaks_data$Geneid <- sub("-.*", ".0", kaks_data$Geneid)
head(green_genes)
head(red_genes)
head(yellow_genes)
# 创建数据框
green_data <- data.frame(Geneid = green_genes, Color = "Green")
red_data <- data.frame(Geneid = red_genes, Color = "Red")
yellow_data <- data.frame(Geneid = yellow_genes, Color = "Yellow")

# 合并所有颜色数据
genes_data <- rbind(green_data, red_data, yellow_data)

# 合并 Ka/Ks 数据
combined_data <- merge(genes_data, ase_kaks, by = "Geneid", all.x = TRUE)

# 转换 Ka/Ks 为数值类型
combined_data$KaKs <- as.numeric(combined_data$KaKs)
combined_data <- combined_data[!is.na(combined_data$KaKs), ]  # 去除无效值
head(combined_data)
str(combined_data)
library(ggplot2)

ggplot(combined_data, aes(x = Color, y = KaKs, fill = Color)) +
  geom_violin(trim = FALSE) +  # 绘制小提琴图，不裁剪尾部
  geom_boxplot(width = 0.1, outlier.colour = "black", outlier.size = 1) +  # 添加箱线图以显示异常值
  scale_fill_manual(values = c("Red" = "lightcoral", "Green" = "lightgreen", "Yellow" = "lightyellow")) +  # 自定义颜色
  labs(
    title = "Ka/Ks Distribution Across Different Tissues",
    x = "Tissue Color",
    y = "Ka/Ks Value"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    axis.title = element_text(size = 14, face = "bold"),
    axis.text = element_text(size = 12),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

ggplot(combined_data, aes(x = Color, y = KaKs, fill = Color)) +
  geom_violin(trim = FALSE) +  # 绘制小提琴图，不裁剪尾部
  scale_fill_manual(values = c("Red" = "lightcoral", "Green" = "lightgreen", "Yellow" = "lightyellow")) +  # 自定义颜色
  labs(
    title = "Ka/Ks Distribution Across Different Tissues",
    x = "Tissue Color",
    y = "Ka/Ks Value"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    axis.title = element_text(size = 14, face = "bold"),
    axis.text = element_text(size = 12),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

ggplot(combined_data, aes(x = Color, y = KaKs, fill = Color)) +
  geom_violin(trim = FALSE) +
  scale_fill_manual(values = c("Red" = "red", "Green" = "green", "Yellow" = "yellow")) +
  labs(
    title = "Ka/Ks Distribution Across Different Tissues",
    x = "Tissue Color",
    y = "Ka/Ks Value"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    axis.title = element_text(size = 14, face = "bold"),
    axis.text = element_text(size = 12),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )
