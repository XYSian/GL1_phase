library(DESeq2)

count_data_new <- hap1_exp[, c("Geneid", "Red1.x", "Red2.x", "Red3.x", 
                               "Green1.x", "Green2.x", "Green3.x", 
                               "Yellow1.x", "Yellow2.x", "Yellow3.x")]
count_data_new <- hap2_exp[, c("Geneid", "Red1.y", "Red2.y", "Red3.y", 
                               "Green1.y", "Green2.y", "Green3.y", 
                               "Yellow1.y", "Yellow2.y", "Yellow3.y")]
# 设置 Geneid 为行名
rownames(count_data_new) <- count_data_new$Geneid
count_data_new <- count_data_new[,-1]

# 将所有值转换为整数
count_data_new <- round(count_data_new)

# 创建条件信息（颜色信息）
conditions_new <- factor(c(rep("Red", 3), rep("Green", 3), rep("Yellow", 3)))

# 创建 DESeqDataSet 对象
dds_new <- DESeqDataSetFromMatrix(countData = as.matrix(count_data_new),
                                  colData = data.frame(condition = conditions_new),
                                  design = ~ condition)

# 运行差异表达分析
dds_new <- DESeq(dds_new)

# 提取比较组结果
res_Red_vs_Green_new <- results(dds_new, contrast = c("condition", "Red", "Green"))
res_Red_vs_Yellow_new <- results(dds_new, contrast = c("condition", "Red", "Yellow"))
res_Green_vs_Yellow_new <- results(dds_new, contrast = c("condition", "Green", "Yellow"))

# 查看红绿比较的显著差异基因
res_Red_vs_Green_sig_new <- subset(res_Red_vs_Green_new, padj < 0.05 & abs(log2FoldChange) >= 1)
res_Red_vs_Green_sig_new
res_Red_vs_Green_sig_new$RowNames <- rownames(res_Red_vs_Green_sig_new)
res_Red_vs_Green_sig_new <- as.data.frame(res_Red_vs_Green_sig_new)

res_Red_vs_Green_sig_new <- merge(res_Red_vs_Green_sig_new, merged_anno, by.x = "RowNames", by.y = "Geneid", all.x = TRUE)
res_Red_vs_Green_sig_new <- merge(res_Red_vs_Green_sig_new, wide_data[, c("Geneid", "Haplotype2")], 
                                   by.x = "RowNames", by.y = "Geneid", all.x = TRUE)
write.xlsx(as.data.frame(res_Red_vs_Green_sig_new), 
           file = "res_Red_vs_Green_sig_new.xlsx", 
           rowNames = TRUE)
# 查看红白比较的显著差异基因
res_Red_vs_Yellow_sig_new <- subset(res_Red_vs_Yellow_new, padj < 0.05 & abs(log2FoldChange) >= 1)
res_Red_vs_Yellow_sig_new
res_Red_vs_Yellow_sig_new$RowNames <- rownames(res_Red_vs_Yellow_sig_new)
res_Red_vs_Yellow_sig_new <- as.data.frame(res_Red_vs_Yellow_sig_new)

res_Red_vs_Yellow_sig_new <- merge(res_Red_vs_Yellow_sig_new, merged_anno, by.x = "RowNames", by.y = "Geneid", all.x = TRUE)
res_Red_vs_Yellow_sig_new <- merge(res_Red_vs_Yellow_sig_new, wide_data[, c("Geneid", "Haplotype2")], 
                               by.x = "RowNames", by.y = "Geneid", all.x = TRUE)
write.xlsx(as.data.frame(res_Red_vs_Yellow_sig_new), 
           file = "res_Red_vs_Yellow_sig_new.xlsx", 
           rowNames = TRUE)
# 查看绿白比较的显著差异基因
res_Green_vs_Yellow_sig_new <- subset(res_Green_vs_Yellow_new, padj < 0.05 & abs(log2FoldChange) >= 1)
res_Green_vs_Yellow_sig_new
res_Green_vs_Yellow_sig_new$RowNames <- rownames(res_Green_vs_Yellow_sig_new)
res_Green_vs_Yellow_sig_new <- as.data.frame(res_Green_vs_Yellow_sig_new)

res_Green_vs_Yellow_sig_new <- merge(res_Green_vs_Yellow_sig_new, merged_anno, by.x = "RowNames", by.y = "Geneid", all.x = TRUE)
res_Green_vs_Yellow_sig_new <- merge(res_Green_vs_Yellow_sig_new, wide_data[, c("Geneid", "Haplotype2")], 
                                   by.x = "RowNames", by.y = "Geneid", all.x = TRUE)
write.xlsx(as.data.frame(res_Green_vs_Yellow_sig_new), 
           file = "res_Green_vs_Yellow_sig_new.xlsx", 
           rowNames = TRUE)
# 绘制火山图
library(ggplot2)
ggplot(as.data.frame(res_Red_vs_Green_new), aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point() +
  theme_minimal() +
  ggtitle("Red vs Green Volcano Plot")

install.packages("openxlsx")
library(openxlsx)

# 将 count_data_new 写入 Excel 文件
write.xlsx(count_data_new, file = "count_data_new.xlsx", rowNames = TRUE)
