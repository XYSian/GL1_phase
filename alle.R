# 读取注释数据
merged_anno0 <- read.table("merged0_final", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
colnames(merged_anno0) <- c("Haplotype1", "AnotherColumn", "GeneAnnotation")

merged_anno1 <- read.table("merged1.final", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
colnames(merged_anno1) <- c("Haplotype1", "AnotherColumn", "GeneAnnotation")

# 读取等位基因数据
file_path99 <- "E:\\rstudio\\work\\ASE\\alle\\01-01.xlsx"
allele_data99 <- read_excel(file_path99)
colnames(allele_data99) <- c("Haplotype1", "Haplotype2")

file_path100 <- "E:\\rstudio\\work\\ASE\\alle\\02-02.xlsx"
allele_data100 <- read_excel(file_path100)
colnames(allele_data100) <- c("Haplotype1", "Haplotype2")
# 读取数据（假设数据已经在环境中）
# allele_data99 <- read_excel("path_to_allele_data99.xlsx")
# merged_anno0 <- read.table("path_to_merged_anno0.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
# merged_anno1 <- read.table("path_to_merged_anno1.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# 合并 allele_data99 和 merged_anno0，根据 Haplotype1 添加 AnotherColumn0
allele_data99_annotated_0 <- merge(allele_data99, merged_anno0, by.x = "Haplotype1", by.y = "Haplotype1", all.x = TRUE)
colnames(allele_data99_annotated_0)[which(names(allele_data99_annotated_0) == "AnotherColumn")] <- "AnotherColumn0"
# 合并 allele_data99 和 merged_anno1，根据 Haplotype2 添加 AnotherColumn1
allele_data99_annotated_1 <- merge(allele_data99, merged_anno1, by.x = "Haplotype2", by.y = "Haplotype1", all.x = TRUE)
colnames(allele_data99_annotated_1)[which(names(allele_data99_annotated_1) == "AnotherColumn")] <- "AnotherColumn1"
# 合并两个注释结果
final_annotated_data99 <- merge(allele_data99_annotated_0, allele_data99_annotated_1[, c("Haplotype2", "AnotherColumn1")], by = "Haplotype2", all.x = TRUE)
# 查看结果
print(final_annotated_data99)



allele_data100_annotated_0 <- merge(allele_data100, merged_anno0, by.x = "Haplotype1", by.y = "Haplotype1", all.x = TRUE)
colnames(allele_data100_annotated_0)[which(names(allele_data100_annotated_0) == "AnotherColumn")] <- "AnotherColumn0"
# 合并 allele_data100 和 merged_anno1，根据 Haplotype2 添加 AnotherColumn1
allele_data100_annotated_1 <- merge(allele_data100, merged_anno1, by.x = "Haplotype2", by.y = "Haplotype1", all.x = TRUE)
colnames(allele_data100_annotated_1)[which(names(allele_data100_annotated_1) == "AnotherColumn")] <- "AnotherColumn1"
# 合并两个注释结果
final_annotated_data100 <- merge(allele_data100_annotated_0, allele_data100_annotated_1[, c("Haplotype2", "AnotherColumn1")], by = "Haplotype2", all.x = TRUE)
print(final_annotated_data100)
