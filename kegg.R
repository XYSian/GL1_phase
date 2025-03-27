kegg <- read.table("query.ko.plant.txt", sep="\t", header=FALSE, fill=TRUE, quote="", stringsAsFactors=FALSE)

kegg_cleaned <- kegg %>%
  filter(V2 != "")
# 查看结果
print(kegg_cleaned)
# 提取 kegg_cleaned 的第二列（KEGG ID）
kegg_ids <- kegg_cleaned$V2
# 将 KEGG ID 输出到文件中
write.table(kegg_ids, file = "kegg_ids.txt", quote = FALSE, sep = "\t", row.names = FALSE, col.names = FALSE)
#install.packages("KEGGREST")

library(KEGGREST)
print(kegg_ids)
# 初始化一个空列表来存储结果
kegg_names <- list()
library(clusterProfiler)

x <- bitr_kegg(kegg_ids,"kegg","Path","ko") 

print(x)
ko_ids <- x$Path
print(ko_ids)
ko_ids <- unique(ko_ids)
# 将 ko_ids 列表分批
batch_size <- 100
ko_ids_list <- split(ko_ids, ceiling(seq_along(ko_ids) / batch_size))
head(ko_ids_list)

#用 lapply 批量查询
y <- lapply(ko_ids_list, function(batch) keggList(batch))

print(y)
head(y)

kegg_df <- data.frame(
  KO_ID = names(y),
  KO_Name = unlist(y),
  stringsAsFactors = FALSE
)
# 去掉 rownames 中的前缀 "1."
rownames(kegg_df) <- gsub("^\\d+\\.", "", rownames(kegg_df))

# 将 KO_ID 列替换为 rownames
kegg_df$KO_ID <- rownames(kegg_df)

# 查看结果
print(kegg_df)
library(dplyr)

#result <- left_join(x, kegg_df, by = c("Path" = "KO_ID"))
kegg_term2gene <- data.frame(kegg_cleaned$V2,kegg_cleaned$V1)
#gene <- as.factor(kegg_cleaned$V1)
#head(gene)
head(kegg_term2name)

# 生成 koid2name 数据框
koid2name <- kegg_df %>%
  rename(KO_ID = KO_ID, KO_Name = KO_Name)
print(koid2name)

kegg_to_ko <- x %>%
  rename(KEGG_ID = kegg, KO_ID = Path)

# 2. 将 KEGG ID 映射到基因 ID
kegg_to_gene <- kegg_term2gene %>%
  rename(KEGG_ID = kegg_cleaned.V2, Gene_ID = kegg_cleaned.V1)

# 3. 将 KO ID 映射到基因 ID
koid_to_gene <- kegg_to_ko %>%
  left_join(kegg_to_gene, by = "KEGG_ID") %>%
  select(KO_ID, Gene_ID)


# 从 x 数据框中提取 KEGG_ID 和 KO_ID 的映射
kegg_to_ko <- x %>%
  select(KEGG_ID = kegg, KO_ID = Path)

print(koid_to_gene)
head(koid2name)
print(koid2name)
print(gene_green)

z <- enricher(
  gene = Green_vs_Yellow, 
  TERM2GENE = koid_to_gene, 
  TERM2NAME = koid2name,
  pvalueCutoff = 0.05, 
  pAdjustMethod = "BH"
)

# 查看结果
print(z)
barplot(z, showCategory = 20)+ ggtitle("KEGG Pathway Enrichment Analysis Of Yellow Leaf DEA")

# 绘制条形图
write.csv(as.data.frame(z),"Green_vs_Yellow.kegg.csv",row.names = F)
zz <- read.csv("red.green.kegg.csv")
print(zz)
# 绘制条形图并添加标题
top10 <- head(zz[order(zz$p.adjust), ], 16)

# 根据 p.adjust 排序 Description（显著性越高越靠前）
top10$Description <- factor(top10$Description, levels = top10$Description[order(top10$p.adjust, decreasing = TRUE)])

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

