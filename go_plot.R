library(ggplot2)
#install.packages("viridis")
library(viridis)
# 绘制 dotplot 并提取数据
p <- dotplot(x, showCategory=20)

# 修改图形对象以添加颜色渐变
p + scale_color_gradient(low = "blue", high = "red") +
  theme_minimal() + 
  ggtitle("GO Enrichment Analysis of Yellow Leaf DEA") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

 
p <- dotplot(x, showCategory=20, color = "pvalue")

# 修改图形对象以添加颜色渐变
p + scale_color_gradient(low = "blue", high = "red") +
  theme_minimal() + 
  ggtitle("GO Enrichment Analysis of Red_vs_Green") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

#########################################################################3# 创建 cnetplot############################
p <- cnetplot(x, categorySize="pvalue")

# 自定义图形的主题和标签字体大小
p + theme(
  plot.margin = margin(1, 1, 1, 1, "cm"),  # 调整图形边距
  axis.text = element_text(size = 8),       # x 和 y 轴标签字体大小
  axis.title = element_text(size = 10),     # x 和 y 轴标题字体大小
  legend.text = element_text(size = 8),     # 图例文本字体大小
  legend.title = element_text(size = 1),   # 图例标题字体大小
  plot.title = element_text(size = 12),     # 图形标题字体大小
  # 调整节点标签（基因）字体大小
  text = element_text(size = 5)              # 调整图形中的所有文本元素
