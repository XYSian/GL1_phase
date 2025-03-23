ggplot(collinearity_data, aes(x = Genome_Position_Mb, y = Position, color = as.factor(group))) +
  geom_point(size = 2) +
  labs(x = "Genome Position (Mb)", y = "Linkage Map Position",
       title = "Collinearity between Genome and Linkage Map") +
  theme_minimal() +
  scale_color_manual(values = colors) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  facet_wrap(~ group, scales = "free")  # 按group分面，自由调整坐标轴


ggplot(collinearity_data, aes(x = Genome_Position_Mb, y = cM_accumulated, color = as.factor(group))) +
  geom_point(size = 2) +
  labs(x = "Genome Position (Mb)", y = "Cumulative Genetic Distance (cM)",
       title = "Collinearity between Genome and Linkage Map") +
  theme_minimal() +
  scale_color_manual(values = colors) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  facet_wrap(~ group, scales = "free_y")  # 每个组分开显示，自由调整y轴刻度
