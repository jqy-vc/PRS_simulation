# 1. 加载必要的包
library(ggplot2)
library(dplyr)
library(readxl) 
library(tidyr) 

# --- 设定参数 ---
excel_file_path <- "/data2/qiuyue/PRS_simulation/PRS_LD_Overlap_Dprime_Simulation_Results.xlsx"
H2_BENCHMARK <- 0.5 
BAR_WIDTH <- 0.7        # 柱子的宽度
DODGE_WIDTH <- 0.8      # position_dodge 的总宽度

# 2. 读取原始数据文件
data_raw <- read_excel(excel_file_path, sheet = 1) %>% 
  as.data.frame() 

# 3. 数据处理和重塑 
observed_summary <- data_raw %>%
  mutate(across(c(Dprime_Min, Dprime_Max, Observed_R_Squared, Predicted_R_Squared), as.numeric)) %>%
  group_by(Dprime_Min, Dprime_Max, Predicted_R_Squared) %>%
  summarise(
    R_squared_Value = mean(Observed_R_Squared, na.rm = TRUE), 
    R_squared_SE = sd(Observed_R_Squared, na.rm = TRUE) / sqrt(n()), 
    .groups = 'drop'
  ) %>%
  mutate(R_squared_Type = "mean_rsq_obs") 

predicted_summary <- observed_summary %>%
  select(Dprime_Min, Dprime_Max, Predicted_R_Squared) %>% 
  distinct() %>% 
  mutate(
    R_squared_Value = Predicted_R_Squared, 
    R_squared_SE = 0,                       
    R_squared_Type = "mean_pred_rsq"          
  )

plot_data_long <- bind_rows(observed_summary, predicted_summary) %>%
  mutate(combo_label = paste0("Dprim=", Dprime_Min, "-", Dprime_Max)) %>%
  arrange(Dprime_Min) %>%
  mutate(combo_label = factor(combo_label, levels = unique(combo_label)))


# --- 4. 生成图表 (已移除 geom_col 中的 'stat' 参数) ---

# 预先计算用于 annotate() 的 x 轴位置
max_x_position <- length(levels(plot_data_long$combo_label))

plot_dprime_bar_final <- ggplot(plot_data_long, aes(x = combo_label, y = R_squared_Value, fill = R_squared_Type)) +
  
  # A. 绘制并排柱状图 (已移除 stat = "identity")
  geom_bar(stat = "identity", position = position_dodge(width = 0.8)) +
  
# B. 添加误差棒 (只在蓝色柱子上)
geom_errorbar(
  data = filter(plot_data_long),
  aes(
    ymin = R_squared_Value - R_squared_SE, 
    ymax = R_squared_Value + R_squared_SE,
    group = R_squared_Type   # 🔑 关键：把 fill 映射也带上，和柱子保持一致
  ),
  position = position_dodge(width = DODGE_WIDTH),
  width = 0.2,
  color = "darkgray"
)+
  
  # C. 添加数值标签
  geom_text(
    aes(label = round(R_squared_Value, 4), group = R_squared_Type), 
    vjust = -0.5,
    position = position_dodge(width = DODGE_WIDTH), 
    size = 3.0, 
    color = "black"
  ) +
  
  # D. 添加 h2 基准虚线
  geom_hline(
    yintercept = H2_BENCHMARK, 
    linetype = "dashed", 
    color = "blue", 
    linewidth = 0.8
  ) +
  
  # E. 添加 h2 标签 (使用 annotate)
  annotate(
    "text", 
    x = max_x_position + 0.5, 
    y = H2_BENCHMARK, 
    label = paste0("h2 = ", H2_BENCHMARK),
    hjust = 1.0, vjust = -0.5, color = "blue", size = 5
  ) +
  
  # F. 标签和主题设置
  labs(
    title = expression(paste("Observed vs. Predicted ", R^2, " by LD D' Range")), 
    x = expression(paste("LD D' Range")),
    y = expression(paste("PRS ", R^2, " Value")),
    fill = "R-squared Type"
  ) +
  
  # G. 颜色和图例设置
  scale_fill_manual(
    values = c("mean_rsq_obs" = "skyblue", "mean_pred_rsq" = "lightcoral"),
    labels = c("Theoritical R-squared", "Observed R-squared")
  ) +
  
  # H. 主题调整
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1), 
    legend.position = "bottom",
    plot.title = element_text(hjust = 0.5, face = "bold")
  ) +
  
  # I. 调整 Y 轴范围
  coord_cartesian(
    ylim = c(0, max(plot_data_long$R_squared_Value + plot_data_long$R_squared_SE) * 1.1+0.2)
  )

# 5. 显示图形
print(plot_dprime_bar_final)

# 保存到 PDF 文件
ggsave(
  filename = "/data2/qiuyue/PRS_simulation/plot_dprime_bar_final.pdf",  # 输出文件名
  plot = plot_dprime_bar_final,           
  width = 10, height = 6,                
  device = cairo_pdf                      
)
