# --- 1. 加载包 ---
library(ggplot2)
library(dplyr)
library(readxl)
library(tidyr)
library(stringr)

# --- 2. 设定参数 ---
excel_file_path <- "/data2/qiuyue/PRS_simulation/PRS_LD_Overlap_N_W_Dprime_Simulation_Results.xlsx"
output_dir <- "/data2/qiuyue/PRS_simulation/"   # 输出目录
H2_BENCHMARK <- 0.5
BAR_WIDTH <- 0.7
DODGE_WIDTH <- 0.8

# --- 3. 读取数据 ---
data_raw <- read_excel(excel_file_path, sheet = 1) %>%
  as.data.frame() %>%
  mutate(
    across(c(Dprime_Min, Dprime_Max, Observed_R_Squared, Predicted_R_Squared), as.numeric),
    combo_label = paste0("n1=",N_Train , "_w=", W_Overlap)
  )
head(data_raw)
# --- 4. 获取所有 D′ 区间 ---
dprime_ranges <- data_raw %>%
  distinct(Dprime_Min, Dprime_Max) %>%
  arrange(Dprime_Min)

# --- 5. 为每个 D′ 区间绘制图并保存 (替换为 bas_with_overlap_h_plot 风格的误差棒处理) ---
for (i in seq_len(nrow(dprime_ranges))) {

  dmin <- dprime_ranges$Dprime_Min[i]
  dmax <- dprime_ranges$Dprime_Max[i]

  subset_data <- data_raw %>%
    filter(Dprime_Min == dmin, Dprime_Max == dmax) %>%
    mutate(
      N_Train = as.numeric(N_Train),
      W_Overlap = as.numeric(W_Overlap)
    )

  # Observed 汇总
  observed_summary <- subset_data %>%
    group_by(N_Train, W_Overlap) %>%
    summarise(
      R_squared_Value = mean(Observed_R_Squared, na.rm = TRUE),
      R_squared_SD = sd(Observed_R_Squared, na.rm = TRUE),
      n_rep = dplyr::n(),
      .groups = "drop"
    ) %>%
    arrange(N_Train, W_Overlap) %>%
    mutate(
      R_squared_Type = "mean_rsq_obs",
      combo_label = paste0("n1=", N_Train, "\nw=", W_Overlap)
    )

  # Predicted 汇总（保持列名一致，sd 若为 NA 则置 0）
  predicted_summary <- subset_data %>%
    group_by(N_Train, W_Overlap) %>%
    summarise(
      R_squared_Value = mean(Predicted_R_Squared, na.rm = TRUE),
      R_squared_SD = sd(Predicted_R_Squared, na.rm = TRUE),
      n_rep = dplyr::n(),
      .groups = "drop"
    ) %>%
    arrange(N_Train, W_Overlap) %>%
    mutate(
      R_squared_Type = "mean_pred_rsq",
      combo_label = paste0("n1=", N_Train, "\nw=", W_Overlap)
    )
  predicted_summary$R_squared_SD[is.na(predicted_summary$R_squared_SD)] <- 0

  # 合并并固定 x 顺序（n1 增、w 增）
  combo_levels <- observed_summary %>%
    arrange(N_Train, W_Overlap) %>%
    pull(combo_label) %>% unique()

  plot_data_long <- bind_rows(observed_summary, predicted_summary) %>%
    mutate(
      combo_label = factor(combo_label, levels = combo_levels),
      R_squared_Type = factor(R_squared_Type, levels = c("mean_rsq_obs", "mean_pred_rsq")),
      # 统一用 SD 字段，绘制误差棒用 sd/sqrt(n)
      R_squared_SE = R_squared_SD / sqrt(pmax(1, n_rep)),
      label_text = ifelse(R_squared_Type == "mean_rsq_obs",
                          paste0("Obs: ", sprintf("%.4f", R_squared_Value)),
                          paste0("Pred: ", sprintf("%.4f", R_squared_Value)))
    )

  # dodge 保持一致，group 明确
  dodge_pos <- position_dodge2(width = DODGE_WIDTH, preserve = "single")

  # 画图：每根柱子都有误差棒，标签带前缀并与对应柱对齐
  p <- ggplot(plot_data_long, aes(x = combo_label, y = R_squared_Value,
                                  fill = R_squared_Type, group = R_squared_Type)) +
    geom_col(stat = "identity", position = dodge_pos, width = BAR_WIDTH) +
    geom_errorbar(aes(ymin = pmax(0, R_squared_Value - R_squared_SE),
                      ymax = R_squared_Value + R_squared_SE),
                  position = position_dodge(width = 0.65), width = 0.2, color = "darkgray", size = 0.5) +
    geom_text(aes(label = label_text),
              position = dodge_pos, vjust = -0.5, size = 3, color = "black") +
    geom_hline(yintercept = H2_BENCHMARK, linetype = "dashed", color = "blue", linewidth = 1) +
    annotate("text", x = Inf, y = H2_BENCHMARK, label = paste0("h2 = ", H2_BENCHMARK),
             hjust = 1.2, vjust = -0.6, color = "blue", size = 4, inherit.aes = FALSE) +
    labs(
      title = bquote("Observed vs. Predicted " ~ R^2 ~ " (D'" ~ " range = " ~ .(dmin) ~ "-" ~ .(dmax) ~ ")"),
      subtitle = "n2: 500, m: 1000, n_simulation: 30",
      x = "Simulation parameters (n1 ascending, w ascending)",
      y = expression(R^2 ~ " Value"),
      fill = "R-squared Type"
    ) +
    scale_fill_manual(values = c("mean_pred_rsq" = "lightcoral", "mean_rsq_obs" = "skyblue"),
                      labels = c("Predicted R-squared", "Observed R-squared")) +
    theme_minimal(base_size = 13) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "bottom",
      plot.title = element_text(hjust = 0.5, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5)
    ) +
    coord_cartesian(ylim = c(0, max(plot_data_long$R_squared_Value + plot_data_long$R_squared_SE, na.rm = TRUE) * 1.10))

  # 保存
  output_file <- file.path(output_dir, paste0("plot_Dprime_", dmin, "-", dmax, ".pdf"))
  ggsave(filename = output_file, plot = p, width = 10, height = 6, device = cairo_pdf)

  message("✅ 已保存图像：", output_file)
}

message("🎯 所有 D′ 区间图已绘制完成并保存到 ", output_dir)
