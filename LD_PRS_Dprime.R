# ----------------------------------------------------
# 0. 参数设置 and LD 定义
source("/data2/qiuyue/PRS_simulation/source.R")

m <- 1000  # 标记总数 (M)
n1 <- 1000  # 训练集总大小
n2 <- 500   # 测试集总大小
w_overlap <- 0.25 # 重叠比例
h2 <- 0.5  # 遗传力
n_overlap <- floor(n2 * w_overlap) # 重叠样本数
n_total <- n1 + n2 - n_overlap # 实际需要的总样本量 (N)

# 基础等位基因频率 (Freq) - 在整个实验中保持不变
set.seed(42)
freq <- runif(m, min = 0.2, max = 0.5)

# --- 定义 LD 情景参数 ---
# 8 种 D' 分布区间：[min, max]
dprime_params <- list(
  c(0.01, 0.02), c(0.05, 0.1), c(0.1, 0.2), c(0.2, 0.3),
  c(0.5, 0.6), c(0.7, 0.8), c(0.8, 0.9), c(0.7, 0.95)
)
num_simulations <- 30 # 每种情景的模拟次数

# ----------------------------------------------------
# 1. 实验循环主体
# ----------------------------------------------------

# 存储结果的数据框：新增 Predicted_R_Squared 字段
results_df <- data.frame(
  Dprime_Min = numeric(),
  Dprime_Max = numeric(),
  Simulation_Run = integer(),
  Observed_R_Squared = numeric(), # 重命名 R_Squared 以区分
  Predicted_R_Squared = numeric(),
  stringsAsFactors = FALSE
)

# ----------------------------------------------------
# 2. 预先计算理论预估值 (predRsq)
# ----------------------------------------------------
# 注意：该公式是基于样本量和遗传力，不随 LD 变化，因此只需计算一次。
# 公式：predRsq = (h2 + w*m/n1)^2 / (w*h2*m/n1 + (h2 + m/n1) * (1 + w*m/n1))
numerator <- (h2 + w_overlap * m / n1)^2
denominator <- (w_overlap * h2 * m / n1) +
  (h2 + m / n1) * (1 + w_overlap * m / n1)
fixed_pred_rsq <- numerator / denominator

cat(sprintf("\n根据公式计算的理论预估 R² (固定值)：%.4f\n", fixed_pred_rsq))


# 外部循环：遍历 8 种 D' 分布参数
for (p_idx in seq_along(dprime_params)) {

  min_d <- dprime_params[[p_idx]][1]
  max_d <- dprime_params[[p_idx]][2]

  cat(sprintf("\n--- 正在运行 D' 分布: U(%.2f, %.2f) ---\n",
    min_d, max_d))

  # 内部循环：对每种 LD 情景进行 30 次模拟
  for (run in seq_len(num_simulations)) {

    set.seed(420 + run)

    # 每次模拟重新生成 D' 和 LD 矩阵
    Dprime_current <- runif(m, min_d, max_d)
    current_ld <- DprimetoD(freq, Dprime_current)

    x_matrix_raw <- GenerateLDGeno(
      freq = freq, ld = current_ld, N = n_total)
    x_matrix <- scale(x_matrix_raw) # Z-score 标准化

    x_train <- x_matrix[1:n1, ]
    start_test <- n1 + 1 - n_overlap
    x_test <- x_matrix[start_test:n_total, ]

    beta_true <- rnorm(m, mean = 0, sd = 1)
    G <- x_matrix %*% beta_true
    var_G <- var(G)

    var_E <- var_G / h2 * (1 - h2)
    ei <- rnorm(n_total, mean = 0, sd = sqrt(var_E))
    y <- G + ei

    Y_train <- y[1:n1]
    y_test <- y[start_test:n_total]

    individual_snp_betas <- numeric(m)
    for (i in seq_len(m)) {
      current_model <- lm(Y_train ~ x_train[, i])
      individual_snp_betas[i] <- coef(current_model)[2]
    }

    # ----------------------------------------------------
    # 6. PRS 计算与评估
    # ----------------------------------------------------

    y_hat <- x_test %*% individual_snp_betas

    cor_obs <- cor(y_hat, y_test, method = "pearson")
    Observed_R_Squared <- cor_obs^2

    # 存储结果：注意顺序与新增字段
    results_df[nrow(results_df) + 1, ] <- list(
      min_d,
      max_d,
      run,
      Observed_R_Squared,
      fixed_pred_rsq
    )
  }
}

# ----------------------------------------------------
# 7. 最终结果和分析 (包含 Excel 保存)
# ----------------------------------------------------

print("--- 样本重叠和不同 LD 强度 (D' 分布) 对简单 PRS 预测 R² 的影响 ---")
print(results_df)

# 计算并打印每种 LD 情景的平均 R²
# 注意：我们现在计算 Observed_R_Squared 的平均值，并包含固定预估值
summary_results <- aggregate(
  Observed_R_Squared ~ Dprime_Min + Dprime_Max,
  data = results_df,
  FUN = mean
)

# 添加理论预估 R² 列 (它是一个固定值)
summary_results$Predicted_R_Squared <- fixed_pred_rsq
colnames(summary_results) <- c(
  "Dprime_Min",
  "Dprime_Max",
  "Mean_Observed_R_Squared",
  "Predicted_R_Squared"
)

print("\n--- 每种 D' 情景的平均 R² 及理论预估 R² ---")
print(summary_results)

# ----------------------------------------------------
# 保存结果到 Excel
# ----------------------------------------------------
library(writexl)

output_file <- "/data2/qiuyue/PRS_simulation/PRS_LD_Overlap_Dprime_Simulation_Results.xlsx"

write_xlsx(
  list(
    Raw_Runs = results_df,
    Summary_Averages = summary_results
  ),
  path = output_file
)

cat(sprintf("\n--- 结果已成功保存到 Excel 文件： %s ---\n",
  output_file))
