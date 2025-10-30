# ----------------------------------------------------
# 0. 参数设置 and LD 定义
# ----------------------------------------------------
source("/data2/qiuyue/PRS_simulation/source.R")

# 固定参数
m <- 1000 # 标记总数 (M) SNP n个数
n2 <- 500 # 测试集总大小 (固定)
h2 <- 0.5 # 遗传力 (固定)
num_simulations <- 30 # 每种情景的模拟次数

# 遍历参数集合
n1_params <- c(1000, 2000, 5000, 10000, 100000) # 训练集大小
w_params <- c(0.1, 0.2, 0.4, 0.8) # 样本重叠比例
dprime_params <- list(
	c(0.01, 0.02), c(0.05, 0.1), c(0.1, 0.2), c(0.2, 0.3),
	c(0.5, 0.6), c(0.7, 0.8), c(0.8, 0.9), c(0.7, 0.95)
)

# 基础等位基因频率 (Freq) - 在整个实验中保持不变
set.seed(42)
freq <- runif(m, min = 0.2, max = 0.5)

# ----------------------------------------------------
# 1. 存储结果的数据框：新增 n1 和 w_overlap 字段
# ----------------------------------------------------
results_df <- data.frame(
	N_Train = numeric(),
	W_Overlap = numeric(),
	Dprime_Min = numeric(),
	Dprime_Max = numeric(),
	Simulation_Run = integer(),
	Observed_R_Squared = numeric(),
	Predicted_R_Squared = numeric(),
	stringsAsFactors = FALSE
)

# ----------------------------------------------------
# 2. 三层嵌套循环主体 (n1 -> w -> D')
# ----------------------------------------------------

# 第一层循环：遍历训练集大小 n1
for (n1 in n1_params) {

cat(sprintf("\n============= N_TRAIN = %d =============\n", n1))

# 第二层循环：遍历重叠比例 w
for (w_overlap in w_params) {

n_overlap <- floor(n2 * w_overlap) # 根据当前 w 计算重叠样本数
n_total <- n1 + n2 - n_overlap # 实际需要的总样本量

cat(sprintf("\n----- W_OVERLAP = %.2f (N_OVERLAP = %d) -----\n",
						w_overlap, n_overlap))

# 预先计算理论预估值 (predRsq)
# 公式：predRsq = (h2 + w*m/n1)^2 /
#        (w*h2*m/n1 + (h2 + m/n1) * (1 + w*m/n1))
numerator <- (h2 + w_overlap * m / n1)^2
denominator <- (w_overlap * h2 * m / n1) +
	(h2 + m / n1) * (1 + w_overlap * m / n1)
fixed_pred_rsq <- numerator / denominator

cat(sprintf(" 理论预估 R² (固定值)：%.4f\n", fixed_pred_rsq))

# 第三层循环：遍历 D' 分布参数
for (p_idx in seq_along(dprime_params)) {

min_d <- dprime_params[[p_idx]][1]
max_d <- dprime_params[[p_idx]][2]

cat(sprintf("\n--- D' 分布: U(%.2f, %.2f) ---\n", min_d, max_d))

# 第四层循环：重复模拟
for (run in seq_len(num_simulations)) {

set.seed(420 + run) # 设置新的随机种子以确保每轮独立

# ----------------------------------------------------
# 3. 基因型生成及 PRS 评估
# ----------------------------------------------------

# 每次模拟重新生成 D' 和 LD 矩阵
Dprime_current <- runif(m, min_d, max_d)
current_ld <- DprimetoD(freq, Dprime_current)

x_matrix_raw <- GenerateLDGeno(freq = freq, ld = current_ld,
																N = n_total)
x_matrix <- scale(x_matrix_raw) # Z-score 标准化

# 划分数据集
x_train <- x_matrix[1:n1, ]
start_test <- n1 + 1 - n_overlap
x_test <- x_matrix[start_test:n_total, ]

# 模拟表型
beta_true <- rnorm(m, mean = 0, sd = 1)
G <- x_matrix %*% beta_true
var_G <- var(G)
var_E <- var_G / h2 * (1 - h2)
ei <- rnorm(n_total, mean = 0, sd = sqrt(var_E))
y <- G + ei

Y_train <- y[1:n1]
y_test <- y[start_test:n_total]

# GWAS (在训练集上)
individual_snp_betas <- numeric(m)
for (i in seq_len(m)) {
current_model <- lm(Y_train ~ x_train[, i])
individual_snp_betas[i] <- coef(current_model)[2]
}

# PRS 计算与 R^2 评估
y_hat <- x_test %*% individual_snp_betas
cor_obs <- cor(y_hat, y_test, method = "pearson")
Observed_R_Squared <- cor_obs^2

# 存储结果
results_df[nrow(results_df) + 1, ] <- list(
	n1, # N_Train
	w_overlap, # W_Overlap
	min_d,
	max_d,
	run,
	Observed_R_Squared,
	fixed_pred_rsq
)
}
}
}
}

# ----------------------------------------------------
# 4. 最终结果和分析 (包含 Excel 保存)
# ----------------------------------------------------

# 计算每种组合 (n1, w, D') 的平均 R²
summary_results <- aggregate(
	Observed_R_Squared ~ N_Train + W_Overlap + Dprime_Min + Dprime_Max,
	data = results_df,
	FUN = mean
)
colnames(summary_results)[5] <- "Mean_Observed_R_Squared"

# 重新计算 Predicted_R_Squared，因为它是 n1 和 w 的函数
summary_results$Predicted_R_Squared <- numeric(nrow(summary_results))

for (i in seq_len(nrow(summary_results))) {
current_n1 <- summary_results$N_Train[i]
current_w <- summary_results$W_Overlap[i]

# 重新计算公式：predRsq = (h2 + w*m/n1)^2 /
#        (w*h2*m/n1 + (h2 + m/n1) * (1 + w*m/n1))
numerator <- (h2 + current_w * m / current_n1)^2
denominator <- (current_w * h2 * m / current_n1) +
	(h2 + m / current_n1) * (1 + current_w * m / current_n1)
summary_results$Predicted_R_Squared[i] <- numerator / denominator
}

print("--- 多参数遍历：平均 R² 及理论预估值 ---")
print(summary_results)


# ----------------------------------------------------
# 保存结果到 Excel
# ----------------------------------------------------
library(writexl)

output_file <- "/data2/qiuyue/PRS_simulation/PRS_LD_Overlap_N_W_Dprime_Simulation_Results.xlsx"

write_xlsx(
	list(
		Raw_Runs = results_df,
		Summary_Averages = summary_results
	),
	path = output_file
)

cat(sprintf("\n--- 结果已成功保存到 Excel 文件： %s ---\n", output_file))
