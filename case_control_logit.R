library(stats)
library(pROC)

# 1. 参数设置与总体基因型生成
num_loci <- 100
min_maf <- 0.05
max_maf <- 0.5
n <- 100000
N_cases <- 100
N_controls <- 100
alpha_intercept <- -3
target_OR <- 1

set.seed(1)
p_pop_true <- runif(num_loci, min_maf, max_maf)

true_beta_j <- log(target_OR)
true_betas <- rep(true_beta_j, num_loci)

# 1.1 生成总人群基因型数据 (G)
genotypes <- matrix(NA, nrow = n, ncol = num_loci)
for (i in seq_len(num_loci)) {
  genotypes[, i] <- rbinom(n, 2, p_pop_true[i])
}
# 这里的 PRS 一般不需要标准化，但如果您想保持代码结构，可以继续使用
genotypes = apply(genotypes, 2, scale) 

# 2. 基于 Logistic 模型生成疾病状态 (Y)
LP_genetic <- genotypes %*% true_betas
LP_total <- alpha_intercept + LP_genetic
P_i_true <- exp(LP_total) / (1 + exp(LP_total))

set.seed(456)
Y <- rbinom(n, 1, P_i_true)

cat(sprintf("人群模拟流行率: %f (理论值 exp(alpha)/(1+exp(alpha)) 约为: %f)\n",
            mean(Y), exp(alpha_intercept) / (1 + exp(alpha_intercept))))
cat(sprintf("人群中Case数: %d\n", sum(Y)))

# 3. 病例-对照抽样
set.seed(123)
case_indices <- which(Y == 1)
control_indices <- which(Y == 0)

if (length(case_indices) < N_cases) {
  stop("病例数量不足！请增大总样本量 n 或减小 alpha_intercept。")
}
if (length(control_indices) < N_controls) {
  stop("对照数量不足！请增大总样本量 n 或减小 alpha_intercept。")
}

# 4. 模拟 GWAS 并构建 PRS (重复 30 次独立模拟)
n_sims <- 30
all_roc_objects <- vector("list", n_sims)
all_sim_auc_roc_results <- numeric(n_sims)

for (sim_iter in seq_len(n_sims)) {
  set.seed(100 + sim_iter)

  sampled_case_indices <- sample(case_indices, N_cases)
  sampled_control_indices <- sample(control_indices, N_controls)

  genotypes_sampled <- rbind(
    genotypes[sampled_case_indices, ],
    genotypes[sampled_control_indices, ]
  )
  Y_sampled <- c(Y[sampled_case_indices], Y[sampled_control_indices])

  genotypes_sampled_scaled <- apply(genotypes_sampled, 2, scale)
  colnames(genotypes_sampled_scaled) <- paste0("SNP_", seq_len(num_loci))
  data_train_df <- data.frame(Y = Y_sampled, genotypes_sampled_scaled)

  # GWAS：单SNP回归，稳健取beta
  individual_snp_betas <- numeric(num_loci)
  for (i in seq_len(num_loci)) {
    snp_name <- paste0("SNP_", i)
    formula <- as.formula(paste("Y ~", snp_name))
    current_model <- tryCatch(
      lm(formula, data = data_train_df),
      error = function(e) NULL
    )
    if (is.null(current_model)) {
      individual_snp_betas[i] <- 0
    } else {
      coefs <- coef(current_model)
      if (length(coefs) >= 2 && !is.na(coefs[2])) {
        individual_snp_betas[i] <- as.numeric(coefs[2])
      } else {
        individual_snp_betas[i] <- 0
      }
    }
  }

  # 构建 PRS 并计算 AUC
  y_hat <- as.numeric(genotypes_sampled_scaled %*% individual_snp_betas)
  roc_curve <- tryCatch(
    roc(Y_sampled, y_hat, quiet = TRUE),
    error = function(e) NULL
  )
  all_roc_objects[[sim_iter]] <- roc_curve
  all_sim_auc_roc_results[sim_iter] <- if (!is.null(roc_curve)) auc(roc_curve) else NA_real_
}

# 5. 结果分析与绘图
T_simplified <- sqrt(num_loci * (N_controls + N_cases) /
                     (2 * N_cases * N_controls))
AUC_theory <- pnorm(T_simplified)

pdf("/data2/qiuyue/PRS_simulation/all_ROC_curves_100_100_100.pdf", width = 8, height = 6)
valid_rocs <- Filter(Negate(is.null), all_roc_objects)
if (length(valid_rocs) > 0) {
  plot(valid_rocs[[1]], main = "All ROC Curves for Simulated PRS",
       print.auc = FALSE, col = rgb(0.1, 0.1, 0.8, 0.2), lwd = 1,
       legacy.axes = TRUE)
  lines(c(0, 1), c(0, 1), col = "gray", lty = 2)

  if (length(valid_rocs) > 1) {
    for (i in seq(2, length(valid_rocs))) {
      lines(valid_rocs[[i]], col = rgb(0.1, 0.1, 0.8, 0.2), lwd = 1)
    }
  }

  mean_auc_roc <- mean(all_sim_auc_roc_results, na.rm = TRUE)
  sd_auc_roc <- sd(all_sim_auc_roc_results, na.rm = TRUE)

  text_label <- paste0("Mean AUC (PRS): ", round(mean_auc_roc, 4), "\n",
                       "AUC SD (PRS): ", round(sd_auc_roc, 4), "\n",
                       "Simplified Theoretical AUC: ", round(AUC_theory, 3), "\n")
  text(0.6, 0.3, text_label, col = "red", cex = 0.9, font = 2)
} else {
  message("没有生成任何ROC曲线对象，无法绘图。")
}
dev.off()

cat("\n--- Simulation Results Summary ---\n")
cat("AUC for each simulation: ", all_sim_auc_roc_results, "\n")
cat("Mean AUC: ", mean(all_sim_auc_roc_results, na.rm = TRUE), "\n")
cat("AUC SD: ", sd(all_sim_auc_roc_results, na.rm = TRUE), "\n")
cat("Simplified Theoretical AUC: ", AUC_theory, "\n")
