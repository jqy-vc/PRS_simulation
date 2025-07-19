library(stats) # For pnorm (Phi function)

# --- 模拟参数设置 ---
num_simulations <- 500 # 运行多少次模拟
N_cases <- 1000       # 病例数
N_controls <- 2000    # 对照数
num_loci <- 1000      # SNP 总数
min_maf <- 0.05       # 最小等位基因频率
max_maf <- 0.5       # 最大等位基因频率
train_ratio <- 1  # 训练集比例，如果为1.0则不分测试集
lambda<-num_loci*(N_controls+N_cases)/(2*N_controls*N_cases)

T <- sqrt(num_loci*(N_controls+N_cases)/(2*N_cases*N_controls))
print(T)
AUC_theory <- pnorm(T)


# --- 存储每次模拟的AUC结果的向量 ---
all_sim_auc_results <- numeric(num_simulations)
# --- 模拟循环开始 ---
for (s in 1:num_simulations) {
  # 1. 定义真实 SNP 参数
  # 总体等位基因频率 (作为对照组频率的基准)
  p_pop_true <- runif(num_loci, min_maf, max_maf)
  
  # p_cl_true 就是 p_pop_true
  p_cl_true <- p_pop_true
  q_cl_true <- 1 - p_cl_true
  
  # 所有SNP都将有效应
  p_cs_true <- numeric(num_loci) 

  for (i in 1:num_loci) { # 遍历所有 SNP
    p_cs_true[i] <- p_cl_true[i]
  }
  
  # 真实 OR 将从 p_cs_true 和 p_cl_true 派生
  q_cs_true <- 1 - p_cs_true
  true_OR_derived <- (p_cs_true * q_cl_true) / (p_cl_true * q_cs_true)
  
  # 2. 生成基因型数据 (病例和对照分开生成)
  genotypes_controls <- matrix(NA, nrow = N_controls, ncol = num_loci)
  for (i in 1:num_loci) {
    genotypes_controls[, i] <- rbinom(N_controls, 2, p_cl_true[i])
  }

  genotypes_cases <- matrix(NA, nrow = N_cases, ncol = num_loci)
  for (i in 1:num_loci) {
    genotypes_cases[, i] <- rbinom(N_cases, 2, p_cs_true[i])
  }

  X_combined <- rbind(genotypes_cases, genotypes_controls)
  Y_combined <- c(rep(1, N_cases), rep(0, N_controls))

  # 3. 划分训练集和测试集
  if (train_ratio < 1.0) {
    case_indices <- which(Y_combined == 1)
    control_indices <- which(Y_combined == 0)

    train_case_indices <- sample(case_indices, size = floor(N_cases * train_ratio))
    test_case_indices <- setdiff(case_indices, train_case_indices)

    train_control_indices <- sample(control_indices, size = floor(N_controls * train_ratio))
    test_control_indices <- setdiff(control_indices, train_control_indices)

    train_indices <- c(train_case_indices, train_control_indices)
    test_indices <- c(test_case_indices, test_control_indices)

    X_train <- X_combined[train_indices, , drop = FALSE]
    Y_train <- Y_combined[train_indices]

    X_test <- X_combined[test_indices, , drop = FALSE]
    Y_test <- Y_combined[test_indices]
    
    N_cases_train <- sum(Y_train == 1)
    N_controls_train <- sum(Y_train == 0)
  } else {
    X_train <- X_combined
    Y_train <- Y_combined
    X_test <- X_combined
    Y_test <- Y_combined
    
    N_cases_train <- N_cases
    N_controls_train <- N_controls
  }

  # --- 4. 从训练集估计效应 (beta_hat) ---
  X_train_cases <- X_train[Y_train == 1, , drop = FALSE]
  X_train_controls <- X_train[Y_train == 0, , drop = FALSE]

  p_cs_hat <- colMeans(X_train_cases) / 2
  q_cs_hat <- 1 - p_cs_hat

  p_cl_hat <- colMeans(X_train_controls) / 2
  q_cl_hat <- 1 - p_cl_hat

  p_cs_hat[p_cs_hat == 0] <- 1e-9
  p_cs_hat[p_cs_hat == 1] <- 1 - 1e-9
  q_cs_hat[q_cs_hat == 0] <- 1e-9
  q_cs_hat[q_cs_hat == 1] <- 1 - 1e-9

  p_cl_hat[p_cl_hat == 0] <- 1e-9
  p_cl_hat[p_cl_hat == 1] <- 1 - 1e-9
  q_cl_hat[q_cl_hat == 0] <- 1e-9
  q_cl_hat[q_cl_hat == 1] <- 1 - 1e-9

  OR_hat <- (p_cs_hat * q_cl_hat) / (p_cl_hat * q_cs_hat)
  beta_hat <- log(OR_hat)

  #beta_hat[is.infinite(beta_hat)] <- sign(beta_hat[is.infinite(beta_hat)]) * 1000
  beta_hat[is.nan(beta_hat)] <- 0

  # --- 5. 计算多基因风险评分 (PRS) ---
  PRS_test <- as.vector(X_test %*% beta_hat)

  # --- 6. 计算 Ds, T, AUC ---
  Ds <- sum(2 * (p_cs_hat - p_cl_hat) *beta_hat)
  print(Ds)
  variance_terms <- (1 / (2 * N_cases_train * p_cs_hat)) +
                        (1 / (2 * N_cases_train * q_cs_hat)) +
                        (1 / (2 * N_controls_train * p_cl_hat)) +
                        (1 / (2 * N_controls_train * q_cl_hat))
  #print(variance_terms)
  print(lambda*2)
  sigma_R2_cs <- sum(2 * p_cs_hat * q_cs_hat * variance_terms)
  print(sigma_R2_cs)  
  sigma_R2_cl <- sum(2 * p_cl_hat * q_cl_hat * variance_terms)
  print(sigma_R2_cl)
  T_val <- Ds / sqrt(sigma_R2_cl+sigma_R2_cs)
  print(T_val)
  
  AUC_val <- pnorm(T_val)
  # 将本次模拟的AUC结果记录下来
  all_sim_auc_results[s] <- AUC_val

  # 直接打印本次模拟的结果
  cat(paste0("Simulation ", s, " AUC: ", round(AUC_val, 4), "\n"))
}

# --- 计算并打印汇总统计 ---
cat("\n--- Simulation Results Summary ---\n")
print(paste("AUC_theory:", AUC_theory))
print(paste("lambda:", lambda))
print(paste("Average AUC across simulations:", mean(all_sim_auc_results)))
print(paste("Standard Deviation of AUC:", sd(all_sim_auc_results)))
