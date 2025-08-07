#####此程序生成病例过程中不参考基因的影响，流行率决定函数截距
library(stats) # For pnorm (Phi function)
library(pROC)
#####需要先生成基因型 随后生成表型
num_loci <- 1000     # SNP 总数
min_maf <- 0.05      # 最小等位基因频率
max_maf <- 0.5       # 最大等位基因频率
n <-1000000 #总样本量
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
beta <- log(true_OR_derived)

#生成基因型数据
genotypes <- matrix(NA, nrow = n, ncol = num_loci)
for (i in 1:num_loci) {
  genotypes[, i] <- rbinom(n, 2, p_pop_true[i])
}
genotypes = apply(genotypes, 2, scale)  # 每列归一化

N_cases <- 1000    # 病例数
N_controls <- 2000    # 对照数

set.seed(123)
# 存储所有ROC曲线对象
all_roc_objects <- list()
all_sim_auc_roc_results <- numeric(30)  # 存储30次模拟的AUC值
for (sim_iter in 1:30) {  # 进行30次模拟
    
    # 确定每个样本是否患病
    alpha = -5  #The intercept \alpha was determined by the given prevalence (that is, the case-control ratios)
    p <- exp(alpha) / (1 + exp(alpha))  # 对于第i个病人，患病的概率时p
    
    p_i <- runif(n, min = 0, max = 1)
    # 根据p_i和p值判断每个样本是否患病
    y <- ifelse(p_i < p, 1, 0)  # 如果p_i大于p则为1（患病），否则为0（未患病）
    
    # 确认病例和对照的索引
    case_indices <- which(y == 1)
    control_indices <- which(y == 0)
    
    # 随机抽取指定数量的病例和对照索引
    sampled_case_indices <- sample(case_indices, N_cases)
    sampled_control_indices <- sample(control_indices, N_controls)
    
    # 提取对应的基因型和表型数据
    genotypes_cases <- genotypes[sampled_case_indices, ]
    genotypes_controls <- genotypes[sampled_control_indices, ]
    
    y_cases <- y[sampled_case_indices]
    y_controls <- y[sampled_control_indices]
    
    # 合并病例和对照数据
    genotypes_sampled <- rbind(genotypes_cases, genotypes_controls)
    y_sampled <- c(y_cases, y_controls)
    
    # 检查结果
    table(y_sampled)  # 应该显示 100 个病例和 100 个对照
    
    colnames(genotypes_sampled) <- paste0("Col_", 1:num_loci)  # 给每列命名
    data_train_df <- data.frame(Y = y_sampled, genotypes_sampled)
    
    individual_snp_models <- list()
    individual_snp_betas <- numeric(num_loci)
    for (i in 1:num_loci) {
      snp_i <- colnames(genotypes_sampled)[i]
      formula <- paste("Y ~", snp_i) #"Y ~ SNP1"
      current_model <- lm(as.formula(formula), data = data_train_df)
      individual_snp_models[[snp_i]] <- current_model
      individual_snp_betas[i] <- coef(current_model)[snp_i] # Get the coefficient for the SNP
    }
    
    # 计算y_hat
    y_hat <- genotypes_sampled %*% individual_snp_betas
    
    # 计算y_hat和y_sampled的AUC
    roc_curve <- roc(y_sampled, y_hat)
    all_roc_objects[[sim_iter]] <- roc_curve
    all_sim_auc_roc_results[sim_iter] <- auc(roc_curve)  # 存储每次模拟的AUC值
    
}

#lambda<-num_loci*(N_controls+N_cases)/(2*N_controls*N_cases)
T <- sqrt(num_loci*(N_controls+N_cases)/(2*N_cases*N_controls))
AUC_theory <- pnorm(T)

# 将ROC曲线绘制到PDF文件
pdf("/data2/qiuyue/PRS_simulation/all_ROC_curves.pdf", width = 8, height = 6)
# 绘制所有ROC曲线
if (length(all_roc_objects) > 0) {
  # 初始化绘图区域，绘制第一条ROC曲线
  plot(all_roc_objects[[1]], main = "all ROC curves", print.auc = FALSE, 
       col = rgb(0.1, 0.1, 0.8, 0.2), # 半透明蓝色
       lwd = 1, legacy.axes = TRUE)
  lines(c(0, 1), c(0, 1), col = "gray", lty = 2) # 添加对角线
  
  # 逐一添加剩余的ROC曲线
  for (i in 2:length(all_roc_objects)) {
    lines(all_roc_objects[[i]], col = rgb(0.1, 0.1, 0.8, 0.2), lwd = 1)
  }
  
  # （可选）计算并展示AUC的均值和标准差
  mean_auc_roc <- mean(all_sim_auc_roc_results)
  sd_auc_roc <- sd(all_sim_auc_roc_results)
  
  text_label <- paste0("mean AUC (ROC): ", round(mean_auc_roc, 4), "\n",
                       "AUC sd (ROC): ", round(sd_auc_roc, 4), "\n",
                       "theoritical AUC (ROC): ", round(AUC_theory, 3),"\n")
  text(0.6, 0.3, text_label, col = "red", cex = 0.9, font = 2)
} else {
  message("没有生成任何ROC曲线对象，无法绘图。")
}
dev.off()

# 输出每次模拟的AUC
cat("\n--- Simulation Results Summary ---\n")
cat("AUC for each simulation: ", all_sim_auc_roc_results, "\n")
cat("Mean AUC: ", mean(all_sim_auc_roc_results), "\n")
cat("AUC SD: ", sd(all_sim_auc_roc_results), "\n")

