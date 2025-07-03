m <- 5000
n1 = 1000 #training set size
n2 = 500 #testing set size
w = 0.25 #propertio of overlap
n_overlap = n2*w 
n <- n1 + n2 - n_overlap  #actual sample size
#set.seed(123)  # For reproducibility

# 每列的概率
p <- runif(m, min = 0.1, max = 0.9)

# 生成 n x m 的矩阵，每个元素是两个独立 Bernoulli(p_j) 的和
generate_row <- function(p) {
  rbinom(length(p), 1, p) + rbinom(length(p), 1, p)
}

# 使用 replicate 生成 n 行， row(x_matrix) 是 n， col(x_matrix) 是 m
x_matrix <- t(replicate(n, generate_row(p)))
cor(colMeans(x_matrix),p)

x_matrix = apply(x_matrix, 2, scale)  # 每列归一化  
# 查看结果维度
dim(x_matrix)  # 应该是 1000 x 10000
#print(x_matrix)

x_train<-x_matrix[1:n1,]
x_test<-x_matrix[(n1+1-n_overlap):n,]

beta<- rnorm(m, mean = 0, sd = 1)  # 生成 m 个 beta 系数

G <- x_matrix %*% beta
#print(G)
print(var(G))
#sum(beta * beta * 2 * p * (1 - p))
G<-0 #没有遗传效应

ei <- rnorm(n, mean = 0, sd = 1)
y <- G + ei
print(var(y))

y_test <- y[(n1+1-n_overlap):n]

colnames(x_train)
dim(x_train)
colnames(x_train) <- paste0("Col_", 1:m)  # 给每列命名
Y <- y[1:n1]  # y-training 
colnames(x_train)


data_train_df <- data.frame(Y = Y, x_train)
#head(data_train_df)

individual_snp_models <- list()
individual_snp_betas <- numeric(m)
for (i in 1:m) {
  snp_i <- colnames(x_train)[i]
  formula <- paste("Y ~", snp_i) #"Y ~ SNP1"
  current_model <- lm(as.formula(formula), data = data_train_df)
  individual_snp_models[[snp_i]] <- current_model
  individual_snp_betas[i] <- coef(current_model)[snp_i] # Get the coefficient for the SNP
  
#  cat(paste0("\nModel for ", snp_i, ":\n"))
#  print(summary(current_model)$coefficients) # Print summary of coefficients for each model
}

y_hat <- x_test %*% individual_snp_betas
var(y_hat)

#print(y_hat)
#print(y_test)

# 计算 y_hat 和 y 的相关性
corObs <- cor(y_hat, y_test, method = "pearson")
# 输出相关性
print(corObs)

predRsq= w^2*m/n1/(1+w*m/n1)
print(paste(predRsq,corObs^2))
