## NOTE (added for the public release):
## Absolute paths that pointed at a personal OneDrive folder have been replaced
## by the placeholder variable `LEGACY_DATA_DIR` below. Set it to a directory
## containing the input files before running anything in this script.
##
## This file is retained for historical reference only. It is NOT part of the
## supported `itsne` package, is excluded from the CRAN tarball, and is not
## covered by the package tests. Use the exported package functions instead.

LEGACY_DATA_DIR <- Sys.getenv("ITSNE_LEGACY_DATA_DIR", unset = ".")

library(mvtnorm)
library(cluster)
library(ggplot2)
library(Rtsne)
library(ggpubr) 
set.seed(12345)
##########################################
#               模擬驗證1  遠            #
#########################################
set.seed(12345)
n1<-300
n2<-300
n3<-300
C1 <- c( 6, 0)
C2 <- c(-6, 0)
C3 <- c( 0, 6)
# 定義共變異矩陣（單位矩陣）
Sigma <- diag(2)
group1_data <- rmvnorm(n1, mean = C1, sigma = Sigma)
group2_data <- rmvnorm(n2, mean = C2, sigma = Sigma)
group3_data <- rmvnorm(n3, mean = C3, sigma = Sigma)
# 合併數據並添加群標籤
true_data <- as.data.frame(rbind(group1_data, group2_data, group3_data))
names(true_data) <- c("X1", "X2")
true_data$Group <- factor(rep(c("G1", "G2", "G3"), c(n1, n2, n3)))

# 顯示真實數據的前幾行
head(true_data)

# 繪製散點圖並標上群標籤
plot(true_data$X1, true_data$X2, 
     col = true_data$Group,  # 使用 Group 欄位來決定顏色
     pch = 16,               # 實心圓點
     main = "真實數據的散點圖 (X1-X2)",
     xlab = "X1",
     ylab = "X2")

k<-10
total_intervals <- k* 3
interval_data_intra <- data.frame(IntervalID = 1:total_intervals)
current_id <- 1
# 步驟 2: 對每個真實群體（G1, G2, G3）進行獨立處理
for (g_name in levels(true_data$Group)) {
  
  # 找出屬於這個真實群體的點
  group_points <- true_data[true_data$Group == g_name, ]
  
  # 在該群體內部使用 k-means 進行分群
  kmeans_intra_group <- kmeans(group_points[, c("X1", "X2")], 
                               centers = k, 
                               nstart = 25)
  cluster_labels_intra <- kmeans_intra_group$cluster
  
  # 步驟 3: 聚合每個子群內的點，得到區間
  for (i in 1:k) {
    # 找出屬於這個子群的所有點
    sub_cluster_points <- group_points[cluster_labels_intra == i, ]
    
    # 聚合為 min-max 區間
    min_x1 <- min(sub_cluster_points$X1)
    max_x1 <- max(sub_cluster_points$X1)
    min_x2 <- min(sub_cluster_points$X2)
    max_x2 <- max(sub_cluster_points$X2)
    
    # 將結果儲存
    interval_data_intra[current_id, "X1_Lower"] <- min_x1
    interval_data_intra[current_id, "X1_Upper"] <- max_x1
    interval_data_intra[current_id, "X2_Lower"] <- min_x2
    interval_data_intra[current_id, "X2_Upper"] <- max_x2
    interval_data_intra[current_id, "TrueGroup"] <- g_name # 記錄原始群體標籤
    
    # 更新區間 ID
    current_id <- current_id + 1
  }
}

# 使用 ggplot2 繪製區間矩形


############################################
#        模擬1遠加入其他的干擾變數         #
#########################################
noise_dims<-3
noise_data <- matrix(rnorm(nrow(true_data) * noise_dims, mean = 0, sd = 1), 
                     nrow(true_data), 
                     noise_dims)
full_data <- cbind(true_data[, c("X1", "X2")], noise_data)
names(full_data) <- c("X1", "X2", paste0("X", 3:5))
full_data$Group<- true_data$Group

k_per_group <- 10

interval_data_noise <- data.frame() 

# 追蹤當前資料列
current_row <- 1

# 對每個真實群體（G1, G2, G3）進行獨立處理
for (g_name in levels(full_data$Group)) {
  
  # 找出屬於這個真實群體的點
  group_points <- full_data[full_data$Group == g_name, ]
  
  # 檢查群體點數是否足夠
  if (nrow(group_points) < k_per_group) {
    warning(paste("群體", g_name, "的點數少於", k_per_group, "，無法進行 k-means 聚合。"))
    next
  }
  
  # 在該群體內部使用 k-means 進行分群，只使用數值維度
  kmeans_intra_group <- kmeans(group_points[, !names(group_points) %in% "Group"],
                               centers = k_per_group,
                               nstart = 25)
  cluster_labels_intra <- kmeans_intra_group$cluster
  
  # 聚合每個子群內的點，得到區間
  for (i in 1:k_per_group) {
    # 找出屬於這個子群的所有高維點
    sub_cluster_points <- group_points[cluster_labels_intra == i, ]
    
    # 聚合為 min-max 區間
    min_vals <- apply(sub_cluster_points[, !names(sub_cluster_points) %in% "Group"], 2, min)
    max_vals <- apply(sub_cluster_points[, !names(sub_cluster_points) %in% "Group"], 2, max)
    
    # 建立一個新列
    new_row <- data.frame()
    
    # 將結果儲存到新列
    for (j in 1:5) {
      var_name <- paste0("X", j)
      new_row[1, paste0(var_name, "_Lower")] <- min_vals[j]
      new_row[1, paste0(var_name, "_Upper")] <- max_vals[j]
    }
    
    # 記錄原始群體標籤
    new_row[1, "TrueGroup"] <- g_name
    
    # 將新列合併到資料框中
    interval_data_noise <- rbind(interval_data_noise, new_row)  }
}

##########################################
#               模擬驗證1近              #
#########################################
set.seed(12345)
n1<-300
n2<-300
n3<-300
C1 <- c( 0, 1)
C2 <- c(-6, 0)
C3 <- c( 0, 6)
# 定義共變異矩陣（單位矩陣）
Sigma <- diag(2)
group1_data <- rmvnorm(n1, mean = C1, sigma = Sigma)
group2_data <- rmvnorm(n2, mean = C2, sigma = Sigma)
group3_data <- rmvnorm(n3, mean = C3, sigma = Sigma)
# 合併數據並添加群標籤
true_data <- as.data.frame(rbind(group1_data, group2_data, group3_data))
names(true_data) <- c("X1", "X2")
true_data$Group <- factor(rep(c("G1", "G2", "G3"), c(n1, n2, n3)))

# 顯示真實數據的前幾行
head(true_data)

# 繪製散點圖並標上群標籤
plot(true_data$X1, true_data$X2, 
     col = true_data$Group,  # 使用 Group 欄位來決定顏色
     pch = 16,               # 實心圓點
     main = "真實數據的散點圖 (X1-X2)",
     xlab = "X1",
     ylab = "X2")

k<-10
total_intervals <- k* 3
interval_data_intra_c <- data.frame(IntervalID = 1:total_intervals)
current_id <- 1
# 步驟 2: 對每個真實群體（G1, G2, G3）進行獨立處理
for (g_name in levels(true_data$Group)) {
  
  # 找出屬於這個真實群體的點
  group_points <- true_data[true_data$Group == g_name, ]
  
  # 在該群體內部使用 k-means 進行分群
  kmeans_intra_group <- kmeans(group_points[, c("X1", "X2")], 
                               centers = k, 
                               nstart = 25)
  cluster_labels_intra <- kmeans_intra_group$cluster
  
  # 步驟 3: 聚合每個子群內的點，得到區間
  for (i in 1:k) {
    # 找出屬於這個子群的所有點
    sub_cluster_points <- group_points[cluster_labels_intra == i, ]
    
    # 聚合為 min-max 區間
    min_x1 <- min(sub_cluster_points$X1)
    max_x1 <- max(sub_cluster_points$X1)
    min_x2 <- min(sub_cluster_points$X2)
    max_x2 <- max(sub_cluster_points$X2)
    
    # 將結果儲存
    interval_data_intra_c[current_id, "X1_Lower"] <- min_x1
    interval_data_intra_c[current_id, "X1_Upper"] <- max_x1
    interval_data_intra_c[current_id, "X2_Lower"] <- min_x2
    interval_data_intra_c[current_id, "X2_Upper"] <- max_x2
    interval_data_intra_c[current_id, "TrueGroup"] <- g_name # 記錄原始群體標籤
    
    # 更新區間 ID
    current_id <- current_id + 1
  }
}


#########################################
#        模擬1近加入其他的干擾變數      #
#########################################
noise_dims<-3
noise_data <- matrix(rnorm(nrow(true_data) * noise_dims, mean = 0, sd = 1), 
                     nrow(true_data), 
                     noise_dims)
full_data <- cbind(true_data[, c("X1", "X2")], noise_data)
names(full_data) <- c("X1", "X2", paste0("X", 3:5))
full_data$Group<- true_data$Group

k_per_group <- 10

interval_data_noise_c <- data.frame() 

# 追蹤當前資料列
current_row <- 1

# 對每個真實群體（G1, G2, G3）進行獨立處理
for (g_name in levels(full_data$Group)) {
  
  # 找出屬於這個真實群體的點
  group_points <- full_data[full_data$Group == g_name, ]
  
  # 檢查群體點數是否足夠
  if (nrow(group_points) < k_per_group) {
    warning(paste("群體", g_name, "的點數少於", k_per_group, "，無法進行 k-means 聚合。"))
    next
  }
  
  # 在該群體內部使用 k-means 進行分群，只使用數值維度
  kmeans_intra_group <- kmeans(group_points[, !names(group_points) %in% "Group"],
                               centers = k_per_group,
                               nstart = 25)
  cluster_labels_intra <- kmeans_intra_group$cluster
  
  # 聚合每個子群內的點，得到區間
  for (i in 1:k_per_group) {
    # 找出屬於這個子群的所有高維點
    sub_cluster_points <- group_points[cluster_labels_intra == i, ]
    
    # 聚合為 min-max 區間
    min_vals <- apply(sub_cluster_points[, !names(sub_cluster_points) %in% "Group"], 2, min)
    max_vals <- apply(sub_cluster_points[, !names(sub_cluster_points) %in% "Group"], 2, max)
    
    # 建立一個新列
    new_row <- data.frame()
    
    # 將結果儲存到新列
    for (j in 1:5) {
      var_name <- paste0("X", j)
      new_row[1, paste0(var_name, "_Lower")] <- min_vals[j]
      new_row[1, paste0(var_name, "_Upper")] <- max_vals[j]
    }
    
    # 記錄原始群體標籤
    new_row[1, "TrueGroup"] <- g_name
    
    # 將新列合併到資料框中
    interval_data_noise_c  <- rbind(interval_data_noise_c , new_row)
  }
}









###result_1模擬1遠2維######

interval_data_intra <- interval_data_intra[, names(interval_data_intra) != "IntervalID"]

vm_tsne_2 <- vm_tsne_projection (interval_data_intra ,id_col = "TrueGroup",
                                 dims = 2,perplexity =30,theta = 0.5,
                                 normalize = FALSE,max_iter = 1000,eta = 200)
qm_tsne_2 <- qm_tsne_projection (interval_data_intra, id_col = "TrueGroup", m =4,dims = 2,
                                 perplexity = 30,theta = 0.5,
                                 eta = 5,max_iter = 1000)
lU2_tsne_2<- lu1_tsne_ab_projection(data = interval_data_intra,id_col = "TrueGroup",dims = 2,
                                    perplexity = 5,alpha = 0.5, penalty_lambda = 0.2,
                                    learning_rate =200,max_iter = 1000,
                                    initial_P_gain = 12,momentum = 0.5,
                                    final_momentum = 0.8,mom_switch_iter = 250,
                                    seed = 12345,verbose = TRUE)
cr1_tsne_2 <- cr1_tsne_wass_projection(data = interval_data_intra,id_col = "TrueGroup",
                                       dims = 2,perplexity = 5,lambda = 1.0,
                                       eta = 5,max_iter = 1000,EE = 12.0,
                                       T_EE = 250,momentum = 0.1,
                                       final_momentum = 0.1,mom_switch_iter = 250,
                                       verbose = TRUE)



methods_list <- list(vm_tsne = vm_tsne_2, qm_tsne = qm_tsne_2, cr1_tsne = cr1_tsne_2$embedding,lu2_tsne = lU2_tsne_2$embedding)

#LCMC
lcmc_tables <- compute_lcmc_tables(high_df = interval_data_intra,
                                   methods_list = methods_list,
                                   k_range = 1:26,
                                   gamma = 0.5)
lcmc_w_df_2 <- lcmc_tables$LCMC_wasserstein
lcmc_h_df_2 <- lcmc_tables$LCMC_hausdorff
lcmc_m_df_2 <- lcmc_tables$LCMC_minkowski
name_map <- c(
  cr1_tsne = "I-tSNE(CR)",
  qm_tsne  = "I-tSNE(QM)",
  lu2_tsne = "I-tSNE(MM)",
  vm_tsne  = "I-tSNE(VM)"
)




###result_2模擬1遠noise######
set.seed(12345)
interval_data_noise <- interval_data_noise[, names(interval_data_noise) != "IntervalID"]

# 或者使用賦值為 NULL (但這樣會直接修改原始變數，除非您想保留原始變數，否則不推薦)
# interval_data_intra$IntervalID <- NULL
#interval_df<-standardize_interval_data(interval_data_intra,id_col = "TrueGroup" ,method = 1)

vm_tsne <- vm_tsne_projection (interval_data_noise ,id_col = "TrueGroup",
                               dims = 2,perplexity =300,theta = 0.5,
                               normalize = FALSE,max_iter = 1000,eta = 200)
qm_tsne <- qm_tsne_projection (interval_data_noise, id_col = "TrueGroup", m =5,dims = 2,
                               perplexity = 30,theta = 0.5, eta = 5,max_iter = 1000)


lU2_tsne<- lu1_tsne_ab_projection(data = interval_data_noise,id_col = "TrueGroup",dims = 2,
                                  perplexity = 5,alpha = 0.5, penalty_lambda = 0.1,
                                  learning_rate =200,max_iter = 1000,
                                  initial_P_gain = 12,momentum = 0.5,
                                  final_momentum = 0.8,mom_switch_iter = 250,
                                  seed = 12345,verbose = TRUE)
cr1_tsne <- cr1_tsne_wass_projection(data = interval_data_noise,id_col = "TrueGroup",
                                     dims = 2,perplexity = 5,lambda = 1.0,
                                     eta = 5,max_iter = 1000,EE = 12.0,
                                     T_EE = 250,momentum = 0.1,
                                     final_momentum = 0.1,mom_switch_iter = 250,
                                     verbose = TRUE)

methods_list <- list(vm_tsne = vm_tsne, qm_tsne = qm_tsne, cr1_tsne = cr1_tsne$embedding,lu2_tsne = lU2_tsne$embedding)

#LCMC
lcmc_tables <- compute_lcmc_tables(high_df = interval_data_noise,
                                   methods_list = methods_list,
                                   k_range = 1:26,
                                   gamma = 0.5)
lcmc_w_df <- lcmc_tables$LCMC_wasserstein
lcmc_h_df <- lcmc_tables$LCMC_hausdorff
lcmc_m_df <- lcmc_tables$LCMC_minkowski
name_map <- c(
  cr1_tsne = "I-tSNE(CR)",
  qm_tsne  = "I-tSNE(QM)",
  lu2_tsne = "I-tSNE(MM)",
  vm_tsne  = "I-tSNE(VM)"
)








##########結果視覺化圖##################
g1<-plot_interval_tsne_grid_2x4(
  vm_tsne_2d = vm_tsne_2,
  qm_tsne_2d = qm_tsne_2,
  lU2_tsne_2d = lU2_tsne_2$embedding,
  cr1_tsne_2d = cr1_tsne_2$embedding,
  vm_tsne_noise = vm_tsne,
  qm_tsne_noise = qm_tsne,
  lU2_tsne_noise = lU2_tsne$embedding,
  cr1_tsne_noise = cr1_tsne$embedding,
  
)
p1_lcmc <- plot_lcmc_2x3(
  lcmc_w_df_2, lcmc_h_df_2, lcmc_m_df_2,
  lcmc_w_df, lcmc_h_df, lcmc_m_df,
  method_name_map = name_map
)
###result_模擬1近2維######

interval_data_intra_c <- interval_data_intra_c[, names(interval_data_intra_c) != "IntervalID"]

vm_tsne_2 <- vm_tsne_projection (interval_data_intra_c ,id_col = "TrueGroup",
                                 dims = 2,perplexity =30,theta = 0.5,
                                 normalize = FALSE,max_iter = 1000,eta = 200)
qm_tsne_2 <- qm_tsne_projection (interval_data_intra_c, id_col = "TrueGroup", m =4,dims = 2,
                                 perplexity = 30,theta = 0.5,
                                 
                                 eta = 5,max_iter = 1000)
lU2_tsne_2<- lu1_tsne_ab_projection(data = interval_data_intra_c,id_col = "TrueGroup",dims = 2,
                                    perplexity = 5,alpha = 0.5, penalty_lambda = 0.2,
                                    learning_rate =200,max_iter = 1000,
                                    initial_P_gain = 12,momentum = 0.5,
                                    final_momentum = 0.8,mom_switch_iter = 250,
                                    seed = 12345,verbose = TRUE)
cr1_tsne_2 <- cr1_tsne_wass_projection(data = interval_data_intra_c,id_col = "TrueGroup",
                                       dims = 2,perplexity = 5,lambda = 1.0,
                                       eta = 5,max_iter = 1000,EE = 12.0,
                                       T_EE = 250,momentum = 0.1,
                                       final_momentum = 0.1,mom_switch_iter = 250,
                                       verbose = TRUE)

methods_list <- list(vm_tsne = vm_tsne_2, qm_tsne = qm_tsne_2, cr1_tsne = cr1_tsne_2$embedding,lu2_tsne = lU2_tsne_2$embedding)

#LCMC
lcmc_tables <- compute_lcmc_tables(high_df = interval_data_intra_c,
                                   methods_list = methods_list,
                                   k_range = 1:26,
                                   gamma = 0.5)
lcmc_w_df_2 <- lcmc_tables$LCMC_wasserstein
lcmc_h_df_2 <- lcmc_tables$LCMC_hausdorff
lcmc_m_df_2 <- lcmc_tables$LCMC_minkowski
name_map <- c(
  cr1_tsne = "I-tSNE(CR)",
  qm_tsne  = "I-tSNE(QM)",
  lu2_tsne = "I-tSNE(MM)",
  vm_tsne  = "I-tSNE(VM)"
)

###result_2模擬1近noise######

interval_data_noise_c <- interval_data_noise_c[, names(interval_data_noise_c) != "IntervalID"]


vm_tsne <- vm_tsne_projection (interval_data_noise_c ,id_col = "TrueGroup",
                               dims = 2,perplexity =200,theta = 0.5,
                               normalize = FALSE,max_iter = 1000,eta = 200)
qm_tsne <- qm_tsne_projection (interval_data_noise_c, id_col = "TrueGroup", m =9,dims = 2,
                               perplexity = 20,theta = 0.5,
                               eta = 5,max_iter = 1000)
lU2_tsne<- lu1_tsne_ab_projection(data = interval_data_noise_c,id_col = "TrueGroup",dims = 2,
                                  perplexity = 5,alpha = 0.5, penalty_lambda = 0.2,
                                  learning_rate =200,max_iter = 1000,
                                  initial_P_gain = 12,momentum = 0.5,
                                  final_momentum = 0.8,mom_switch_iter = 250,
                                  seed = 12345,verbose = TRUE)
cr1_tsne <- cr1_tsne_wass_projection(data = interval_data_noise_c,id_col = "TrueGroup",
                                     dims = 2,perplexity = 5,lambda = 1.0,
                                     eta = 5,max_iter = 1000,EE = 4.0,
                                     T_EE = 250,momentum = 0.1,
                                     final_momentum = 0.3,mom_switch_iter = 250,
                                     verbose = TRUE)



methods_list <- list(vm_tsne = vm_tsne, qm_tsne = qm_tsne, cr1_tsne = cr1_tsne$embedding,lu2_tsne = lU2_tsne$embedding)

#LCMC
lcmc_tables <- compute_lcmc_tables(high_df = interval_data_noise_c,
                                   methods_list = methods_list,
                                   k_range = 1:26,
                                   gamma = 0.5)
lcmc_w_df <- lcmc_tables$LCMC_wasserstein
lcmc_h_df <- lcmc_tables$LCMC_hausdorff
lcmc_m_df <- lcmc_tables$LCMC_minkowski
name_map <- c(
  cr1_tsne = "I-tSNE(CR)",
  qm_tsne  = "I-tSNE(QM)",
  lu2_tsne = "I-tSNE(MM)",
  vm_tsne  = "I-tSNE(VM)"
)

##################近的畫圖######################
g2<-plot_interval_tsne_grid_2x4(
  vm_tsne_2d = vm_tsne_2,
  qm_tsne_2d = qm_tsne_2,
  lU2_tsne_2d = lU2_tsne_2$embedding,
  cr1_tsne_2d = cr1_tsne_2$embedding,
  vm_tsne_noise = vm_tsne,
  qm_tsne_noise = qm_tsne,
  lU2_tsne_noise = lU2_tsne$embedding,
  cr1_tsne_noise = cr1_tsne$embedding
)
p2_lcmc <- plot_lcmc_2x3(
  lcmc_w_df_2, lcmc_h_df_2, lcmc_m_df_2,
  lcmc_w_df, lcmc_h_df, lcmc_m_df,
  method_name_map = name_map
)
library(gridExtra)
library(grid)
combined_plot <- arrangeGrob(
  g1,
  nullGrob(),
  g2,
  ncol = 3,
  widths = c(1, 0.08, 1)
)

grid.newpage()
grid.draw(combined_plot)


pdf(file.path(LEGACY_DATA_DIR, "sim1.pdf"), width = 10, height = 5)
grid.newpage()
grid.draw(combined_plot)
dev.off()
combined_plot <- arrangeGrob(
  p1_lcmc,
  nullGrob(),
  p2_lcmc,
  ncol = 3,
  widths = c(1, 0.08, 1)
)

grid.newpage()
grid.draw(combined_plot)


pdf(file.path(LEGACY_DATA_DIR, "Lcmc_sim1.pdf"), width = 10, height = 5)
grid.newpage()
grid.draw(combined_plot)
dev.off()



##########2維原始圖################
# 1. 確保座標範圍一致 (同時考慮兩個資料集的極值)
x_limits <- c(
  min(interval_data_intra$X1_Lower, interval_data_intra_c$X1_Lower, na.rm = TRUE), 
  max(interval_data_intra$X1_Upper, interval_data_intra_c$X1_Upper, na.rm = TRUE)
)
y_limits <- c(
  min(interval_data_intra$X2_Lower, interval_data_intra_c$X2_Lower, na.rm = TRUE), 
  max(interval_data_intra$X2_Upper, interval_data_intra_c$X2_Upper, na.rm = TRUE)
)

p1 <- ggplot(interval_data_intra, 
             aes(xmin = X1_Lower, xmax = X1_Upper,
                 ymin = X2_Lower, ymax = X2_Upper, 
                 color = TrueGroup)) +
  geom_rect(fill = NA, linewidth = 1) + 
  labs(x = "X1", y = "X2", color = "Group") +
  scale_x_continuous(limits = x_limits) + 
  scale_y_continuous(limits = y_limits) + 
  theme_minimal(base_size = 16) +
  theme(
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 16)
  ) +
  coord_fixed(ratio = 1)

# 3. 繪製 p2 (設定一樣保持原樣)
p2 <- ggplot(interval_data_intra_c, 
             aes(xmin = X1_Lower, xmax = X1_Upper,
                 ymin = X2_Lower, ymax = X2_Upper, 
                 color = TrueGroup)) +
  geom_rect(fill = NA, linewidth = 1) + 
  labs(x = "X1", y = "X2", color = "Group") +
  scale_x_continuous(limits = x_limits) + 
  scale_y_continuous(limits = y_limits) + 
  theme_minimal(base_size = 16) +
  theme(
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 16)
  ) +
  coord_fixed(ratio = 1)

# 4. 關鍵合併步驟：提取共用圖例
combined_plot <- ggarrange(
  p1, p2, 
  ncol = 2, nrow = 1,
  common.legend = TRUE,   # <--- 關鍵 1：自動提取並共用一個圖例
  legend = "right"        # <--- 關鍵 2：將提取出來的圖例統一放在最右側
)

# 顯示圖片
print(combined_plot)
ggsave(file.path(LEGACY_DATA_DIR, "simu1_2D_raw.pdf"), plot = combined_plot, width = 10, height = 4.5)





#################################
#          模擬驗證2            #
#            遠的2維            #
#################################
set.seed(12345)
n_per_cluster <- 10
n <- 3 * n_per_cluster
K_signal <- 2

sep <- 2.0
mu <- rbind(
  c(-sep,  0),
  c( 0 ,  sep),
  c( sep, -sep)
)

sigma_center <- 0.20
sigma_range  <- 0.10

# 區分主體群 (C1, C2, C3)
cluster_labels <- rep(1:3, each = n_per_cluster)

data_out <- data.frame(Group = cluster_labels)

# 建立儲存區間上下界的資料框
data_LU <- as.data.frame(matrix(NA, n, K_signal * 2))
names(data_LU) <- c("X1_Lower", "X2_Lower", "X1_Upper", "X2_Upper")

# 迴圈生成每一筆資料
for (i in 1:n) {
  g <- cluster_labels[i]
  
  # a) 生成中心點
  C_signal <- rnorm(K_signal, mu[g, ], sigma_center)
  
  # b) 生成半徑 r，並確保 r > 0
  R_signal <- pmax(0, rnorm(K_signal, mean = 0.5, sd = sigma_range))
  
  # c) 根據區間算術得到上下界
  lower_vec <- C_signal - R_signal
  upper_vec <- C_signal + R_signal
  
  # d) 存入資料框
  data_LU[i, paste0("X", 1:K_signal, "_Lower")] <- lower_vec
  data_LU[i, paste0("X", 1:K_signal, "_Upper")] <- upper_vec
}

# 最終資料集：遠的2維
sim_data_far_2d <- cbind(data_out, data_LU)
head(sim_data_far_2d)





#################################
#       模擬驗證2：遠的noise    #
#################################

# 參數設定
K_noise <- 3
K_dims  <- K_signal + K_noise

# 建立儲存雜訊區間上下界的資料框
noise_data_LU <- as.data.frame(matrix(NA, n, K_noise * 2))
names(noise_data_LU) <- c(
  paste0("X", (K_signal + 1):K_dims, "_Lower"),
  paste0("X", (K_signal + 1):K_dims, "_Upper")
)

# 迴圈生成雜訊維度
for (i in 1:n) {
  # a) 生成雜訊維度的中心點，服從 N(0,1)
  C_noise <- rnorm(K_noise, mean = 0, sd = 1)
  
  # b) 生成雜訊維度的半徑
  R_noise <- pmax(0, rnorm(K_noise, mean = 0.5, sd = 1))
  
  # c) 根據區間算術得到雜訊上下界
  lower_vec_noise <- C_noise - R_noise
  upper_vec_noise <- C_noise + R_noise
  
  # d) 存入資料框
  noise_data_LU[i, paste0("X", (K_signal + 1):K_dims, "_Lower")] <- lower_vec_noise
  noise_data_LU[i, paste0("X", (K_signal + 1):K_dims, "_Upper")] <- upper_vec_noise
}

# 最終資料集：遠的noise
sim_data_far_noise <- cbind(sim_data_far_2d, noise_data_LU)
head(sim_data_far_noise)

#################################
#          模擬驗證2           #
#            近的2維            #
#################################
set.seed(12345)
n_per_cluster <- 10
n <- 3 * n_per_cluster
K_signal <- 2

sep <- 1.5
mu <- rbind(
  c(-sep,  0),
  c( 0 ,  sep),
  c( 0 , -sep)
)

sigma_center <- 0.20
sigma_range  <- 0.10

# 區分主體群 (C1, C2, C3)
cluster_labels <- rep(1:3, each = n_per_cluster)

data_out <- data.frame(Group = cluster_labels)

# 建立儲存區間上下界的資料框
data_LU <- as.data.frame(matrix(NA, n, K_signal * 2))
names(data_LU) <- c("X1_Lower", "X2_Lower", "X1_Upper", "X2_Upper")

# 迴圈生成每一筆資料
for (i in 1:n) {
  g <- cluster_labels[i]
  
  # a) 生成中心點
  C_signal <- rnorm(K_signal, mu[g, ], sigma_center)
  
  # b) 生成半徑 r，並確保 r > 0
  R_signal <- pmax(0, rnorm(K_signal, mean = 0.5, sd = sigma_range))
  
  # c) 根據區間算術得到上下界
  lower_vec <- C_signal - R_signal
  upper_vec <- C_signal + R_signal
  
  # d) 存入資料框
  data_LU[i, paste0("X", 1:K_signal, "_Lower")] <- lower_vec
  data_LU[i, paste0("X", 1:K_signal, "_Upper")] <- upper_vec
}

# 最終資料集：近的2維
sim_data_near_2d <- cbind(data_out, data_LU)
head(sim_data_near_2d)



#################################
#       模擬驗證2：近的noise    #
#################################

# 參數設定
K_noise <- 3
K_dims  <- K_signal + K_noise

# 建立儲存雜訊區間上下界的資料框
noise_data_LU <- as.data.frame(matrix(NA, n, K_noise * 2))
names(noise_data_LU) <- c(
  paste0("X", (K_signal + 1):K_dims, "_Lower"),
  paste0("X", (K_signal + 1):K_dims, "_Upper")
)

# 迴圈生成雜訊維度
for (i in 1:n) {
  # a) 生成雜訊維度的中心點
  C_noise <- rnorm(K_noise, mean = 0, sd = 1)
  
  # b) 生成雜訊維度的半徑
  R_noise <- pmax(0, rnorm(K_noise, mean = 0.5, sd = 1))
  
  # c) 根據區間算術得到雜訊上下界
  lower_vec_noise <- C_noise - R_noise
  upper_vec_noise <- C_noise + R_noise
  
  # d) 存入資料框
  noise_data_LU[i, paste0("X", (K_signal + 1):K_dims, "_Lower")] <- lower_vec_noise
  noise_data_LU[i, paste0("X", (K_signal + 1):K_dims, "_Upper")] <- upper_vec_noise
}

# 最終資料集：近的noise
sim_data_near_noise <- cbind(sim_data_near_2d, noise_data_LU)
head(sim_data_near_noise)


###result_1######
lU2_tsne_2<- lu1_tsne_ab_projection(data =sim_data_far_2d,id_col = "Group",dims = 2,
                                    perplexity = 5,alpha = 0.5, penalty_lambda = 0.2,
                                    learning_rate =200,max_iter = 1000,
                                    initial_P_gain = 12,momentum = 0.5,
                                    final_momentum = 0.8,mom_switch_iter = 250,
                                    seed = 12345,verbose = TRUE)



vm_tsne_2 <- vm_tsne_projection (sim_data_far_2d,id_col = "Group",
                                 dims = 2,perplexity =30,theta = 0.5,
                                 normalize = FALSE,max_iter = 1000,eta = 200)

qm_tsne_2 <- qm_tsne_projection (sim_data_far_2d, id_col = "Group", m =5,dims = 2,
                                 perplexity = 40,theta = 0.5,
                                 eta = 200,max_iter = 1000)


cr1_tsne_2 <- cr1_tsne_wass_projection(data = sim_data_far_2d,id_col =  "Group",
                                       dims = 2,perplexity = 5,lambda = 1.0,
                                       eta = 5,max_iter = 1000,EE = 12.0,
                                       T_EE = 250,momentum = 0.5,
                                       final_momentum = 0.8,mom_switch_iter = 250,
                                       verbose = TRUE)


#LCMC
methods_list <- list(vm_tsne = vm_tsne_2, qm_tsne = qm_tsne_2, cr1_tsne = cr1_tsne_2$embedding,lu2_tsne = lU2_tsne_2$embedding)

lcmc_tables <- compute_lcmc_tables(high_df = sim_data_far_2d,
                                   methods_list = methods_list,
                                   k_range = 1:26,
                                   gamma = 0.5)
lcmc_w_df_2 <- lcmc_tables$LCMC_wasserstein
lcmc_h_df_2 <- lcmc_tables$LCMC_hausdorff
lcmc_m_df_2 <- lcmc_tables$LCMC_minkowski
name_map <- c(
  cr1_tsne = "I-tSNE(CR)",
  qm_tsne  = "I-tSNE(QM)",
  lu2_tsne = "I-tSNE(MM)",
  vm_tsne  = "I-tSNE(VM)"
)


###result_2######
lU2_tsne<- lu1_tsne_ab_projection(data = sim_data_far_noise,id_col = "Group",dims 
                                  = 2,
                                  perplexity = 5,alpha = 0.5, penalty_lambda = 0.1,
                                  learning_rate =200,max_iter = 1000,
                                  initial_P_gain = 12,momentum = 0.5,
                                  final_momentum = 0.8,mom_switch_iter = 250,
                                  seed = 12345,verbose = TRUE)



vm_tsne <- vm_tsne_projection (sim_data_far_noise,id_col = "Group",
                               dims = 2,perplexity =50,theta = 0.5,
                               normalize = FALSE,max_iter = 1000,eta = 200)

qm_tsne <- qm_tsne_projection (sim_data_far_noise, id_col = "Group", m =5,dims = 2,
                               perplexity = 30,theta = 0.5,
                               
                               eta = 200,max_iter = 1000)


cr1_tsne <- cr1_tsne_wass_projection(data = sim_data_far_noise,id_col =  "Group",
                                     dims = 2,perplexity = 5,lambda = 1.0,
                                     eta = 5,max_iter = 1000,EE = 12.0,
                                     T_EE = 250,momentum = 0.5,
                                     final_momentum = 0.8,mom_switch_iter = 250,
                                     verbose = TRUE)


#LCMC
methods_list <- list(vm_tsne = vm_tsne, qm_tsne = qm_tsne, cr1_tsne = cr1_tsne$embedding,lu2_tsne = lU2_tsne$embedding)

lcmc_tables <- compute_lcmc_tables(high_df = sim_data_far_noise,
                                   methods_list = methods_list,
                                   k_range = 1:26,
                                   gamma = 0.5)
lcmc_w_df <- lcmc_tables$LCMC_wasserstein
lcmc_h_df <- lcmc_tables$LCMC_hausdorff
lcmc_m_df <- lcmc_tables$LCMC_minkowski
name_map <- c(
  cr1_tsne = "I-tSNE(CR)",
  qm_tsne  = "I-tSNE(QM)",
  lu2_tsne = "I-tSNE(MM)",
  vm_tsne  = "I-tSNE(VM)"
)







#######遠的視覺化圖########################
g1<-plot_interval_tsne_grid_2x4(
  vm_tsne_2d = vm_tsne_2,
  qm_tsne_2d = qm_tsne_2,
  lU2_tsne_2d = lU2_tsne_2$embedding,
  cr1_tsne_2d = cr1_tsne_2$embedding,
  vm_tsne_noise = vm_tsne,
  qm_tsne_noise = qm_tsne,
  lU2_tsne_noise = lU2_tsne$embedding,
  cr1_tsne_noise = cr1_tsne$embedding,
  
)
p1_lcmc <- plot_lcmc_2x3(
  lcmc_w_df_2, lcmc_h_df_2, lcmc_m_df_2,
  lcmc_w_df, lcmc_h_df, lcmc_m_df,
  method_name_map = name_map
)
###result_3######


lU2_tsne_2<- lu1_tsne_ab_projection(data = sim_data_near_2d,id_col = "Group",dims = 2,
                                    perplexity = 5,alpha = 0.5, penalty_lambda = 0.15,
                                    learning_rate =200,max_iter = 1000,
                                    initial_P_gain = 12,momentum = 0.5,
                                    final_momentum = 0.8,mom_switch_iter = 250,
                                    seed = 12345,verbose = TRUE)


vm_tsne_2 <- vm_tsne_projection (sim_data_near_2d,id_col = "Group",
                                 dims = 2,perplexity =30,theta = 0.5,
                                 normalize = FALSE,max_iter = 1000,eta = 200)

qm_tsne_2 <- qm_tsne_projection (sim_data_near_2d, id_col = "Group", m =5,dims = 2,
                                 perplexity = 40,theta = 0.5,
                                 eta = 200,max_iter = 1000)



cr1_tsne_2 <- cr1_tsne_wass_projection(data = sim_data_near_2d,id_col =  "Group",
                                       dims = 2,perplexity = 5,lambda = 1.0,
                                       eta = 5,max_iter = 1000,EE = 12.0,
                                       T_EE = 250,momentum = 0.1,
                                       final_momentum = 0.1,mom_switch_iter = 250,
                                       verbose = TRUE)


methods_list <- list(vm_tsne = vm_tsne_2, qm_tsne = qm_tsne_2, cr1_tsne = cr1_tsne_2$embedding,lu2_tsne = lU2_tsne_2$embedding)

#LCMC
lcmc_tables <- compute_lcmc_tables(high_df = sim_data_near_2d,
                                   methods_list = methods_list,
                                   k_range = 1:26,
                                   gamma = 0.5)
lcmc_w_df_2 <- lcmc_tables$LCMC_wasserstein
lcmc_h_df_2 <- lcmc_tables$LCMC_hausdorff
lcmc_m_df_2 <- lcmc_tables$LCMC_minkowski
name_map <- c(
  cr1_tsne = "I-tSNE(CR)",
  qm_tsne  = "I-tSNE(QM)",
  lu2_tsne = "I-tSNE(MM)",
  vm_tsne  = "I-tSNE(VM)"
)


###result_4######

lU2_tsne<- lu1_tsne_ab_projection(data = sim_data_near_noise,id_col = "Group",dims 
                                  = 2,
                                  perplexity = 5,alpha = 0.5, penalty_lambda = 0.06,
                                  learning_rate =100,max_iter = 1000,
                                  initial_P_gain = 12,momentum = 0.5,
                                  final_momentum = 0.8,mom_switch_iter = 250,
                                  seed = 12345,verbose = TRUE)



vm_tsne <- vm_tsne_projection (sim_data_near_noise,id_col = "Group",
                               dims = 2,perplexity =300,theta = 0.5,
                               normalize = FALSE,max_iter = 1000,eta = 2000)

qm_tsne <- qm_tsne_projection (sim_data_near_noise, id_col = "Group", m =7,dims = 2,
                               perplexity = 50,theta = 0.5,
                               eta = 200,max_iter = 1000)



cr1_tsne <- cr1_tsne_wass_projection(data =sim_data_near_noise,id_col =  "Group",
                                     dims = 2,perplexity = 10,lambda = 1.0,
                                     eta = 5,max_iter = 1000,EE = 12.0,
                                     T_EE = 250,momentum = 0.5,
                                     final_momentum = 0.5,mom_switch_iter = 250,
                                     verbose = TRUE)



methods_list <- list(vm_tsne = vm_tsne, qm_tsne = qm_tsne, cr1_tsne = cr1_tsne$embedding,lu2_tsne = lU2_tsne$embedding)

#LCMC
lcmc_tables <- compute_lcmc_tables(high_df = sim_data_near_noise,
                                   methods_list = methods_list,
                                   k_range = 1:26,
                                   gamma = 0.5)
lcmc_w_df <- lcmc_tables$LCMC_wasserstein
lcmc_h_df <- lcmc_tables$LCMC_hausdorff
lcmc_m_df <- lcmc_tables$LCMC_minkowski
name_map <- c(
  cr1_tsne = "I-tSNE(CR)",
  qm_tsne  = "I-tSNE(QM)",
  lu2_tsne = "I-tSNE(MM)",
  vm_tsne  = "I-tSNE(VM)"
)


###########近的視覺化圖############
g2<-plot_interval_tsne_grid_2x4(
  vm_tsne_2d = vm_tsne_2,
  qm_tsne_2d = qm_tsne_2,
  lU2_tsne_2d = lU2_tsne_2$embedding,
  cr1_tsne_2d = cr1_tsne_2$embedding,
  vm_tsne_noise = vm_tsne,
  qm_tsne_noise = qm_tsne,
  lU2_tsne_noise = lU2_tsne$embedding,
  cr1_tsne_noise = cr1_tsne$embedding
)
p2_lcmc <- plot_lcmc_2x3(
  lcmc_w_df_2, lcmc_h_df_2, lcmc_m_df_2,
  lcmc_w_df, lcmc_h_df, lcmc_m_df,
  method_name_map = name_map
)

combined_plot <- arrangeGrob(
  g1,
  nullGrob(),
  g2,
  ncol = 3,
  widths = c(1, 0.08, 1)
)

grid.newpage()
grid.draw(combined_plot)


pdf(file.path(LEGACY_DATA_DIR, "sim2.pdf"), width = 10, height = 5)
grid.newpage()
grid.draw(combined_plot)
dev.off()
combined_plot <- arrangeGrob(
  p1_lcmc,
  nullGrob(),
  p2_lcmc,
  ncol = 3,
  widths = c(1, 0.08, 1)
)

grid.newpage()
grid.draw(combined_plot)


pdf(file.path(LEGACY_DATA_DIR, "Lcmc_sim2.pdf"), width = 10, height = 5)
grid.newpage()
grid.draw(combined_plot)
dev.off()
#########################畫2d圖######################

# 1. 確保座標範圍一致 (同時考慮遠的2維、近的2維兩個資料集的極值)
x_limits <- c(
  min(sim_data_far_2d$X1_Lower, sim_data_near_2d$X1_Lower, na.rm = TRUE), 
  max(sim_data_far_2d$X1_Upper, sim_data_near_2d$X1_Upper, na.rm = TRUE)
)

y_limits <- c(
  min(sim_data_far_2d$X2_Lower, sim_data_near_2d$X2_Lower, na.rm = TRUE), 
  max(sim_data_far_2d$X2_Upper, sim_data_near_2d$X2_Upper, na.rm = TRUE)
)

# 2. Group 轉成字串
sim_data_far_2d$Group  <- as.character(sim_data_far_2d$Group)
sim_data_near_2d$Group <- as.character(sim_data_near_2d$Group)

# 3. 繪製遠的2維
p1 <- ggplot(
  sim_data_far_2d,
  aes(
    xmin = X1_Lower, xmax = X1_Upper,
    ymin = X2_Lower, ymax = X2_Upper,
    color = Group
  )
) +
  geom_rect(fill = NA, linewidth = 1) +
  labs(x = "X1", y = "X2", color = "Group") +
  scale_x_continuous(limits = x_limits) +
  scale_y_continuous(limits = y_limits) +
  theme_minimal(base_size = 16) +
  theme(
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 16)
  ) +
  coord_fixed(ratio = 1)

# 4. 繪製近的2維
p2 <- ggplot(
  sim_data_near_2d,
  aes(
    xmin = X1_Lower, xmax = X1_Upper,
    ymin = X2_Lower, ymax = X2_Upper,
    color = Group
  )
) +
  geom_rect(fill = NA, linewidth = 1) +
  labs(x = "X1", y = "X2", color = "Group") +
  scale_x_continuous(limits = x_limits) +
  scale_y_continuous(limits = y_limits) +
  theme_minimal(base_size = 16) +
  theme(
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 16)
  ) +
  coord_fixed(ratio = 1)

# 5. 合併圖形，共用圖例
combined_plot <- ggarrange(
  p1, p2,
  ncol = 2, nrow = 1,
  common.legend = TRUE,
  legend = "right"
)

combined_plot

ggsave(file.path(LEGACY_DATA_DIR, "simu2_2D_raw.pdf"), plot = combined_plot, width = 10, height = 4.5)
