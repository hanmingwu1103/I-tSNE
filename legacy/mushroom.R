library(Rtsne)
library(ggplot2)
library(dplyr)
library(tidyr)
library(stringr)
library(RSDA)
library(mdsOpt)
library(smds)
library(reshape2)
library(dataSDA)

####畫圖######
plot_interval_projection <- function(data,
                                     title = "Interval Projection",
                                     xlab = "Dimension 1",
                                     ylab = "Dimension 2",
                                     base_size = 12,
                                     palette = c("auto", "hcl", "hue", "brewer_set2", "brewer_set3")) {
  palette <- match.arg(palette)
  library(ggplot2)
  
  # 確保 Group 是 factor，並取得類別數（排除 NA）
  if (!"Group" %in% names(data)) stop("data 必須包含 Group 欄位")
  data$Group <- factor(data$Group)
  lvls <- levels(data$Group)
  n_groups <- length(lvls)
  
  # 安全防呆：若沒有任何群組（或全 NA），給一個假群組避免色盤長度為 0
  if (n_groups == 0) {
    data$Group <- factor("All")
    lvls <- levels(data$Group)
    n_groups <- 1
  }
  
  # 依情況自動選色：>8 顏色時改用 hcl.colors，避免 brewer 的上限
  cols <- switch(palette,
                 auto = if (n_groups <= 8) RColorBrewer::brewer.pal(n_groups, "Set2")
                 else grDevices::hcl.colors(n_groups, "Dark 3"),
                 hcl  = grDevices::hcl.colors(n_groups, "Dark 3"),
                 hue  = scales::hue_pal()(n_groups),
                 brewer_set2 = {
                   if (n_groups > 8) warning("Set2 最多 8 色；已自動回退到 hcl.colors。")
                   if (n_groups <= 8) RColorBrewer::brewer.pal(n_groups, "Set2")
                   else grDevices::hcl.colors(n_groups, "Dark 3")
                 },
                 brewer_set3 = {
                   if (n_groups > 12) warning("Set3 最多 12 色；已自動回退到 hcl.colors。")
                   if (n_groups <= 12) RColorBrewer::brewer.pal(n_groups, "Set3")
                   else grDevices::hcl.colors(n_groups, "Dark 3")
                 }
  )
  
  ggplot(
    data,
    aes(xmin = Dim1_Lower, xmax = Dim1_Upper,
        ymin = Dim2_Lower, ymax = Dim2_Upper,
        colour = Group)
  ) +
    geom_rect(fill = NA, linewidth = 0.7, show.legend = TRUE) +
    coord_equal() +  # 固定比例，避免方框變形
    scale_colour_manual(values = setNames(cols, lvls), drop = FALSE) +
    labs(title = title, x = xlab, y = ylab, colour = "Group") +
    theme_classic(base_size = base_size) +
    theme(
      legend.position = "bottom",
      plot.title = element_text(hjust = 0.5, face = "bold")
    )
}

####facedata####
set.seed(12345)
facedata <- read.csv("C:/Users/may/OneDrive - National ChengChi University/meeting/data/facedata.csv",
                     fileEncoding = "UTF-8", stringsAsFactors = FALSE)

interval_df <- facedata[, c("species",
                            grep("_(Lower|Upper)$", names(facedata), value = TRUE))]

interval_df<-standardize_interval_data(interval_df,id_col = "species" ,method = 1)
vm_pca <- vm_pca_projection(interval_df, id_col = "species")

qm_pca <- qm_pca_projection(interval_df, m = 4, id_col = "species")
cr_pca <- cr_pca_projection(interval_df, id_col = "species")
imds <- run_imds(interval_df, id_col = "species")
vm_tsne <- vm_tsne_projection (interval_df,id_col = "species",
                               dims = 2,perplexity =50,theta = 0.5,
                               normalize = FALSE,max_iter = 1000,eta = 200)
qm_tsne <- qm_tsne_projection (interval_df, id_col = "species", m =5,dims = 2,
                               perplexity = 50,theta = 0.5,
                              
                               eta = 200,max_iter = 1000)

lU2_tsne<- lu1_tsne_ab_projection(data = interval_df,id_col = "species",dims = 2,
                                  perplexity = 5,alpha = 0.5, penalty_lambda =0.07,
                                  learning_rate =5,max_iter = 1000,
                                  initial_P_gain = 4,momentum = 0.5,
                                  final_momentum = 0.8,mom_switch_iter = 250,
                                  seed = 12345,verbose = TRUE)
plot_interval_projection(
  data =lU2_tsne$embedding,
  title = "I-tSNE(MM)",
  xlab = "Dim1",
  ylab = "Dim2"
)
cr1_tsne <- cr1_tsne_wass_projection(data = interval_df,id_col = "species",
                                dims = 2,perplexity = 5,lambda = 1.0,
                                eta = 5,max_iter = 1000,EE = 1.0,
                                T_EE = 250,momentum = 0.5,
                                final_momentum = 0.8,mom_switch_iter = 250,
                                verbose = TRUE)

plot_interval_pca_tsne_grid(
  vm_pca = vm_pca,
  qm_pca = qm_pca,
  cr_pca = cr_pca,
  imds   = imds,
  vm_tsne = vm_tsne,
  qm_tsne = qm_tsne,
  cr1_tsne_embedding = cr1_tsne$embedding,
  lU2_tsne_embedding = lU2_tsne$embedding,
  out_dir  = "C:/Users/may/OneDrive - National ChengChi University/meeting/論文latex/中文/Interval-tSNE_Wu&Wu_20260305/figure",                 # 你指定資料夾
  out_file = "facedata.pdf"   # 你指定英文檔名
)



methods_list <- list(vm_pca = vm_pca, qm_pca = qm_pca, cr_pca = cr_pca,
                     vm_tsne = vm_tsne, qm_tsne = qm_tsne, imds = imds, cr1_tsne = cr1_tsne$embedding,
                   lu2_tsne = lU2_tsne$embedding)

#LCMC
lcmc_tables <- compute_lcmc_tables(high_df = interval_df,
                                   methods_list = methods_list,
                                   k_range = 1:26,
                                   gamma = 0.5)
lcmc_w_df <- lcmc_tables$LCMC_wasserstein
lcmc_h_df <- lcmc_tables$LCMC_hausdorff
lcmc_m_df <- lcmc_tables$LCMC_minkowski
name_map <- c(
  vm_pca="IPCA(VM)",
  qm_pca ="IPCA(QM)",
  cr_pca="IPCA(CR)",
  imds="IMDS",
  cr1_tsne = "I-tSNE(CR)",
  qm_tsne  = "I-tSNE(QM)",
  lu2_tsne = "I-tSNE(MM)",
  vm_tsne  = "I-tSNE(VM)"
)
plot_lcmc_1x3(lcmc_w_df, lcmc_h_df, lcmc_m_df,
              out_dir  = "C:/Users/may/OneDrive - National ChengChi University/meeting/論文latex/中文/Interval-tSNE_Wu&Wu_20260305/figure",                 # 你指定資料夾
              out_file = "facedata_lcmc.pdf",width = 10,
              height = 6,
              save = TRUE,method_name_map = name_map)




####digit####
set.seed(12345)
digit <- read.csv("C:/Users/may/OneDrive - National ChengChi University/meeting/data/digits_interval_pca.csv",
                  fileEncoding = "UTF-8", stringsAsFactors = FALSE)

digit <- read.csv("C:/Users/may/OneDrive - National ChengChi University/meeting/data/digits_interval.csv",
                  fileEncoding = "UTF-8", stringsAsFactors = FALSE)

interval_df <- digit[, c("label",
                         grep("_(Lower|Upper)$", names(digit), value = TRUE))]


interval_df<-standardize_interval_data(interval_df,id_col = "label" ,method = 3)

vm_pca <- vm_pca_projection(interval_df, id_col = "label")

plot_interval_projection(
  data = vm_pca,
  title = "VMPCA",
  xlab = "Dim1",
  ylab = "Dim2"
)
#cm_pca <- cm_pca_projection(interval_df, id_col = "species")
qm_pca <- qm_pca_projection(interval_df, m = 5, id_col = "label")
plot_interval_projection(
  data = qm_pca,
  title = "QMPCA",
  xlab = "Dim1",
  ylab = "Dim2"
)
cr_pca <- cr_pca_projection(interval_df, id_col = "label")
plot_interval_projection(
  data = cr_pca,
  title = "IPCA(CR)",
  xlab = "Dim1",
  ylab = "Dim2"
)
imds <- run_imds(interval_df, id_col = "label")
plot_interval_projection(
  data = imds,
  title = "IMDS",
  xlab = "Dim1",
  ylab = "Dim2"
)


vm_tsne <- vm_tsne_projection (interval_df,id_col = "label",
                               dims = 2,perplexity =400,theta = 0.5,
                               normalize = FALSE,max_iter = 1000,eta = 200)
plot_interval_projection(
  data =vm_tsne,
  title = "I-tSNE(VM)",
  xlab = "Dim1",
  ylab = "Dim2"
)
qm_tsne <- qm_tsne_projection (interval_df, id_col = "label", m =5,dims = 2,
                               perplexity = 100,theta = 0.5,
                               
                               eta = 200,max_iter = 1000)
plot_interval_projection(
  data =qm_tsne,
  title = "I-tSNE(QM)",
  xlab = "Dim1",
  ylab = "Dim2"
)
lU2_tsne<- lu1_tsne_ab_projection(data = interval_df,id_col = "label",dims = 2,
                                  perplexity = 5,alpha = 0.5, penalty_lambda =0.4,
                                  learning_rate =200,max_iter = 1000,
                                  initial_P_gain = 12,momentum = 0.5,
                                  final_momentum = 0.8,mom_switch_iter = 250,
                                  seed = 12345,verbose = TRUE)
lU2_tsne<- lu1_tsne_ab_projection(data = interval_df,id_col = "label",dims = 2,
                                  perplexity = 5,alpha = 0.5, penalty_lambda =2,
                                  learning_rate =5,max_iter = 1000,
                                  initial_P_gain = 12,momentum = 0.5,
                                  final_momentum = 0.8,mom_switch_iter = 250,
                                  seed = 12345,verbose = TRUE)
lU2_tsne<- lu1_tsne_ab_projection(data = interval_df,id_col = "label",dims = 2,
                                  perplexity = 30,alpha = 0.5, penalty_lambda =1.8,
                                  learning_rate =5,max_iter = 1000,
                                  initial_P_gain = 4,momentum = 0.5,
                                  final_momentum = 0.8,mom_switch_iter = 250,
                                  seed = 12345,verbose = TRUE)


plot_interval_projection(
  data =lU2_tsne$embedding,
  title = "I-tSNE(MM)",
  xlab = "Dim1",
  ylab = "Dim2"
)
cr1_tsne <- cr1_tsne_wass_projection(data = interval_df,id_col = "label",
                                     dims = 2,perplexity = 50,lambda = 1.0,
                                     eta = 5,max_iter = 1000,EE = 12.0,
                                     T_EE = 250,momentum = 0.5,
                                     final_momentum = 0.8,mom_switch_iter = 250,
                                     verbose = TRUE)

plot_interval_projection(
  data =cr1_tsne$embedding,
  title = "I-tSNE(CR)",
  xlab = "Dim1",
  ylab = "Dim2"
)

plot_interval_pca_tsne_grid(
  vm_pca = vm_pca,
  qm_pca = qm_pca,
  cr_pca = cr_pca,
  imds   = imds,
  vm_tsne = vm_tsne,
  qm_tsne = qm_tsne,
  cr1_tsne_embedding = cr1_tsne$embedding,
  lU2_tsne_embedding = lU2_tsne$embedding,
  out_dir  = "C:/Users/may/OneDrive - National ChengChi University/meeting/論文latex/中文/Interval-tSNE_Wu&Wu_20260305/figure",                 # 你指定資料夾
  out_file = "digit.pdf"   # 你指定英文檔名
)


methods_list <- list(vm_pca = vm_pca, qm_pca = qm_pca, cr_pca = cr_pca,
                     qm_tsne = qm_tsne, imds = imds, vm_tsne = vm_tsne,cr1_tsne = cr1_tsne$embedding,
                     lu2_tsne = lU2_tsne$embedding)

#LCMC
lcmc_tables <- compute_lcmc_tables(high_df = interval_df,
                                   methods_list = methods_list,
                                   k_range = 1:100,
                                   gamma = 0.5)
lcmc_w_df <- lcmc_tables$LCMC_wasserstein
lcmc_h_df <- lcmc_tables$LCMC_hausdorff
lcmc_m_df <- lcmc_tables$LCMC_minkowski
name_map <- c(
  vm_pca='IPCA(VM)',
  qm_pca ="IPCA(QM)",
  cr_pca="IPCA(CR)",
  imda="IMDS",
  vm_tsne="I-tSNE(VM)",
  cr1_tsne = "I-tSNE(CR)",
  qm_tsne  = "I-tSNE(QM)",
  lu2_tsne = "I-tSNE(MM)"
)
plot_lcmc_1x3(lcmc_w_df, lcmc_h_df, lcmc_m_df,
              out_dir  = "C:/Users/may/OneDrive - National ChengChi University/meeting/論文latex/中文/Interval-tSNE_Wu&Wu_20260305/figure",                 # 你指定資料夾
              out_file = "diigit_lcmc.pdf",width = 10,
              height = 6,
              save = TRUE,method_name_map = name_map)




########irisdata########
iris<-read.csv("C:/Users/may/OneDrive - National ChengChi University/meeting/data/iris.csv",
               fileEncoding = "UTF-8", stringsAsFactors = FALSE)
interval_df <- iris[, c("Species",
                        grep("_(Lower|Upper)$", names(iris), value = TRUE))]
interval_df<-standardize_interval_data(interval_df,id_col = "Species", method = 1)
vm_pca <- vm_pca_projection(interval_df, id_col = "Species")
plot_interval_projection(
  data = vm_pca,
  title = "VMPCA",
  xlab = "Dim1",
  ylab = "Dim2"
)

qm_pca <- qm_pca_projection(interval_df, m = 4, id_col = "Species")
plot_interval_projection(
  data = qm_pca,
  title = "QMPCA",
  xlab = "Dim1",
  ylab = "Dim2"
)
cr_pca <- cr_pca_projection(interval_df, id_col = "Species")
plot_interval_projection(
  data = cr_pca,
  title = "CRPCA",
  xlab = "Dim1",
  ylab = "Dim2"
)
imds <- run_imds(interval_df, id_col = "Species")
plot_interval_projection(
  data = imds,
  title = "imds",
  xlab = "Dim1",
  ylab = "Dim2"
)
vm_tsne <- vm_tsne_projection (interval_df,id_col = "Species",
                               dims = 2,perplexity = 50,theta = 0.5,max_iter = 1000,eta = 200)
plot_interval_projection(
  data = vm_tsne,
  title = "It-sne(VM)",
  xlab = "Dim1",
  ylab = "Dim2"
)
qm_tsne <- qm_tsne_projection (interval_df, id_col = "Species", m =5,dims = 2,
                               perplexity = 50,theta = 0.5,
                               eta = 200,max_iter = 1000)
plot_interval_projection(
  data = qm_tsne,
  title = "It-sne(QM)",
  xlab = "Dim1",
  ylab = "Dim2"
)


lU2_tsne<- lu1_tsne_ab_projection(data = interval_df,id_col = "Species",dims = 2,
                                  perplexity = 5,alpha = 0.5, penalty_lambda =0.1,
                                  learning_rate =200,max_iter = 1000,
                                  initial_P_gain = 12,momentum = 0.5,
                                  final_momentum = 0.8,mom_switch_iter = 250,
                                  seed = 12345,verbose = TRUE)

plot_interval_projection(
  data = lU2_tsne$embedding,
  title = "It-sne(LU2)",
  xlab = "Dim1",
  ylab = "Dim2"
)



cr1_tsne <- cr1_tsne_wass_projection(data = interval_df,id_col = "Species",
                                     dims = 2,perplexity = 5,lambda = 1.0,
                                     eta = 10,max_iter = 1000,EE = 12.0,
                                     T_EE = 250,momentum = 0.5,
                                     final_momentum = 0.8,mom_switch_iter = 250,
                                     verbose = TRUE)

plot_interval_projection(
  data = cr1_tsne$embedding,
  title = "It-sne(CR1)",
  xlab = "Dim1",
  ylab = "Dim2"
)
plot_interval_pca_tsne_grid(
  vm_pca = vm_pca,
  qm_pca = qm_pca,
  cr_pca = cr_pca,
  imds   = imds,
  vm_tsne = vm_tsne,
  qm_tsne = qm_tsne,
  cr1_tsne_embedding = cr1_tsne$embedding,
  lU2_tsne_embedding = lU2_tsne$embedding,
  out_dir  = "C:/Users/may/OneDrive - National ChengChi University/meeting/論文latex/中文/Interval-tSNE_Wu&Wu_20260305/figure",                 # 你指定資料夾
  out_file = "iris_simu.pdf"   # 你指定英文檔名
)



methods_list <- list(vm_pca = vm_pca, qm_pca = qm_pca, cr_pca = cr_pca,
                     vm_tsne = vm_tsne, qm_tsne = qm_tsne, imds = imds, cr1_tsne = cr1_tsne$embedding,
                     lu2_tsne = lU2_tsne$embedding)

#LCMC
lcmc_tables <- compute_lcmc_tables(high_df = interval_df,
                                   methods_list = methods_list,
                                   k_range = 1:26,
                                   gamma = 0.5)
lcmc_w_df <- lcmc_tables$LCMC_wasserstein
lcmc_h_df <- lcmc_tables$LCMC_hausdorff
lcmc_m_df <- lcmc_tables$LCMC_minkowski
name_map <- c(
  vm_pca="IPCA(VM)",
  qm_pca ="IPCA(QM)",
  cr_pca="IPCA(CR)",
  imds="IMDS",
  cr1_tsne = "I-tSNE(CR)",
  qm_tsne  = "I-tSNE(QM)",
  lu2_tsne = "I-tSNE(MM)",
  vm_tsne  = "I-tSNE(VM)"
)
plot_lcmc_1x3(lcmc_w_df, lcmc_h_df, lcmc_m_df,
              out_dir  = "C:/Users/may/OneDrive - National ChengChi University/meeting/論文latex/中文/Interval-tSNE_Wu&Wu_20260305/figure",                 # 你指定資料夾
              out_file = "iris_simu_lcmc.pdf",width = 10,
              height = 6,
              save = TRUE,method_name_map = name_map)







#######Cars#####
car<-read.csv("C:/Users/may/OneDrive - National ChengChi University/meeting/data/Cars_data.csv",
              fileEncoding = "UTF-8", stringsAsFactors = FALSE)
interval_df <- car[, c("class",
                       grep("_(Lower|Upper)$", names(car), value = TRUE))]
interval_df<-standardize_interval_data(interval_df,id_col = "class", method = 1)

vm_pca <- vm_pca_projection(interval_df, id_col = "class")
plot_interval_projection(
  data = vm_pca,
  title = "VMPCA",
  xlab = "Dim1",
  ylab = "Dim2"
)


qm_pca <- qm_pca_projection(interval_df, m = 4, id_col = "class")
plot_interval_projection(
  data = qm_pca,
  title = "QMPCA",
  xlab = "Dim1",
  ylab = "Dim2"
)

cr_pca <- cr_pca_projection(interval_df, id_col = "class")
plot_interval_projection(
  data = cr_pca,
  title = "CRPCA",
  xlab = "Dim1",
  ylab = "Dim2"
)
imds <- run_imds(interval_df, id_col = "class")
plot_interval_projection(
  data = imds,
  title = "imds",
  xlab = "Dim1",
  ylab = "Dim2"
)
vm_tsne <- vm_tsne_projection (interval_df,id_col = "class",
                               dims = 2,perplexity =50,theta = 0.5,
                               normalize = FALSE,max_iter = 1000,eta = 200)
plot_interval_projection(
  data = vm_tsne,
  title = "It-sne(VM)",
  xlab = "Dim1",
  ylab = "Dim2"
)
qm_tsne <- qm_tsne_projection (interval_df, id_col = "class", m =4,dims = 2,
                               perplexity = 30,theta = 0.5,
                               
                               eta = 200,max_iter = 1000)
plot_interval_projection(
  data = qm_tsne,
  title = "It-sne(QM)",
  xlab = "Dim1",
  ylab = "Dim2"
)


lU2_tsne<- lu1_tsne_ab_projection(data = interval_df,id_col = "class",dims = 2,
                                  perplexity = 5,alpha = 0.5, penalty_lambda =0.07,
                                  learning_rate =200,max_iter = 1000,
                                  initial_P_gain = 12,momentum = 0.5,
                                  final_momentum = 0.8,mom_switch_iter = 250,
                                  seed = 12345,verbose = TRUE)
plot_interval_projection(
  data = lU2_tsne$embedding,
  title = "It-sne(LU2)",
  xlab = "Dim1",
  ylab = "Dim2"
)
cr1_tsne <- cr1_tsne_wass_projection(data = interval_df,id_col = "class",
                                     dims = 2,perplexity = 5,lambda = 1.0,
                                     eta = 5,max_iter = 1000,EE = 1.0,
                                     T_EE = 250,momentum = 0.5,
                                     final_momentum = 0.8,mom_switch_iter = 250,
                                     verbose = TRUE)



plot_interval_projection(
  data = cr1_tsne$embedding,
  title = "It-sne(CR1)",
  xlab = "Dim1",
  ylab = "Dim2"
)
plot_interval_pca_tsne_grid(
  vm_pca = vm_pca,
  qm_pca = qm_pca,
  cr_pca = cr_pca,
  imds   = imds,
  vm_tsne = vm_tsne,
  qm_tsne = qm_tsne,
  cr1_tsne_embedding = cr1_tsne$embedding,
  lU2_tsne_embedding = lU2_tsne$embedding,
  out_dir  = "C:/Users/may/OneDrive - National ChengChi University/meeting/論文latex/中文/Interval-tSNE_Wu&Wu_20260305/figure",                 # 你指定資料夾
  out_file = "car.pdf"   # 你指定英文檔名
)

methods_list <- list(vm_pca = vm_pca, qm_pca = qm_pca, cr_pca = cr_pca,
                     vm_tsne = vm_tsne, qm_tsne = qm_tsne, imds = imds, cr1_tsne = cr1_tsne$embedding,
                     lu2_tsne = lU2_tsne$embedding)

#LCMC
lcmc_tables <- compute_lcmc_tables(high_df = interval_df,
                                   methods_list = methods_list,
                                   k_range = 1:26,
                                   gamma = 0.5)
lcmc_w_df <- lcmc_tables$LCMC_wasserstein
lcmc_h_df <- lcmc_tables$LCMC_hausdorff
lcmc_m_df <- lcmc_tables$LCMC_minkowski
name_map <- c(
  vm_pca="IPCA(VM)",
  qm_pca ="IPCA(QM)",
  cr_pca="IPCA(CR)",
  imds="IMDS",
  cr1_tsne = "I-tSNE(CR)",
  qm_tsne  = "I-tSNE(QM)",
  lu2_tsne = "I-tSNE(MM)",
  vm_tsne  = "I-tSNE(VM)"
)
plot_lcmc_1x3(lcmc_w_df, lcmc_h_df, lcmc_m_df,
              out_dir  = "C:/Users/may/OneDrive - National ChengChi University/meeting/論文latex/中文/Interval-tSNE_Wu&Wu_20260305/figure",                 # 你指定資料夾
              out_file = "car_lcmc.pdf",width = 10,
              height = 6,
              save = TRUE,method_name_map = name_map)


######chinatemp####
china<-read.csv("C:/Users/may/OneDrive - National ChengChi University/meeting/data/ChinaTemp_data.csv",
                fileEncoding = "UTF-8", stringsAsFactors = FALSE)
interval_df <- china[, c("GeoReg",
                         grep("_(Lower|Upper)$", names(china), value = TRUE))]
interval_df<-standardize_interval_data(interval_df,id_col = "GeoReg" ,method = 1)
vm_pca <- vm_pca_projection(interval_df, id_col = "GeoReg")
plot_interval_projection(
  data = vm_pca,
  title = "VMPCA",
  xlab = "Dim1",
  ylab = "Dim2"
)


qm_pca <- qm_pca_projection(interval_df, m = 4, id_col = "GeoReg")
plot_interval_projection(
  data = qm_pca,
  title = "QMPCA",
  xlab = "Dim1",
  ylab = "Dim2"
)

cr_pca <- cr_pca_projection(interval_df, id_col = "GeoReg")
plot_interval_projection(
  data = cr_pca,
  title = "CRPCA",
  xlab = "Dim1",
  ylab = "Dim2"
)
imds <- run_imds(interval_df, id_col = "GeoReg")
plot_interval_projection(
  data = imds,
  title = "imds",
  xlab = "Dim1",
  ylab = "Dim2"
)
vm_tsne <- vm_tsne_projection (interval_df,id_col = "GeoReg",
                               dims = 2,perplexity =1000,theta = 0.5,
                               normalize = FALSE,max_iter = 1000,eta = 200)
plot_interval_projection(
  data = vm_tsne,
  title = "It-sne(VM)",
  xlab = "Dim1",
  ylab = "Dim2"
)
qm_tsne <- qm_tsne_projection (interval_df, id_col = "GeoReg", m =4,dims = 2,
                               perplexity = 100,theta = 0.5,
                           
                               eta = 200,max_iter = 1000)
plot_interval_projection(
  data = qm_tsne,
  title = "It-sne(QM)",
  xlab = "Dim1",
  ylab = "Dim2"
)



lU2_tsne<- lu1_tsne_ab_projection(data = interval_df,id_col = "GeoReg",dims = 2,
                                  perplexity = 200,alpha = 0.5, penalty_lambda = 2.4,
                                  learning_rate = 200,max_iter = 1000,
                                  initial_P_gain = 12,momentum = 0.5,
                                  final_momentum = 0.8,mom_switch_iter = 250,
                                  seed = 12345,verbose = TRUE)
plot_interval_projection(
  data = lU2_tsne$embedding,
  title = "It-sne(LU2)",
  xlab = "Dim1",
  ylab = "Dim2"
)
cr1_tsne <- cr1_tsne_wass_projection(data = interval_df,id_col = "GeoReg",
                                     dims = 2,perplexity = 200,lambda = 1.0,
                                     eta = 10,max_iter = 1000,EE = 12.0,
                                     T_EE = 250,momentum = 0.5,
                                     final_momentum = 0.8,mom_switch_iter = 250,
                                     verbose = TRUE)



plot_interval_projection(
  data = cr1_tsne$embedding,
  title = "It-sne(CR1)",
  xlab = "Dim1",
  ylab = "Dim2"
)
plot_interval_pca_tsne_grid(
  vm_pca = vm_pca,
  qm_pca = qm_pca,
  cr_pca = cr_pca,
  imds   = imds,
  vm_tsne = vm_tsne,
  qm_tsne = qm_tsne,
  cr1_tsne_embedding = cr1_tsne$embedding,
  lU2_tsne_embedding = lU2_tsne$embedding,
  out_dir  = "C:/Users/may/OneDrive - National ChengChi University/meeting/論文latex/中文/Interval-tSNE_Wu&Wu_20260305/figure",                 # 你指定資料夾
  out_file = "chinatemp.pdf"   # 你指定英文檔名
)



methods_list <- list(vm_pca = vm_pca, qm_pca = qm_pca, cr_pca = cr_pca,
                     vm_tsne = vm_tsne, qm_tsne = qm_tsne, imds = imds, cr1_tsne = cr1_tsne$embedding,
                     lu2_tsne = lU2_tsne$embedding)

#LCMC
lcmc_tables <- compute_lcmc_tables(high_df = interval_df,
                                   methods_list = methods_list,
                                   k_range = 1:400,
                                   gamma = 0.5)
lcmc_w_df <- lcmc_tables$LCMC_wasserstein
lcmc_h_df <- lcmc_tables$LCMC_hausdorff
lcmc_m_df <- lcmc_tables$LCMC_minkowski
name_map <- c(
  vm_pca="IPCA(VM)",
  qm_pca ="IPCA(QM)",
  cr_pca="IPCA(CR)",
  imda="IMDS",
  cr1_tsne = "I-tSNE(CR)",
  qm_tsne  = "I-tSNE(QM)",
  lu2_tsne = "I-tSNE(MM)",
  vm_tsne  = "I-tSNE(VM)"
)
plot_lcmc_1x3(lcmc_w_df, lcmc_h_df, lcmc_m_df,
              out_dir  = "C:/Users/may/OneDrive - National ChengChi University/meeting/論文latex/中文/Interval-tSNE_Wu&Wu_20260305/figure",                 # 你指定資料夾
              out_file = "chinatemp_lcmc.pdf",width = 10,
              height = 6,
              save = TRUE,method_name_map = name_map)














########mushroom#######
mushroom<-read.csv("C:/Users/may/OneDrive - National ChengChi University/meeting/data/mushroom_data.csv",
                   fileEncoding = "UTF-8", stringsAsFactors = FALSE)
interval_df <- mushroom[, c("Edibility",
                            grep("_(Lower|Upper)$", names(mushroom), value = TRUE))]
interval_df<-standardize_interval_data(interval_df,id_col = "Edibility", method = 1)
vm_pca <- vm_pca_projection(interval_df, id_col = "Edibility")
vm_pca$Species <- mushroom$Species
qm_pca <- qm_pca_projection(interval_df, m = 4, id_col = "Edibility")
cr_pca <- cr_pca_projection(interval_df, id_col = "Edibility")
imds <- run_imds(interval_df, id_col = "Edibility")

plot_interval_projection(
  data = vm_pca,
  title = "VMPCA",
  xlab = "Dim1",
  ylab = "Dim2"
)



plot_interval_projection(
  data = qm_pca,
  title = "QMPCA",
  xlab = "Dim1",
  ylab = "Dim2"
)


plot_interval_projection(
  data = cr_pca,
  title = "CRPCA",
  xlab = "Dim1",
  ylab = "Dim2"
)

plot_interval_projection(
  data = imds,
  title = "imds",
  xlab = "Dim1",
  ylab = "Dim2"
)
vm_tsne <- vm_tsne_projection (interval_df,id_col = "Edibility",
                               dims = 2,perplexity =50,theta = 0.5, normalize = FALSE,max_iter = 1000
                               ,eta = 200)
vm_tsne <- vm_tsne_projection (interval_df,id_col = "species",
                               dims = 2,perplexity =50,theta = 0.5,
                               normalize = FALSE,max_iter = 1000,eta = 200)
plot_interval_projection(
  data = vm_tsne,
  title = "It-sne(VM)",
  xlab = "Dim1",
  ylab = "Dim2"
)
qm_tsne <- qm_tsne_projection (interval_df, id_col = "Edibility", m =5,dims = 2,
                               perplexity = 40,theta = 0.5,
                               pca = FALSE,normalize = FALSE,
                               eta = 200,max_iter = 1000)
plot_interval_projection(
  data = qm_tsne,
  title = "It-sne(QM)",
  xlab = "Dim1",
  ylab = "Dim2"
)
lU1_tsne <- lu1_tsne(data =interval_df,id_col = "Edibility",
                     dims = 2,perplexity = 5, alpha = 0.5
                     ,learning_rate = 5,max_iter = 1000,
                     initial_P_gain = 1,momentum = 0.5,
                     final_momentum = 0.8,mom_switch_iter = 250,verbose = TRUE
)
plot_interval_projection(
  data = lU1_tsne$embedding,
  title = "It-sne(LU1)",
  xlab = "Dim1",
  ylab = "Dim2"
)
lU2_tsne<- lu1_tsne_ab_projection(data = interval_df,id_col = "Edibility",dims = 2,
                                  perplexity = 5,alpha = 0.5, penalty_lambda = 0.4,
                                  learning_rate = 200,max_iter = 1000,
                                  initial_P_gain = 12,momentum = 0.5,
                                  final_momentum = 0.8,mom_switch_iter = 250,
                                  seed = 12345,verbose = TRUE)
plot_interval_projection(
  data = lU2_tsne$embedding,
  title = "It-sne(LU2)",
  xlab = "Dim1",
  ylab = "Dim2"
)
cr1_tsne <- cr1_tsne_wass_projection(data = interval_df,id_col = "Edibility",
                                     dims = 2,perplexity = 5,lambda = 1.0,
                                     eta = 5,max_iter = 1000,EE = 4.0,
                                     T_EE = 250,momentum = 0.5,
                                     final_momentum = 0.8,mom_switch_iter = 250,
                                     verbose = TRUE)
plot_interval_projection(
  data = cr1_tsne$embedding,
  title = "It-sne(CR1)",
  xlab = "Dim1",
  ylab = "Dim2"
)
plot_interval_pca_tsne_grid(
  vm_pca = vm_pca,
  qm_pca = qm_pca,
  cr_pca = cr_pca,
  imds   = imds,
  vm_tsne = vm_tsne,
  qm_tsne = qm_tsne,
  cr1_tsne_embedding = cr1_tsne$embedding,
  lU2_tsne_embedding = lU2_tsne$embedding,
  out_dir  = "C:/Users/may/OneDrive - National ChengChi University/meeting/論文latex/中文/Interval-tSNE_Wu&Wu_20260126/figure",                 # 你指定資料夾
  out_file = "mashroom.pdf"   # 你指定英文檔名
)

methods_list <- list(vm_pca = vm_pca, cm_pca = cm_pca, qm_pca = qm_pca, cr_pca = cr_pca,
                     vm_tsne = vm_tsne, qm_tsne = qm_tsne, imds = imds, cr1_tsne = cr1_tsne$embedding,
                     lu1_tsne = lU1_tsne$embedding,lu2_tsne = lU2_tsne$embedding)

lcmc_tables <- compute_lcmc_tables(high_df = interval_df,
                                   methods_list = methods_list,
                                   k_range = 1:23,
                                   gamma = 0.5)
lcmc_w_df <- lcmc_tables$LCMC_wasserstein
lcmc_h_df <- lcmc_tables$LCMC_hausdorff
lcmc_m_df <- lcmc_tables$LCMC_minkowski
# 只顯示 1×3 圖
plot_lcmc_1x3(lcmc_w_df, lcmc_h_df, lcmc_m_df)




######Ablone######
abalone<-read.csv("C:/Users/may/OneDrive - National ChengChi University/meeting/data/Abalone_data.csv",
              fileEncoding = "UTF-8", stringsAsFactors = FALSE)

interval_df <- abalone[, c("subject",
                       grep("_(Lower|Upper)$", names(abalone), value = TRUE))]
interval_df<-standardize_interval_data(interval_df,id_col = "subject", method = 1)

vm_pca <- vm_pca_projection(interval_df, id_col = "subject")
plot_interval_projection(
  data = vm_pca,
  title = "VMPCA",
  xlab = "Dim1",
  ylab = "Dim2"
)


qm_pca <- qm_pca_projection(interval_df, m = 4, id_col = "subject")
plot_interval_projection(
  data = qm_pca,
  title = "QMPCA",
  xlab = "Dim1",
  ylab = "Dim2"
)

cr_pca <- cr_pca_projection(interval_df, id_col = "subject")
plot_interval_projection(
  data = cr_pca,
  title = "CRPCA",
  xlab = "Dim1",
  ylab = "Dim2"
)
imds <- run_imds(interval_df, id_col = "subject")
plot_interval_projection(
  data = imds,
  title = "imds",
  xlab = "Dim1",
  ylab = "Dim2"
)
vm_tsne <- vm_tsne_projection (interval_df,id_col = "subject",
                               dims = 2,perplexity =100,theta = 0.5,
                               normalize = FALSE,max_iter = 1000,eta = 200)
plot_interval_projection(
  data = vm_tsne,
  title = "It-sne(VM)",
  xlab = "Dim1",
  ylab = "Dim2"
)
qm_tsne <- qm_tsne_projection (interval_df, id_col = "subject", m =5,dims = 2,
                               perplexity = 40,theta = 0.5,
                               
                               eta = 200,max_iter = 1000)
plot_interval_projection(
  data = qm_tsne,
  title = "It-sne(QM)",
  xlab = "Dim1",
  ylab = "Dim2"
)

lU2_tsne<- lu1_tsne_ab_projection(data = interval_df,id_col = "subject",dims = 2,
                                  perplexity = 5,alpha = 0.5, penalty_lambda = 0.14,
                                  learning_rate = 200,max_iter = 1000,
                                  initial_P_gain = 12,momentum = 0.5,
                                  final_momentum = 0.8,mom_switch_iter = 250,
                                  seed = 12345,verbose = TRUE)
plot_interval_projection(
  data = lU2_tsne$embedding,
  title = "It-sne(LU2)",
  xlab = "Dim1",
  ylab = "Dim2"
)
cr1_tsne <- cr1_tsne_wass_projection(data = interval_df,id_col = "subject",
                                     dims = 2,perplexity = 5,lambda = 1.0,
                                     eta = 5,max_iter = 1000,EE = 12.0,
                                     T_EE = 250,momentum = 0.5,
                                     final_momentum = 0.8,mom_switch_iter = 250,
                                     verbose = TRUE)



plot_interval_projection(
  data = cr1_tsne$embedding,
  title = "It-sne(CR1)",
  xlab = "Dim1",
  ylab = "Dim2"
)
plot_interval_pca_tsne_grid(
  vm_pca = vm_pca,
  qm_pca = qm_pca,
  cr_pca = cr_pca,
  imds   = imds,
  vm_tsne = vm_tsne,
  qm_tsne = qm_tsne,
  cr1_tsne_embedding = cr1_tsne$embedding,
  lU2_tsne_embedding = lU2_tsne$embedding,
  out_dir  = "C:/Users/may/OneDrive - National ChengChi University/meeting/論文latex/中文/Interval-tSNE_Wu&Wu_20260305/figure",                 # 你指定資料夾
  out_file = "abalone.pdf"   # 你指定英文檔名
)

methods_list <- list(vm_pca = vm_pca, qm_pca = qm_pca, cr_pca = cr_pca,
                     vm_tsne = vm_tsne, qm_tsne = qm_tsne, imds = imds, cr1_tsne = cr1_tsne$embedding,
                     lu2_tsne = lU2_tsne$embedding)

#LCMC
lcmc_tables <- compute_lcmc_tables(high_df = interval_df,
                                   methods_list = methods_list,
                                   k_range = 1:26,
                                   gamma = 0.5)
lcmc_w_df <- lcmc_tables$LCMC_wasserstein
lcmc_h_df <- lcmc_tables$LCMC_hausdorff
lcmc_m_df <- lcmc_tables$LCMC_minkowski
name_map <- c(
  vm_pca="IPCA(VM)",
  qm_pca ="IPCA(QM)",
  cr_pca="IPCA(CR)",
  imds="IMDS",
  cr1_tsne = "I-tSNE(CR)",
  qm_tsne  = "I-tSNE(QM)",
  lu2_tsne = "I-tSNE(MM)",
  vm_tsne  = "I-tSNE(VM)"
)
plot_lcmc_1x3(lcmc_w_df, lcmc_h_df, lcmc_m_df,
              out_dir  = "C:/Users/may/OneDrive - National ChengChi University/meeting/論文latex/中文/Interval-tSNE_Wu&Wu_20260305/figure",                 # 你指定資料夾
              out_file = "abalone_lcmc.pdf",width = 10,
              height = 6,
              save = TRUE,method_name_map = name_map)


######wine######
wine<-read.csv("C:/Users/may/OneDrive - National ChengChi University/meeting/data/wine.csv",
                fileEncoding = "UTF-8", stringsAsFactors = FALSE)
interval_df <- wine[, c("Type",
                         grep("_(Lower|Upper)$", names(wine), value = TRUE))]
interval_df<-standardize_interval_data(interval_df,id_col = "Type", method = 1)
vm_pca <- vm_pca_projection(interval_df, id_col = "Type")
plot_interval_projection(
  data = vm_pca,
  title = "VMPCA",
  xlab = "Dim1",
  ylab = "Dim2"
)
cm_pca <- cm_pca_projection(interval_df, id_col = "Type")
plot_interval_projection(
  data = cm_pca,
  title = "CMPCA",
  xlab = "Dim1",
  ylab = "Dim2"
)

qm_pca <- qm_pca_projection(interval_df, m = 4, id_col = "Type")
plot_interval_projection(
  data = qm_pca,
  title = "QMPCA",
  xlab = "Dim1",
  ylab = "Dim2"
)

cr_pca <- cr_pca_projection(interval_df, id_col = "Type")
plot_interval_projection(
  data = cr_pca,
  title = "CRPCA",
  xlab = "Dim1",
  ylab = "Dim2"
)
imds <- run_imds(interval_df, id_col = "Type")
plot_interval_projection(
  data = imds,
  title = "imds",
  xlab = "Dim1",
  ylab = "Dim2"
)
vm_tsne <- vm_tsne_projection (interval_df,id_col = "Type",
                               dims = 2,perplexity =500,theta = 0.5,pca = FALSE,
                               normalize = FALSE,max_iter = 1000,eta = 200)
plot_interval_projection(
  data = vm_tsne,
  title = "It-sne(VM)",
  xlab = "Dim1",
  ylab = "Dim2"
)
qm_tsne <- qm_tsne_projection (interval_df, id_col = "Type", m =4,dims = 2,
                               perplexity = 40,theta = 0.5,
                               pca = FALSE,normalize = FALSE,
                               eta = 200,max_iter = 1000)
plot_interval_projection(
  data = qm_tsne,
  title = "It-sne(QM)",
  xlab = "Dim1",
  ylab = "Dim2"
)

library(ggplot2)
library(cowplot)


all_groups <- levels(factor(c(vm_pca$Group, cm_pca$Group, qm_pca$Group, cr_pca$Group)))
n_g <- length(all_groups)

cols <- setNames(grDevices::hcl.colors(n_g, palette = "Dark 3"), all_groups)


##########mnist########
mnist<-read.csv("C:/Users/may/OneDrive - National ChengChi University/meeting/data/mnist.csv",
               fileEncoding = "UTF-8", stringsAsFactors = FALSE)

interval_df <- mnist[, c("label",
                        grep("_(Lower|Upper)$", names(mnist), value = TRUE))]

vm_pca <- vm_pca_projection(interval_df, id_col = "label")
plot_interval_projection(
  data = vm_pca,
  title = "VMPCA",
  xlab = "Dim1",
  ylab = "Dim2"
)
cm_pca <- cm_pca_projection(interval_df, id_col = "Type")
plot_interval_projection(
  data = cm_pca,
  title = "CMPCA",
  xlab = "Dim1",
  ylab = "Dim2"
)

qm_pca <- qm_pca_projection(interval_df, m = 4, id_col = "Type")
plot_interval_projection(
  data = qm_pca,
  title = "QMPCA",
  xlab = "Dim1",
  ylab = "Dim2"
)

cr_pca <- cr_pca_projection(interval_df, id_col = "Type")
plot_interval_projection(
  data = cr_pca,
  title = "CRPCA",
  xlab = "Dim1",
  ylab = "Dim2"
)


























