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
#################################
#1997 vmpca                     #
#################################
vm_pca_projection <- function(interval_df, id_col = "Subject", group_prefix_len = 3) {
  
  data_matrix <- interval_df[, !(names(interval_df) %in% id_col), drop = FALSE]
  var_prefixes <- unique(gsub("_(Lower|Upper)$", "", names(data_matrix)))
  
  # 1. 分離 lower / upper
  lower_mat <- as.matrix(data_matrix[, paste0(var_prefixes, "_Lower"), drop = FALSE])
  upper_mat <- as.matrix(data_matrix[, paste0(var_prefixes, "_Upper"), drop = FALSE])
  
  rownames(lower_mat) <- interval_df[[id_col]]
  rownames(upper_mat) <- interval_df[[id_col]]
  
  n <- nrow(lower_mat)
  p <- ncol(lower_mat)
  
  # 2. 生成每筆資料的所有 vertices
  generate_vertices <- function(a, b) {
    expand.grid(lapply(seq_along(a), function(i) c(a[i], b[i]))) |> as.matrix()
  }
  
  all_vertices <- vector("list", n)
  vertex_owner <- vector("list", n)
  vertex_weights <- vector("list", n)
  
  # 1997 paper 的等權版本：
  # 每個 object 權重 p_i = 1/n
  # 每個 vertex 權重 q_ik = p_i / 2^p
  object_weights <- rep(1 / n, n)
  
  for (i in seq_len(n)) {
    verts <- generate_vertices(lower_mat[i, ], upper_mat[i, ])
    q_i <- nrow(verts)  # 理論上是 2^p
    
    all_vertices[[i]] <- verts
    vertex_owner[[i]] <- rep(i, q_i)
    vertex_weights[[i]] <- rep(object_weights[i] / q_i, q_i)
  }
  
  vertices_matrix <- do.call(rbind, all_vertices)
  owner_index <- unlist(vertex_owner)
  w <- unlist(vertex_weights)
  
  colnames(vertices_matrix) <- var_prefixes
  
  # 3. 依 1997 Sommets method 概念：
  #    先用 vertex weights 算 mean，再中心化，再做 covariance PCA
  vertex_mean <- colSums(vertices_matrix * w)
  centered_vertices <- sweep(vertices_matrix, 2, vertex_mean, FUN = "-")
  
  # weighted covariance matrix
  cov_mat <- t(centered_vertices) %*% (centered_vertices * w)
  
  eig <- eigen(cov_mat, symmetric = TRUE)
  loadings <- eig$vectors[, 1:2, drop = FALSE]
  
  # 4. 主成分分數：用中心化後的 vertices 去投影
  proj_vertices <- centered_vertices %*% loadings
  
  # 5. 將每筆資料對應的 vertices 範圍轉成投影區間
  proj_lower <- matrix(NA_real_, n, 2)
  proj_upper <- matrix(NA_real_, n, 2)
  
  for (i in seq_len(n)) {
    idx <- which(owner_index == i)
    proj_i <- proj_vertices[idx, , drop = FALSE]
    proj_lower[i, ] <- apply(proj_i, 2, min)
    proj_upper[i, ] <- apply(proj_i, 2, max)
  }
  
  # 6. 組成資料框（輸出格式和你原本一樣）
  result <- data.frame(
    Group = interval_df[[id_col]],
    Dim1_Lower = proj_lower[, 1],
    Dim1_Upper = proj_upper[, 1],
    Dim2_Lower = proj_lower[, 2],
    Dim2_Upper = proj_upper[, 2]
  )
  
  result$Dim1 <- (result$Dim1_Lower + result$Dim1_Upper) / 2
  result$Dim2 <- (result$Dim2_Lower + result$Dim2_Upper) / 2
  
  return(result)
}

#########################################
#  qm修改                               #
#########################################
qm_pca_projection <- function(data, m = 4, id_col = "Subject") {
  
  df <- as.data.frame(data, stringsAsFactors = FALSE)
  if (!id_col %in% names(df)) stop("找不到識別欄：", id_col)
  
  # ===== 1) 找區間變數前綴 =====
  var_names <- setdiff(names(df), id_col)
  prefixes  <- unique(gsub("_(Lower|Upper)$", "", var_names))
  n <- nrow(df)
  p <- length(prefixes)
  if (p < 2) stop("區間變數數量不足（至少要 2 個變數）")
  
  lower_mat <- sapply(prefixes, function(pf) df[[paste0(pf, "_Lower")]], simplify = "matrix")
  upper_mat <- sapply(prefixes, function(pf) df[[paste0(pf, "_Upper")]], simplify = "matrix")
  colnames(lower_mat) <- colnames(upper_mat) <- prefixes
  
  # ===== 2) Quantile 展開（Ichino：均勻分布假設；含端點 0..m）=====
  expand_quantiles <- function(lower, upper, m) {
    sapply(0:m, function(k) lower + (upper - lower) * (k / m))
  }
  
  expanded_list <- lapply(seq_len(n), function(i) {
    t(expand_quantiles(lower_mat[i, ], upper_mat[i, ], m))  # (m+1) x p
  })
  X <- do.call(rbind, expanded_list)  # {n*(m+1)} x p
  
  # ===== 3) 文獻：PCA on correlation matrix（固定 Pearson）=====
  S <- cor(X, method = "pearson", use = "pairwise.complete.obs")
  S[is.na(S)] <- 0
  diag(S) <- 1
  S <- (S + t(S)) / 2
  
  eig <- eigen(S, symmetric = TRUE)
  V <- eig$vectors[, 1:2, drop = FALSE]   # loadings (p x 2)
  
  # ===== 4) scores：用展開後 z-score（與 correlation PCA 等價）=====
  # 注意：即使你的區間在外面已標準化，Ichino 的 cor-based PCA 在這一步仍會
  # 對「展開後矩陣 X」再做變數層級的 z-score（這是 correlation matrix 的本質）。
  Z <- scale(X, center = TRUE, scale = TRUE)
  coords <- Z %*% V                        # (n*(m+1)) x 2
  
  # ===== 5) 每筆資料取 (m+1) 點在 PC1/PC2 的 min/max（你的視覺化不變）=====
  dim1_bounds <- matrix(coords[, 1], nrow = m + 1, ncol = n)
  dim2_bounds <- matrix(coords[, 2], nrow = m + 1, ncol = n)
  
  dim1_lower <- apply(dim1_bounds, 2, min)
  dim1_upper <- apply(dim1_bounds, 2, max)
  dim2_lower <- apply(dim2_bounds, 2, min)
  dim2_upper <- apply(dim2_bounds, 2, max)
  
  out <- data.frame(
    ID = df[[id_col]],
    Group = df[[id_col]],   # 你原本就是 Group=ID
    Dim1_Lower = dim1_lower,
    Dim1_Upper = dim1_upper,
    Dim2_Lower = dim2_lower,
    Dim2_Upper = dim2_upper,
    check.names = FALSE
  )
  names(out)[1] <- id_col
  
  # attrs（可留可刪；不影響視覺化）
  attr(out, "m") <- m
  attr(out, "corr_matrix_S") <- S
  attr(out, "loadings") <- V
  attr(out, "eigenvalues") <- eig$values[1:2]
  
  return(out)
}



#########################################
#    CRPCA修改                          #
#########################################
cr_pca_projection <- function(data, id_col = "Subject", k = 2, eps = 1e-12) {
  
  # ===== 1) variable prefixes =====
  var_names <- names(data)[names(data) != id_col]
  prefixes <- unique(gsub("_(Lower|Upper)$", "", var_names))
  
  # ===== 2) lower/upper -> midpoints/radii =====
  lower <- sapply(prefixes, function(p) data[[paste0(p, "_Lower")]])
  upper <- sapply(prefixes, function(p) data[[paste0(p, "_Upper")]])
  colnames(lower) <- prefixes
  colnames(upper) <- prefixes
  
  Xc <- (lower + upper) / 2  # midpoints
  Xr <- (upper - lower) / 2  # radii
  N  <- nrow(Xc)
  p  <- ncol(Xc)
  
  if (k > p) stop("k cannot exceed the number of variables (p).")
  
  # ===== 3) Centering (mean interval centering) =====
  Xc <- scale(Xc, center = TRUE, scale = FALSE)
  Xr <- scale(Xr, center = TRUE, scale = FALSE)
  
  # ===== 4) Build Vx and Sigma from diag(Vx), then standardize =====
  # NOTE: Palumbo & Lauro (2003) use 1/N; we follow that to match the paper
  div <- N
  
  Cc  <- (t(Xc) %*% Xc) / div
  Cr  <- (t(Xr) %*% Xr) / div
  Ccr <- (t(Xc) %*% Xr) / div
  
  # Vx = Cc + Cr + |Ccr| + |Ccr'|
  Vx <- Cc + Cr + abs(Ccr) + abs(t(Ccr))
  
  sigma <- sqrt(pmax(diag(Vx), eps))  # avoid division by 0
  Zc <- sweep(Xc, 2, sigma, "/")
  Zr <- sweep(Xr, 2, sigma, "/")
  
  # ===== 5) Separate PCAs for midpoints and radii =====
  Rc <- (t(Zc) %*% Zc) / div
  Rr <- (t(Zr) %*% Zr) / div
  
  eig_c <- eigen(Rc, symmetric = TRUE)
  eig_r <- eigen(Rr, symmetric = TRUE)
  
  Uc <- eig_c$vectors[, 1:k, drop = FALSE]   # midpoint loadings
  Ur <- eig_r$vectors[, 1:k, drop = FALSE]   # radius loadings
  lambdas_c <- eig_c$values[1:k]
  lambdas_r <- eig_r$values[1:k]
  
  # Scores in their own PC spaces
  Score_c <- Zc %*% Uc  # N x k
  Score_r <- Zr %*% Ur  # N x k
  
  # ===== 6) Procrustes rotation in VARIABLE space (p x p), EXACTLY matching the paper =====
  # Paper: SVD of (Zc' Zr) = P Delta Q', then A = Q P'
  M_full  <- t(Zc) %*% Zr          # p x p  (MATCH paper: Xc'Xr)
  sv_full <- svd(M_full)
  P <- sv_full$u
  Q <- sv_full$v
  A_full  <- Q %*% t(P)            # A = Q P'
  
  # rotate radii matrix, then represent radii on MIDPOINT PCs (supplementary points)
  Zr_rot      <- Zr %*% A_full     # N x p
  Score_r_rot <- Zr_rot %*% Uc     # N x k
  
  # ===== 7) Interval reconstruction on PCs =====
  Y_lower <- Score_c - abs(Score_r_rot)
  Y_upper <- Score_c + abs(Score_r_rot)
  
  # ===== 8) Output =====
  out <- data.frame(Group = data[[id_col]])
  
  for (d in 1:k) {
    out[[paste0("Dim", d)]]            <- Score_c[, d]
    out[[paste0("Radius", d)]]         <- abs(Score_r_rot[, d])
    out[[paste0("Dim", d, "_Lower")]]  <- Y_lower[, d]
    out[[paste0("Dim", d, "_Upper")]]  <- Y_upper[, d]
  }
  
  # attrs for reproducibility/debug
  attr(out, "sigma") <- sigma
  attr(out, "Rc_eigenvalues") <- lambdas_c
  attr(out, "Rr_eigenvalues") <- lambdas_r
  attr(out, "Uc_loadings") <- Uc
  attr(out, "Ur_loadings") <- Ur
  attr(out, "rotation_A") <- A_full
  
  return(out)
}


#######################################
#               IMDS                  #
#######################################
run_imds <- function(facedata, dim = 2, seed = 12345,id_col="Subject" ) {
  library(mdsOpt)
  
  # Step 1: 取得變數前綴與中心/半徑矩陣
  var_names <- names(facedata)[-1]
  prefixes <- unique(gsub("_(Lower|Upper)", "", var_names))
  
  center <- sapply(prefixes, function(p) {
    (facedata[[paste0(p, "_Lower")]] + facedata[[paste0(p, "_Upper")]]) / 2
  })
  radius <- sapply(prefixes, function(p) {
    (facedata[[paste0(p, "_Upper")]] - facedata[[paste0(p, "_Lower")]]) / 2
  })
  
  # Step 2: 用 idistBox() 直接計算 interval distance
  IDM <- idistBox(X = center, R = radius)
  
  # Step 3: 執行 IMDS
  set.seed(seed)
  result <- IMDS(IDM, p = dim, eps = 1e-5, maxit = 1000, model = "box",
                 opt.method = "MM", ini = "auto", report = 100,
                 grad.num = FALSE, rel = 0, dil = 1)
  
  # Step 4: 整理輸出資料框
  df <- as.data.frame(result$X)
  colnames(df) <- paste0("Dim", 1:dim)
  df$Group <- facedata[[id_col]]
  df$Width <- 2 * result$R[, 1]
  df$Height <- 2 * result$R[, 2]
  df$Dim1_Lower <- df$Dim1 - df$Width / 2
  df$Dim1_Upper <- df$Dim1 + df$Width / 2
  df$Dim2_Lower <- df$Dim2 - df$Height / 2
  df$Dim2_Upper <- df$Dim2 + df$Height / 2
  return(df)
}








































































































