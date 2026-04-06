###############################
#            畫圖             #
###############################
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
    #coord_equal() +  # 固定比例，避免方框變形
    coord_fixed(ratio = 1)+
    scale_colour_manual(values = setNames(cols, lvls), drop = FALSE) +
    labs(title = title, x = xlab, y = ylab, colour = "Group") +
    theme_classic(base_size = base_size) +
    theme(
      legend.position = "bottom",
      plot.title = element_text(hjust = 0.5, face = "bold")
    )
}




################################
#           vmtsne             #
################################
vm_tsne_projection <- function(interval_df,
                               id_col = "Subject",
                               dims = 2,
                               perplexity = 50,
                               theta = 0.5,
                               normalize = FALSE,
                               max_iter = 1000,
                               eta = 200) {
  
 
  
  # ---- Step 0: strict column validation (error if not compliant)
  if (!is.data.frame(interval_df)) stop("interval_df must be a data.frame.")
  if (!id_col %in% names(interval_df)) stop("id_col not found: ", id_col)
  
  var_names <- setdiff(names(interval_df), id_col)
  
  # only allow interval columns ending with _Lower/_Upper (exclude other columns)
  iv_cols <- grep("_(Lower|Upper)$", var_names, value = TRUE)
  
  # If you want to allow extra non-interval columns (besides id_col), comment this block.
  non_iv_cols <- setdiff(var_names, iv_cols)
  if (length(non_iv_cols) > 0) {
    stop("Found non-interval columns (not ending with _Lower/_Upper): ",
         paste(non_iv_cols, collapse = ", "))
  }
  
  # derive prefixes and sort for stable order
  var_prefixes <- sort(unique(sub("_(Lower|Upper)$", "", iv_cols)))
  
  # must have both Lower and Upper for each prefix
  need_L <- paste0(var_prefixes, "_Lower")
  need_U <- paste0(var_prefixes, "_Upper")
  missing_L <- setdiff(need_L, names(interval_df))
  missing_U <- setdiff(need_U, names(interval_df))
  if (length(missing_L) > 0 || length(missing_U) > 0) {
    stop("Unpaired interval columns. Missing: ",
         paste(c(missing_L, missing_U), collapse = ", "))
  }
  
  # ---- Step 1: extract lower/upper matrices in stable order
  lower_mat <- as.matrix(interval_df[, need_L, drop = FALSE])
  upper_mat <- as.matrix(interval_df[, need_U, drop = FALSE])
  colnames(lower_mat) <- var_prefixes
  colnames(upper_mat) <- var_prefixes
  rownames(lower_mat) <- interval_df[[id_col]]
  rownames(upper_mat) <- interval_df[[id_col]]
  
  # safety: Lower <= Upper
  if (any(lower_mat > upper_mat, na.rm = TRUE)) {
    stop("Found Lower > Upper in the input data (after your preprocessing).")
  }
  
  # ---- Step 2: build all vertices (2^p)
  generate_vertices <- function(lower, upper) {
    p <- length(lower)
    expand.grid(lapply(1:p, function(j) c(lower[j], upper[j]))) |> as.matrix()
  }
  
  vertices_list <- list()
  vertex_owner <- integer(0)
  for (i in 1:nrow(lower_mat)) {
    verts <- generate_vertices(lower_mat[i, ], upper_mat[i, ])
    vertices_list[[i]] <- verts
    vertex_owner <- c(vertex_owner, rep(i, nrow(verts)))
  }
  vertices_all <- do.call(rbind, vertices_list)
  
  # (you currently do not standardize here)
  vertices_scaled <- vertices_all
  
  # ---- Step 3: t-SNE
  set.seed(12345)
  tsne_model <- Rtsne::Rtsne(vertices_scaled,
                             dims = dims,
                             perplexity = perplexity,
                             theta = theta,
                             pca = FALSE,
                             normalize = normalize,
                             max_iter = max_iter,
                             eta = eta,
                             verbose = TRUE,
                             check_duplicates = FALSE)
  
  proj_vertices <- tsne_model$Y
  
  # ---- Step 4: reconstruct interval per sample via owner
  n <- nrow(lower_mat)
  proj_lower <- matrix(NA_real_, n, dims)
  proj_upper <- matrix(NA_real_, n, dims)
  for (i in 1:n) {
    idx <- which(vertex_owner == i)
    proj_i <- proj_vertices[idx, , drop = FALSE]
    proj_lower[i, ] <- apply(proj_i, 2, min)
    proj_upper[i, ] <- apply(proj_i, 2, max)
  }
  
  # ---- Step 5: output
  result <- data.frame(
    Group      = interval_df[[id_col]],
    Dim1_Lower = proj_lower[, 1],
    Dim1_Upper = proj_upper[, 1],
    Dim2_Lower = proj_lower[, 2],
    Dim2_Upper = proj_upper[, 2]
  )
  #result$Dim1 <- (result$Dim1_Lower + result$Dim1_Upper) / 2
  #result$Dim2 <- (result$Dim2_Lower + result$Dim2_Upper) / 2
  
  return(result)
}


################################
#        qm_tsne (strict)      #
################################
qm_tsne_projection <- function(data,
                               id_col = "Subject",
                               m = 4,
                               dims = 2,
                               theta = 0.5,
                               perplexity = 40,
                          
                               eta = 200,
                               max_iter = 1000) {
  
  if (!requireNamespace("Rtsne", quietly = TRUE)) {
    stop("Package 'Rtsne' is required.")
  }
  
  ## ---- Step 0: strict column validation
  if (!is.data.frame(data)) stop("data must be a data.frame.")
  if (!id_col %in% names(data)) stop("id_col not found: ", id_col)
  
  var_names <- setdiff(names(data), id_col)
  
  # 只允許 *_Lower / *_Upper
  iv_cols <- grep("_(Lower|Upper)$", var_names, value = TRUE)
  
  non_iv_cols <- setdiff(var_names, iv_cols)
  if (length(non_iv_cols) > 0) {
    stop("Found non-interval columns (not ending with _Lower/_Upper): ",
         paste(non_iv_cols, collapse = ", "))
  }
  
  # prefix：排序後固定
  prefixes <- sort(unique(sub("_(Lower|Upper)$", "", iv_cols)))
  
  # Lower / Upper 必須成對
  need_L <- paste0(prefixes, "_Lower")
  need_U <- paste0(prefixes, "_Upper")
  missing_L <- setdiff(need_L, names(data))
  missing_U <- setdiff(need_U, names(data))
  if (length(missing_L) > 0 || length(missing_U) > 0) {
    stop("Unpaired interval columns. Missing: ",
         paste(c(missing_L, missing_U), collapse = ", "))
  }
  
  ## ---- Step 1: extract lower / upper matrices
  lower_mat <- as.matrix(data[, need_L, drop = FALSE])
  upper_mat <- as.matrix(data[, need_U, drop = FALSE])
  colnames(lower_mat) <- prefixes
  colnames(upper_mat) <- prefixes
  
  n <- nrow(lower_mat)
  
  # safety: Lower <= Upper
  if (any(lower_mat > upper_mat, na.rm = TRUE)) {
    stop("Found Lower > Upper in input data.")
  }
  
  ## ---- Step 2: quantile expansion + owner
  expand_quantiles <- function(lower, upper, m) {
    sapply(0:m, function(k) lower + (upper - lower) * k / m)
  }
  
  expanded_list <- vector("list", n)
  owner <- integer(0)
  
  for (i in 1:n) {
    Xi <- t(expand_quantiles(lower_mat[i, ], upper_mat[i, ], m))  # (m+1) x p
    expanded_list[[i]] <- Xi
    owner <- c(owner, rep(i, nrow(Xi)))
  }
  
  expanded_matrix <- do.call(rbind, expanded_list)
  
  # 防呆：NaN / Inf
  if (any(!is.finite(expanded_matrix))) {
    stop("Inf / NaN found in expanded_matrix.")
  }
  
  ## ---- Step 3: perplexity sanity check (expanded points)
  N <- nrow(expanded_matrix)
  perp_max <- floor((N - 1) / 3)
  if (perplexity > perp_max) {
    stop("perplexity too large for expanded points: ",
         perplexity, " > floor((N-1)/3) = ", perp_max)
  }
  
  ## ---- Step 4: t-SNE
  set.seed(12345)
  tsne_result <- Rtsne::Rtsne(expanded_matrix,
                              dims = dims,
                              perplexity = perplexity,
                              theta = theta,
                              pca = FALSE,
                              normalize = FALSE,
                              eta = eta,
                              max_iter = max_iter,
                              verbose = TRUE,
                              check_duplicates = FALSE)
  
  tsne_points <- tsne_result$Y
  
  ## ---- Step 5: owner-based interval reconstruction
  proj_lower <- matrix(NA_real_, n, dims)
  proj_upper <- matrix(NA_real_, n, dims)
  
  for (i in 1:n) {
    idx <- which(owner == i)
    Zi  <- tsne_points[idx, , drop = FALSE]
    proj_lower[i, ] <- apply(Zi, 2, min)
    proj_upper[i, ] <- apply(Zi, 2, max)
  }
  
  ## ---- Step 6: output
  result <- data.frame(
    Group      = data[[id_col]],
    Dim1_Lower = proj_lower[, 1],
    Dim1_Upper = proj_upper[, 1],
    Dim2_Lower = proj_lower[, 2],
    Dim2_Upper = proj_upper[, 2]
  )
  
  #result$Dim1 <- (result$Dim1_Lower + result$Dim1_Upper) / 2
  #result$Dim2 <- (result$Dim2_Lower + result$Dim2_Upper) / 2
  
  return(result)
}



###################################
#          LU懲罰項  (修改後)     #
###################################
lu1_tsne_ab_projection <- function(
    data,
    id_col = "Subject",
    dims = 2,
    perplexity = 30,
    alpha = 0.5,
    learning_rate = 50,
    max_iter = 1000,
    initial_P_gain = 1,
    momentum = 0.5,
    final_momentum = 0.8,
    mom_switch_iter = 250,
    penalty_lambda = 1,
    seed = 12345,
    verbose = TRUE,
    
    init_a = NULL,     # n x d
    init_b = NULL,     # n x d
    init_gap = 1e-3    # 當只給 init_a 時使用
) {
  set.seed(seed)
  stopifnot(dims == 2)
  
  # --- 資料處理 ---
  lower_names <- grep("_Lower$", names(data), value = TRUE)
  upper_names <- grep("_Upper$", names(data), value = TRUE)
  bases <- intersect(sub("_Lower$", "", lower_names),
                     sub("_Upper$", "", upper_names))
  if (length(bases) == 0)
    stop("找不到成對的 *_Lower/*_Upper 欄位。")
  
  lower_mat <- as.matrix(data[, paste0(bases, "_Lower")])
  upper_mat <- as.matrix(data[, paste0(bases, "_Upper")])
  n <- nrow(lower_mat); d <- dims
  
  # =====================================================
  # 初始化（可指定，否則用你原本的）
  # =====================================================
  if (!is.null(init_a)) {
    if (!is.matrix(init_a)) init_a <- as.matrix(init_a)
    if (!all(dim(init_a) == c(n, d)))
      stop("init_a 必須是 n x d 矩陣")
    
    a_prime <- init_a
    
    if (!is.null(init_b)) {
      if (!is.matrix(init_b)) init_b <- as.matrix(init_b)
      if (!all(dim(init_b) == c(n, d)))
        stop("init_b 必須是 n x d 矩陣")
      b_prime <- init_b
    } else {
      b_prime <- a_prime + init_gap
    }
    
  } else {
    a_prime <- matrix(rnorm(n * d, 0, 0.01), n, d)
    b_prime <- a_prime + 1e-3
  }
  
  # --- 高維相似度 P ---
  compute_P <- function(X, perplexity = 30) {
    n <- nrow(X)
    D <- as.matrix(dist(X))^2
    target_entropy <- log(perplexity)
    P <- matrix(0, n, n)
    for (i in 1:n) {
      betamin <- -Inf; betamax <- Inf; beta <- 1
      Di <- D[i, -i]
      for (iter in 1:50) {
        Pi <- exp(-Di * beta)
        sumPi <- sum(Pi); if (sumPi == 0) sumPi <- .Machine$double.eps
        Pi <- Pi / sumPi
        H <- -sum(Pi * log(pmax(Pi, .Machine$double.eps)))
        Hdiff <- H - target_entropy
        if (abs(Hdiff) < 1e-5) break
        if (Hdiff > 0) {
          betamin <- beta
          beta <- ifelse(is.infinite(betamax), beta * 2,
                         (beta + betamax) / 2)
        } else {
          betamax <- beta
          beta <- ifelse(is.infinite(betamin), beta / 2,
                         (beta + betamin) / 2)
        }
      }
      P[i, -i] <- Pi
    }
    (P + t(P)) / (2 * n)
  }
  
  P_L <- compute_P(lower_mat, perplexity) * initial_P_gain
  P_U <- compute_P(upper_mat, perplexity) * initial_P_gain
  
  # --- 動量與紀錄 ---
  inc_a <- matrix(0, n, d)
  inc_b <- matrix(0, n, d)
  loss_history <- numeric(max_iter)
  eps <- 1e-12
  
  for (iter in 1:max_iter) {
    # --- 低維 Q ---
    Q_num_L <- 1 / (1 + as.matrix(dist(a_prime))^2)
    Q_num_U <- 1 / (1 + as.matrix(dist(b_prime))^2)
    diag(Q_num_L) <- 0
    diag(Q_num_U) <- 0
    Q_L <- Q_num_L / sum(Q_num_L)
    Q_U <- Q_num_U / sum(Q_num_U)
    
    # --- KL ---
    kl_L <- sum(P_L * log((P_L + eps) / (Q_L + eps)))
    kl_U <- sum(P_U * log((P_U + eps) / (Q_U + eps)))
    
    # =====================================================
    # LaTeX 的 Log-Hinge penalty:
    #   phi = (1/(nd)) * sum log(1 + max(0, S_ik)),  S=a-b
    # =====================================================
    S <- a_prime - b_prime
    hinge <- pmax(0, S)
    penalty_term <- mean(log1p(hinge))  # log(1+hinge)
    
    loss <- alpha * kl_L + (1 - alpha) * kl_U +
      penalty_lambda * penalty_term
    loss_history[iter] <- loss
    
    # --- KL gradients (你原本寫法保留) ---
    grad_a <- matrix(0, n, d)
    grad_b <- matrix(0, n, d)
    for (i in 1:n) {
      SL <- numeric(d); SU <- numeric(d)
      for (j in 1:n) {
        if (i == j) next
        fL <- (P_L[i, j] - Q_L[i, j]) * Q_num_L[i, j]
        fU <- (P_U[i, j] - Q_U[i, j]) * Q_num_U[i, j]
        SL <- SL + fL * (a_prime[i, ] - a_prime[j, ])
        SU <- SU + fU * (b_prime[i, ] - b_prime[j, ])
      }
      grad_a[i, ] <- 4 * alpha * SL
      grad_b[i, ] <- 4 * (1 - alpha) * SU
    }
    
    # =====================================================
    # LaTeX penalty gradients:
    #   dphi/da_ik = (1/(nd)) * I(S>0)/(1+S)
    #   dphi/db_ik = -(1/(nd)) * I(S>0)/(1+S)
    # =====================================================
    I <- (S > 0)
    denom <- 1 + S
    denom[!I] <- 1  # 避免除到無意義的位置
    grad_pen_a <- (I / denom) / (n * d)
    grad_pen_b <- -(I / denom) / (n * d)
    
    grad_a <- grad_a + penalty_lambda * grad_pen_a
    grad_b <- grad_b + penalty_lambda * grad_pen_b
    
    # --- momentum update ---
    mom <- if (iter < mom_switch_iter) momentum else final_momentum
    inc_a <- mom * inc_a - learning_rate * grad_a
    inc_b <- mom * inc_b - learning_rate * grad_b
    a_prime <- a_prime + inc_a
    b_prime <- b_prime + inc_b
    
    if (iter == 100) {
      P_L <- P_L / initial_P_gain
      P_U <- P_U / initial_P_gain
    }
    
    if (verbose && iter %% 50 == 0) {
      cat(sprintf("Iter %4d | Loss %.6f | Mean width %.6f | Viol %.0f\n",
                  iter, loss, mean(b_prime - a_prime), sum(a_prime >= b_prime)))
    }
  }
  
  viol_mat <- (a_prime >= b_prime)
  if (any(viol_mat)) {
    warning(sprintf(
      "最終嵌入中有 %d 個端點違反 Lower < Upper。",
      sum(viol_mat)
    ))
  }
  
  list(
    embedding = data.frame(
      Group = if (!is.null(data[[id_col]])) data[[id_col]] else seq_len(n),
      Dim1_Lower = a_prime[,1],
      Dim2_Lower = a_prime[,2],
      Dim1_Upper = b_prime[,1],
      Dim2_Upper = b_prime[,2]
    ),
    loss_history = loss_history
  )
}



###################################
#           wasserstein-CR        #
###################################
#       CR                        #
###################################
cr1_tsne_wass_projection <- function(data,
                                     id_col = "Subject",
                                     dims = 2,
                                     perplexity = 30,
                                     lambda = 1.0,   # 半徑權重（lambda=1 時就是理論 W2）
                                     eta = 5,
                                     max_iter = 1000,
                                     EE = 4.0,       # early exaggeration 倍數（=1 表關閉）
                                     T_EE = 250,     # early exaggeration 結束迭代
                                     momentum = 0.5,
                                     final_momentum = 0.8,
                                     mom_switch_iter = 250,
                                     tol = 1e-5,
                                     beta_max_iter = 50,
                                     alpha = 1,      # softplus 的全域尺度 α
                                     seed = 12345,
                                     verbose = TRUE,

                                     init_C = NULL,        # n x d 低維中心初始值
                                     init_S = NULL,        # n x d S 參數初始值
                                     init_r = NULL,        # 可選：用 r' 初始值指定（scalar / length n / nxd）
                                     init_inc_C = NULL,    # n x d 動量累積器初始值
                                     init_inc_S = NULL     # n x d 動量累積器初始值
) {
  stopifnot(dims >= 1)
  set.seed(seed)
  
  ## ---------- 1. 取出高維 L/U（嚴謹對齊前綴） ----------
  lower_names <- grep("_Lower$", names(data), value = TRUE)
  upper_names <- grep("_Upper$", names(data), value = TRUE)
  bases <- intersect(sub("_Lower$", "", lower_names),
                     sub("_Upper$", "", upper_names))
  if (!length(bases)) stop("找不到成對的 *_Lower/*_Upper 欄位。")
  
  L_high <- as.matrix(data[, paste0(bases, "_Lower"), drop = FALSE])
  U_high <- as.matrix(data[, paste0(bases, "_Upper"), drop = FALSE])
  n <- nrow(L_high); d <- as.integer(dims)
  
  ## ---------- 2. 高維：以 Wasserstein (W2) CR 距離建立 P ----------
  # W2^2 = ||m_i - m_j||^2 + (lambda^2/3) ||r_i - r_j||^2
  compute_P_W2_CR <- function(L, U,
                              perplexity = 30,
                              lambda = 1.0,
                              tol = 1e-5, max_iter = 50) {
    n <- nrow(L)
    M <- (L + U) / 2   # 中心向量 m
    R <- (U - L) / 2   # 半徑向量 r
    
    # 中心差的平方距離矩陣
    sumM <- rowSums(M^2)
    Dm2  <- outer(sumM, sumM, "+") - 2 * (M %*% t(M))
    
    # 半徑差的平方距離矩陣
    sumR <- rowSums(R^2)
    Dr2  <- outer(sumR, sumR, "+") - 2 * (R %*% t(R))
    
    coef_r <- (lambda^2) / 3
    D2     <- pmax(Dm2 + coef_r * Dr2, 0)
    diag(D2) <- 0
    
    P <- matrix(0, n, n)
    logPerp <- log(perplexity)
    eps <- .Machine$double.eps
    
    for (i in 1:n) {
      betamin <- -Inf; betamax <- Inf; beta <- 1
      Di <- D2[i, -i]
      for (iter2 in 1:max_iter) {
        t  <- -beta * Di
        t  <- t - max(t)                 # 數值穩定 softmax
        Pi <- exp(t)
        s  <- sum(Pi)
        if (!is.finite(s) || s <= eps) { Pi[] <- 1; s <- length(Pi) }
        Pi <- Pi / s
        
        H  <- -sum(Pi * log(pmax(Pi, eps)))
        Hdiff <- H - logPerp
        if (abs(Hdiff) < tol) break
        if (Hdiff > 0) {
          betamin <- beta
          beta   <- ifelse(is.infinite(betamax), beta * 2, (beta + betamax)/2)
        } else {
          betamax <- beta
          beta   <- ifelse(is.infinite(betamin), beta / 2, (beta + betamin)/2)
        }
      }
      P[i, -i] <- Pi
    }
    (P + t(P)) / (2 * n)
  }
  
  P <- compute_P_W2_CR(L_high, U_high,
                       perplexity = perplexity,
                       lambda = lambda,
                       tol = tol, max_iter = beta_max_iter)
  
  ## ---------- 3. 初始化（支援外部輸入；未輸入則用原本預設） ----------
  # 小工具：inverse softplus（對應你原本 log(exp(x)-1)）
  inv_softplus <- function(y) log(pmax(exp(y) - 1, 1e-12))
  
  # 3.1 初始化 C_prime
  if (!is.null(init_C)) {
    if (!is.matrix(init_C)) init_C <- as.matrix(init_C)
    if (!all(dim(init_C) == c(n, d))) {
      stop(sprintf("init_C 維度需為 n x d = %d x %d", n, d))
    }
    C_prime <- init_C
  } else {
    # 你原本：小隨機
    C_prime <- matrix(rnorm(n * d, 0, 1e-4), n, d)
  }
  
  # 3.2 初始化 S_mat（優先順序：init_S > init_r > 你原本預設）
  if (!is.null(init_S)) {
    if (!is.matrix(init_S)) init_S <- as.matrix(init_S)
    if (!all(dim(init_S) == c(n, d))) {
      stop(sprintf("init_S 維度需為 n x d = %d x %d", n, d))
    }
    S_mat <- init_S
    
  } else if (!is.null(init_r)) {
    # init_r 可接受：scalar / length n / nxd
    if (length(init_r) == 1) {
      r0 <- matrix(rep(init_r, n * d), n, d)
    } else {
      r0 <- init_r
      if (!is.matrix(r0)) r0 <- as.matrix(r0)
      
      if (length(r0) == n && is.null(dim(r0))) {
        # length n 向量：每個點一個半徑，broadcast 到 d 維
        r0 <- matrix(rep(as.vector(r0), each = d), n, d)
      }
      
      if (!all(dim(r0) == c(n, d))) {
        stop(sprintf("init_r 必須是 scalar、length n=%d，或 n x d=%d x %d", n, n, d))
      }
    }
    r0 <- pmax(r0, 0)                 # 半徑不可負
    S_mat <- inv_softplus(r0 / alpha) # r' = alpha*softplus(S) 的反推
    
  } else {
    # 你原本預設：由高維平均半徑決定 S_init，再加噪音
    R_high <- (U_high - L_high) / 2
    mean_r_high <- mean(R_high)
    init_r0 <- pmax(mean_r_high, 0.1)            # 避免太小
    S_init  <- inv_softplus(init_r0 / alpha)     # softplus^{-1}(r/α)
    S_mat   <- matrix(rnorm(n * d, mean = S_init, sd = 0.1), n, d)
  }
  
  # 3.3 動量累積器（可外部指定；否則 0）
  if (!is.null(init_inc_C)) {
    if (!is.matrix(init_inc_C)) init_inc_C <- as.matrix(init_inc_C)
    if (!all(dim(init_inc_C) == c(n, d))) stop("init_inc_C 維度需為 n x d")
    inc_C <- init_inc_C
  } else {
    inc_C <- matrix(0, n, d)
  }
  
  if (!is.null(init_inc_S)) {
    if (!is.matrix(init_inc_S)) init_inc_S <- as.matrix(init_inc_S)
    if (!all(dim(init_inc_S) == c(n, d))) stop("init_inc_S 維度需為 n x d")
    inc_S <- init_inc_S
  } else {
    inc_S <- matrix(0, n, d)
  }
  
  loss_history <- numeric(max_iter)
  
  ## 小工具：softplus 跟 sigmoid
  softplus <- function(x) log1p(exp(x))            # 你原本定義
  sigmoid  <- function(x) 1 / (1 + exp(-x))
  
  coef_r_low <- (lambda^2) / 3
  
  ## ---------- 4. 主訓練迴圈（梯度下降 + 動量） ----------
  for (iter in 1:max_iter) {
    r_prime  <- alpha * softplus(S_mat)
    
    sumCp   <- rowSums(C_prime^2)
    dist_c  <- outer(sumCp, sumCp, "+") - 2 * (C_prime %*% t(C_prime))
    
    sumRp   <- rowSums(r_prime^2)
    dist_r  <- outer(sumRp, sumRp, "+") - 2 * (r_prime %*% t(r_prime))
    
    dist2   <- pmax(dist_c + coef_r_low * dist_r, 0)
    diag(dist2) <- 0
    
    num      <- 1 / (1 + dist2); diag(num) <- 0
    Q        <- num / sum(num)
    
    if (EE != 1 && iter <= T_EE) {
      P_eff <- EE * P
    } else {
      P_eff <- P
    }
    
    loss     <- sum(P_eff * log(pmax(P_eff, 1e-12) / pmax(Q, 1e-12)))
    loss_history[iter] <- loss
    
    W        <- (P_eff - Q) * num
    
    grad_C   <- 4 * (C_prime * rowSums(W) - W %*% C_prime)
    grad_r   <- 4 * coef_r_low * (r_prime * rowSums(W) - W %*% r_prime)
    
    dR_dS    <- alpha * sigmoid(S_mat)
    grad_S   <- grad_r * dR_dS
    
    mom      <- ifelse(iter < mom_switch_iter, momentum, final_momentum)
    inc_C    <- mom * inc_C - eta * grad_C
    inc_S    <- mom * inc_S - eta * grad_S
    
    C_prime  <- C_prime + inc_C
    S_mat    <- S_mat  + inc_S
    
    if (verbose && iter %% 50 == 0) {
      cat(sprintf("Iter: %4d  Loss: %.6f  mean(r')=%.4f\n",
                  iter, loss, mean(r_prime)))
    }
  }
  
  ## ---------- 5. 輸出低維區間 ----------
  r_prime <- alpha * softplus(S_mat)
  Lower   <- C_prime - r_prime
  Upper   <- C_prime + r_prime
  
  out <- data.frame(
    Group = if (id_col %in% names(data)) data[[id_col]] else seq_len(n),
    check.names = FALSE
  )
  for (k in 1:d) {
    out[[paste0("Dim", k, "_Lower")]]  <- Lower[, k]
    out[[paste0("Dim", k, "_Upper")]]  <- Upper[, k]
    out[[paste0("Dim", k)]]            <- C_prime[, k]
    out[[paste0("Dim", k, "_Radius")]] <- r_prime[, k]
  }
  
  list(embedding = out, loss_history = loss_history)
}



































































































