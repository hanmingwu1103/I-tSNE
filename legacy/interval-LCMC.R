library(ggplot2)
library(patchwork)

########################################
#    wasserstein、hausdorff、minkowski #
########################################

dist_interval <- function(lower_mat, upper_mat,
                          metric = c("wasserstein","hausdorff","minkowski"),
                          gamma = 0.5, p = 2, lambda = 1/3) {
  metric <- match.arg(metric)
  
  if (metric == "wasserstein") {
    center <- (lower_mat + upper_mat) / 2
    radius <- (upper_mat - lower_mat) / 2
    n <- nrow(center)
    D <- matrix(0, n, n)
    for (i in 1:(n-1)) for (j in (i+1):n) {
      d <- sqrt(sum((center[i,]-center[j,])^2) + (1/3)*sum((radius[i,]-radius[j,])^2))
      D[i,j] <- D[j,i] <- d
    }
    return(D)
  }
  
  ## === Ichino–Yaguchi 距離（以 "minkowski" 名稱呼叫） ===
  if (metric == "minkowski") {
    n <- nrow(lower_mat); p_dim <- ncol(lower_mat)
    D <- matrix(0, n, n)
    for (i in 1:(n-1)) {
      for (j in (i+1):n) {
        d2 <- 0
        for (k in 1:p_dim) {
          a1 <- lower_mat[i, k]; b1 <- upper_mat[i, k]
          a2 <- lower_mat[j, k]; b2 <- upper_mat[j, k]
          len1 <- b1 - a1; len2 <- b2 - a2
          union_len <- max(b1, b2) - min(a1, a2)
          inter_len <- max(0, min(b1, b2) - max(a1, a2))
          phi <- (union_len - inter_len) + gamma * (2 * inter_len - len1 - len2)
          d2 <- d2 + phi^2
        }
        D[i, j] <- D[j, i] <- sqrt(d2)
      }
    }
    return(D)
  }
  
  if (metric == "hausdorff") {
    n <- nrow(lower_mat); p_dim <- ncol(lower_mat)
    D <- matrix(0, n, n)
    for (i in 1:(n-1)) for (j in (i+1):n) {
      d <- sqrt(sum(sapply(1:p_dim, function(k)
        max(abs(upper_mat[i,k]-upper_mat[j,k]),
            abs(lower_mat[i,k]-lower_mat[j,k]))^2 )))
      D[i,j] <- D[j,i] <- d
    }
    return(D)
  }
}

#################################################
# 2) KNN/LCMC：預算滿足 k 的最大值，避免重覆排序#
#################################################
neighbor_orders <- function(D) {
  n <- nrow(D)
  lapply(seq_len(n), function(i) {
    di <- D[i, ]; di[i] <- Inf; order(di)
  })
}

lcmc_from_orders <- function(ord_H, ord_L, k) {
  n <- length(ord_H); k <- min(k, n-1)
  hits <- 0L
  for (i in 1:n) {
    hits <- hits + length(intersect(ord_H[[i]][1:k], ord_L[[i]][1:k]))
  }
  qk <- hits / (n * k)
  qk - k/(n - 1)
}

################################################
# Extract low-dim lower/upper from df
################################################
extract_low_LU <- function(low_df) {
  lower_cols <- grep("_Lower$", names(low_df), value = TRUE)
  upper_cols <- grep("_Upper$", names(low_df), value = TRUE)
  
  if (length(lower_cols) == 0 || length(upper_cols) == 0) {
    stop("找不到任何 *_Lower / *_Upper 欄位，請確認低維資料格式。")
  }
  
  base_lower <- sub("_Lower$", "", lower_cols)
  base_upper <- sub("_Upper$", "", upper_cols)
  bases <- intersect(base_lower, base_upper)
  
  if (length(bases) < 2) {
    stop("可配對的低維區間變數不足（少於 2 個）。目前找到的 base：",
         paste(bases, collapse = ", "))
  }
  
  bases <- sort(bases)[1:2]
  
  lower <- as.matrix(low_df[, paste0(bases, "_Lower")])
  upper <- as.matrix(low_df[, paste0(bases, "_Upper")])
  
  message("使用低維區間欄位：",
          paste(paste0(bases, "_Lower"), collapse = ", "),
          " / ",
          paste(paste0(bases, "_Upper"), collapse = ", "))
  
  list(lower = lower, upper = upper)
}

#####################################################
#               主函式                              #
#####################################################
compute_lcmc_tables <- function(high_df,
                                methods_list,
                                k_range = 1:26,
                                gamma = 0.5, p = 2, lambda = 1/3) {
  lower_cols <- grep("_Lower$", names(high_df), value = TRUE)
  upper_cols <- grep("_Upper$", names(high_df), value = TRUE)
  bases <- intersect(sub("_Lower$", "", lower_cols), sub("_Upper$", "", upper_cols))
  lower_H <- as.matrix(high_df[, paste0(bases, "_Lower")])
  upper_H <- as.matrix(high_df[, paste0(bases, "_Upper")])
  
  n <- nrow(lower_H)
  k_range <- k_range[k_range >= 1 & k_range <= (n-1)]
  if (length(k_range) == 0) stop("k_range 與樣本數不符：需有至少一個 k 落在 [1, n-1]。")
  
  metrics <- c("wasserstein","hausdorff","minkowski")
  out_dfs <- vector("list", length(metrics))
  names(out_dfs) <- c("LCMC_wasserstein", "LCMC_hausdorff", "LCMC_minkowski")
  
  for (m_idx in seq_along(metrics)) {
    metric <- metrics[m_idx]
    
    D_H <- dist_interval(lower_H, upper_H, metric = metric, gamma = gamma, p = p, lambda = lambda)
    ord_H <- neighbor_orders(D_H)
    
    all_rows <- list()
    
    for (method_name in names(methods_list)) {
      low_df <- methods_list[[method_name]]
      lu <- extract_low_LU(low_df)
      
      D_L <- dist_interval(lu$lower, lu$upper, metric = metric, gamma = gamma, p = p, lambda = lambda)
      ord_L <- neighbor_orders(D_L)
      
      lcmc_vals <- vapply(k_range, function(k) lcmc_from_orders(ord_H, ord_L, k), numeric(1))
      all_rows[[method_name]] <- data.frame(
        k = k_range,
        Method = method_name,
        LCMC = lcmc_vals,
        stringsAsFactors = FALSE
      )
    }
    
    out_dfs[[m_idx]] <- do.call(rbind, all_rows)
  }
  
  out_dfs
}

#######################################
# 函數 1：plot_lcmc_single（新增 method_name_map）
#######################################
plot_lcmc_single <- function(lcmc_df, metric_title,
                             method_levels = NULL,
                             method_colors = NULL,
                             method_name_map = NULL) {
  
  # ---- NEW: rename methods for legend ----
  if (!is.null(method_name_map)) {
    lcmc_df$Method <- as.character(lcmc_df$Method)
    idx <- lcmc_df$Method %in% names(method_name_map)
    lcmc_df$Method[idx] <- method_name_map[lcmc_df$Method[idx]]
  }
  
  if (!is.null(method_levels)) {
    lcmc_df$Method <- factor(lcmc_df$Method, levels = method_levels)
  } else {
    lcmc_df$Method <- factor(lcmc_df$Method)
  }
  
  p <- ggplot(lcmc_df, aes(x = k, y = LCMC, color = Method, group = Method)) +
    geom_line(linewidth = 1.0) +
    geom_point(size = 1.8) +
    labs(title = metric_title,
         x = "k (number of neighbors)",
         y = "LCMC") +
    theme_minimal(base_size = 13) +
    theme(
      panel.grid.minor = element_blank(),
      legend.position  = "bottom",
      legend.title     = element_blank(),
      legend.text      = element_text(size = 9),
      legend.spacing.x = unit(4, "pt"),
      plot.title       = element_text(hjust = 0.5, face = "bold", size = 14),
      axis.title       = element_text(face = "bold")
    ) +
    guides(color = guide_legend(nrow = 2, byrow = TRUE))
  
  if (!is.null(method_colors)) {
    p <- p + scale_color_manual(values = method_colors, drop = FALSE)
  }
  
  p
}

#######################################
# 函數 2：plot_lcmc_1x3（新增 method_name_map）
#######################################
plot_lcmc_1x3 <- function(lcmc_w_df, lcmc_h_df, lcmc_m_df,
                          out_dir  = ".",
                          out_file = "LCMC_1x3.pdf",
                          save = TRUE,
                          make_dir = TRUE,
                          width = 15, height = 5, dpi = 300,
                          use_cairo = TRUE,
                          method_name_map = NULL) {
  suppressPackageStartupMessages({
    library(ggplot2)
  })
  
  ## ========= 0) output path =========
  stopifnot(is.character(out_dir), length(out_dir) == 1)
  stopifnot(is.character(out_file), length(out_file) == 1)
  
  if (!dir.exists(out_dir)) {
    if (isTRUE(make_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    else stop("out_dir does not exist: ", out_dir)
  }
  
  if (!grepl("\\.(pdf|png)$", tolower(out_file))) {
    out_file <- paste0(out_file, ".pdf")
  }
  out_path <- file.path(out_dir, out_file)
  
  ## ========= 1) add metric labels =========
  lcmc_w_df$Metric <- "Wasserstein"
  lcmc_h_df$Metric <- "Hausdorff"
  lcmc_m_df$Metric <- "Minkowski / Ichino–Yaguchi"
  
  ## ========= 2) merge =========
  all_df <- rbind(lcmc_w_df, lcmc_h_df, lcmc_m_df)
  
  ## ========= 2.5) NEW: rename methods for legend =========
  if (!is.null(method_name_map)) {
    all_df$Method <- as.character(all_df$Method)
    idx <- all_df$Method %in% names(method_name_map)
    all_df$Method[idx] <- method_name_map[all_df$Method[idx]]
  }
  
  ## ========= 3) consistent method levels & colors =========
  all_methods <- sort(unique(as.character(all_df$Method)))
  all_df$Method <- factor(all_df$Method, levels = all_methods)
  
  method_colors <- setNames(
    grDevices::hcl.colors(length(all_methods), palette = "Dark 3"),
    all_methods
  )
  
  ## ========= 4) plot =========
  p <- ggplot(all_df, aes(x = k, y = LCMC, color = Method, group = Method)) +
    geom_line(linewidth = 1) +
    geom_point(size = 1.8) +
    facet_wrap(~ Metric, nrow = 1) +
    scale_color_manual(values = method_colors, drop = FALSE) +
    labs(x = "k (number of neighbors)", y = "LCMC") +
    theme_minimal(base_size = 13) +
    theme(
      panel.grid.minor = element_blank(),
      legend.position  = "bottom",
      legend.title     = element_blank(),
      legend.text      = element_text(size = 9),
      legend.spacing.x = unit(4, "pt"),
      axis.title       = element_text(face = "bold"),
      strip.text       = element_text(face = "bold", size = 12)
    ) +
    guides(color = guide_legend(nrow = 2, byrow = TRUE))
  
  print(p)
  
  ## ========= 5) save =========
  if (isTRUE(save)) {
    ext <- tolower(tools::file_ext(out_path))
    
    dev_fun <- NULL
    if (ext == "pdf") {
      dev_fun <- if (isTRUE(use_cairo) && exists("cairo_pdf", mode = "function")) cairo_pdf else "pdf"
    } else if (ext == "png") {
      dev_fun <- "png"
    } else {
      stop("Unsupported file extension: ", ext, " (use .pdf or .png)")
    }
    
    ggsave(
      filename = out_path,
      plot = p,
      device = dev_fun,
      width = width,
      height = height,
      units = "in",
      dpi = dpi,
      limitsize = FALSE
    )
    message("Saved: ", out_path)
  }
  
  invisible(p)
}
