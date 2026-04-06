# ============================================================================
# 區間資料標準化函數
# 基於 De Carvalho et al. (2006) 的三種標準化方法
# ============================================================================

# 方法 1: 使用區間中心點離散度的標準化
# Standardization using the dispersion of the interval centers
standardize_method1 <- function(lower, upper) {
  # 計算中心點
  midpoints <- (lower + upper) / 2
  
  # 計算中心點的平均值 m_j
  m_j <- mean(midpoints)
  
  # 計算中心點的離散度 s_j^2
  s_j_squared <- mean((midpoints - m_j)^2)
  s_j <- sqrt(s_j_squared)
  
  # 轉換區間邊界
  lower_prime <- (lower - m_j) / s_j
  upper_prime <- (upper - m_j) / s_j
  
  return(list(
    lower = lower_prime,
    upper = upper_prime,
    m_j = m_j,
    s_j = s_j
  ))
}

# 方法 2: 使用區間邊界離散度的標準化
# Standardization using the dispersion of the interval boundaries
standardize_method2 <- function(lower, upper) {
  # 計算中心點的平均值 m_j (與方法1相同)
  midpoints <- (lower + upper) / 2
  m_j <- mean(midpoints)
  
  # 計算聯合離散度 tilde(s_j)^2
  s_tilde_j_squared <- mean(((lower - m_j)^2 + (upper - m_j)^2) / 2)
  s_tilde_j <- sqrt(s_tilde_j_squared)
  
  # 轉換區間邊界
  lower_prime <- (lower - m_j) / s_tilde_j
  upper_prime <- (upper - m_j) / s_tilde_j
  
  return(list(
    lower = lower_prime,
    upper = upper_prime,
    m_j = m_j,
    s_tilde_j = s_tilde_j
  ))
}

# 方法 3: 使用全域範圍的標準化 (Min-Max 標準化)
# Standardization using the global range
standardize_method3 <- function(lower, upper) {
  # 找出全域的最小下界與最大上界
  Min_j <- min(lower)
  Max_j <- max(upper)
  
  # 計算範圍
  range_j <- Max_j - Min_j
  
  # 轉換區間邊界到 [0, 1]
  lower_prime <- (lower - Min_j) / range_j
  upper_prime <- (upper - Min_j) / range_j
  
  return(list(
    lower = lower_prime,
    upper = upper_prime,
    Min_j = Min_j,
    Max_j = Max_j
  ))
}

# ============================================================================
# 批次處理函數：對整個資料框進行標準化
# ============================================================================

# ============================================================================
# 批次處理函數：可指定分類欄位名稱 id_col（例如 "species", "Group", "Class"...）
# ============================================================================
standardize_interval_data <- function(data, 
                                      id_col = NULL,   # 這裡改成參數
                                      method = 1) {
  if (!method %in% 1:3) {
    stop("方法必須是 1, 2, 或 3")
  }
  
  std_function <- switch(
    as.character(method),
    `1` = standardize_method1,
    `2` = standardize_method2,
    `3` = standardize_method3
  )
  
  col_names <- names(data)
  n <- nrow(data)
  
  # 只抓真正是區間欄位（結尾是 _Lower 或 _Upper）
  interval_cols <- grep("_(Lower|Upper)$", col_names, value = TRUE)
  variable_prefixes <- unique(sub("_(Lower|Upper)$", "", interval_cols))
  
  # 先建立結果 data.frame：如果有給 id_col 且存在，就複製那一欄
  if (!is.null(id_col) && id_col %in% col_names) {
    result <- data.frame(data[id_col], check.names = FALSE)
  } else {
    # 沒指定或找不到，就只給 row_id（但你說一定有分類欄位，就不會用到這裡）
    result <- data.frame(row_id = seq_len(n))
  }
  
  # 對每個變數前綴做標準化
  for (var_prefix in variable_prefixes) {
    lower_col <- paste0(var_prefix, "_Lower")
    upper_col <- paste0(var_prefix, "_Upper")
    
    if (lower_col %in% col_names && upper_col %in% col_names) {
      lower_vec <- data[[lower_col]]
      upper_vec <- data[[upper_col]]
      
      std_result <- std_function(lower_vec, upper_vec)
      
      result[[lower_col]] <- std_result$lower
      result[[upper_col]] <- std_result$upper
    }
  }
  
  return(result)
}
