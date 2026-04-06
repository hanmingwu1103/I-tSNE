############################################
#畫圖實際資料                             #
############################################
plot_interval_pca_tsne_grid <- function(
    vm_pca, qm_pca, cr_pca, imds,
    vm_tsne, qm_tsne, lU2_tsne_embedding, cr1_tsne_embedding,
    out_file = "interval_pca_tsne_comparison.pdf",
    out_dir  = ".",
    make_dir = TRUE,
    panel_in = 2.4,
    pad_pca = 0.02,
    pad_tsne = 0.05,
    legend_height_in = 0.7,
    use_cairo = TRUE,
    print_plot = TRUE,
    legend_ncol = NULL   # NULL = auto; or set to length(all_groups) for 1-row
) {
  stopifnot(is.character(out_file), length(out_file) == 1)
  stopifnot(is.character(out_dir),  length(out_dir)  == 1)
  
  suppressPackageStartupMessages({
    library(ggplot2)
    library(grid)
    library(gtable)
    library(gridExtra)
  })
  
  ## ---- output path
  if (!dir.exists(out_dir)) {
    if (isTRUE(make_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    else stop("out_dir does not exist: ", out_dir)
  }
  out_path <- file.path(out_dir, basename(out_file))
  
  ## ---- themes
  panel_theme <- theme_classic(base_size = 11) +
    theme(
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      plot.margin  = margin(2, 2, 2, 2),
      panel.background = element_blank(),
      panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.6),
      plot.title = element_text(hjust = 0.5, face = "bold", size = 12),
      legend.title = element_text(size = 12, face = "bold"),
      legend.text  = element_text(size = 10)
    )
  
  ## ---- helpers
  get_bounds <- function(df){
    if (all(c("Dim1_Lower","Dim1_Upper","Dim2_Lower","Dim2_Upper") %in% names(df))) {
      c(
        xmin = min(df$Dim1_Lower, na.rm = TRUE),
        xmax = max(df$Dim1_Upper, na.rm = TRUE),
        ymin = min(df$Dim2_Lower, na.rm = TRUE),
        ymax = max(df$Dim2_Upper, na.rm = TRUE)
      )
    } else if (all(c("Lower_X","Upper_X","Lower_Y","Upper_Y") %in% names(df))) {
      c(
        xmin = min(df$Lower_X, na.rm = TRUE),
        xmax = max(df$Upper_X, na.rm = TRUE),
        ymin = min(df$Lower_Y, na.rm = TRUE),
        ymax = max(df$Upper_Y, na.rm = TRUE)
      )
    } else stop("Missing interval columns: need Dim*_Lower/Upper or Lower_*/Upper_*")
  }
  
  square_lims <- function(bounds, pad = 0.02) {
    xr <- bounds["xmax"] - bounds["xmin"]
    yr <- bounds["ymax"] - bounds["ymin"]
    r  <- max(xr, yr) / 2
    cx <- (bounds["xmin"] + bounds["xmax"]) / 2
    cy <- (bounds["ymin"] + bounds["ymax"]) / 2
    r  <- r * (1 + pad)
    list(x = c(cx - r, cx + r), y = c(cy - r, cy + r))
  }
  
  ## ---- global colors
  all_groups <- levels(factor(c(
    vm_pca$Group, qm_pca$Group, cr_pca$Group, imds$Group,
    vm_tsne$Group, qm_tsne$Group, lU2_tsne_embedding$Group, cr1_tsne_embedding$Group
  )))
  cols <- setNames(grDevices::hcl.colors(length(all_groups), "Dark 3"), all_groups)
  
  ## ---- plot builder (NO legend inside panels)
  panel_plot_interval <- function(df, title, lims = NULL) {
    df$Group <- factor(df$Group, levels = all_groups)
    
    if (all(c("Dim1_Lower","Dim1_Upper","Dim2_Lower","Dim2_Upper") %in% names(df))) {
      p <- ggplot(df, aes(colour = Group)) +
        geom_rect(aes(
          xmin = Dim1_Lower, xmax = Dim1_Upper,
          ymin = Dim2_Lower, ymax = Dim2_Upper
        ), fill = NA, linewidth = 0.7)
    } else if (all(c("Lower_X","Upper_X","Lower_Y","Upper_Y") %in% names(df))) {
      p <- ggplot(df, aes(colour = Group)) +
        geom_rect(aes(
          xmin = Lower_X, xmax = Upper_X,
          ymin = Lower_Y, ymax = Upper_Y
        ), fill = NA, linewidth = 0.7)
    } else stop("Missing interval columns for plotting.")
    
    if (!is.null(lims)) {
      p <- p +
        coord_fixed(1) +
        scale_x_continuous(limits = lims$x, expand = expansion(mult = 0)) +
        scale_y_continuous(limits = lims$y, expand = expansion(mult = 0))
    } else {
      p <- p + coord_fixed(1)
    }
    
    p +
      scale_colour_manual(values = cols, drop = FALSE, name = "Group") +
      labs(title = title) +
      panel_theme +
      guides(colour = "none") +            # ★硬關掉 guide
      theme(legend.position = "none")      # ★硬關掉 legend
  }
  
  ## ---- build PCA shared square limits
  b_vm <- get_bounds(vm_pca)
  b_qm <- get_bounds(qm_pca)
  b_cr <- get_bounds(cr_pca)
  b_im <- get_bounds(imds)
  
  bounds_pca_all <- c(
    xmin = min(b_vm["xmin"], b_qm["xmin"], b_cr["xmin"], b_im["xmin"]),
    xmax = max(b_vm["xmax"], b_qm["xmax"], b_cr["xmax"], b_im["xmax"]),
    ymin = min(b_vm["ymin"], b_qm["ymin"], b_cr["ymin"], b_im["ymin"]),
    ymax = max(b_vm["ymax"], b_qm["ymax"], b_cr["ymax"], b_im["ymax"])
  )
  lims_pca_square <- square_lims(bounds_pca_all, pad = pad_pca)
  
  ## ---- build 8 panels
  p1_pca <- panel_plot_interval(vm_pca, "IPCA(VM)", lims_pca_square)
  p2_pca <- panel_plot_interval(qm_pca, "IPCA(QM)", lims_pca_square)
  p3_pca <- panel_plot_interval(cr_pca, "IPCA(CR)", lims_pca_square)
  p4_pca <- panel_plot_interval(imds,   "IMDS",     lims_pca_square)
  
  p1_tsne <- panel_plot_interval(vm_tsne,            "I-tSNE(VM)",  square_lims(get_bounds(vm_tsne),            pad_tsne))
  p2_tsne <- panel_plot_interval(qm_tsne,            "I-tSNE(QM)",  square_lims(get_bounds(qm_tsne),            pad_tsne))
  p3_tsne <- panel_plot_interval(cr1_tsne_embedding, "I-tSNE(CR)", square_lims(get_bounds(cr1_tsne_embedding), pad_tsne))
  p4_tsne <- panel_plot_interval(lU2_tsne_embedding, "I-tSNE(MM)", square_lims(get_bounds(lU2_tsne_embedding), pad_tsne))
  

  ## ---- legend grob (only once)
  carrier <- ggplot(qm_pca, aes(colour = factor(Group, levels = all_groups))) +
    geom_rect(aes(
      xmin = if ("Dim1_Lower" %in% names(qm_pca)) Dim1_Lower else Lower_X,
      xmax = if ("Dim1_Upper" %in% names(qm_pca)) Dim1_Upper else Upper_X,
      ymin = if ("Dim2_Lower" %in% names(qm_pca)) Dim2_Lower else Lower_Y,
      ymax = if ("Dim2_Upper" %in% names(qm_pca)) Dim2_Upper else Upper_Y
    ), fill = NA, linewidth = 0.7) +
    scale_colour_manual(values = cols, drop = FALSE, name = "Group") +
    theme_void() +
    theme(
      legend.position = "bottom",
      legend.title = element_text(size = 12, face = "bold"),
      legend.text  = element_text(size = 10)
    ) +
    guides(colour = guide_legend(
      nrow = 1,
      ncol = legend_ncol,
      byrow = TRUE,
      override.aes = list(linewidth = 1)
    ))
  
  carrier_grob <- ggplotGrob(carrier)
  leg <- gtable_filter(carrier_grob, "guide-box")
  
  ## ---- layout using gridExtra (stable)
  row1 <- arrangeGrob(p1_pca, p2_pca, p3_pca, p4_pca, ncol = 4)
  row2 <- arrangeGrob(p1_tsne, p2_tsne, p3_tsne, p4_tsne, ncol = 4)
  body <- arrangeGrob(row1, row2, ncol = 1, heights = unit.c(unit(1, "null"), unit(1, "null")))
  
  final_grob <- arrangeGrob(
    body,
    leg,
    ncol = 1,
    heights = unit.c(unit(1, "null"), unit(legend_height_in, "in"))
  )
  
  if (isTRUE(print_plot)) {
    grid.newpage()
    grid.draw(final_grob)
  }
  
  ## ---- export (tight PDF)
  pdf_w <- 4 * panel_in
  pdf_h <- 2 * panel_in + legend_height_in
  
  if (isTRUE(use_cairo) && exists("cairo_pdf", mode = "function")) {
    cairo_pdf(filename = out_path, width = pdf_w, height = pdf_h)
    grid.draw(final_grob)
    dev.off()
  } else {
    pdf(file = out_path, width = pdf_w, height = pdf_h)
    grid.draw(final_grob)
    dev.off()
  }
  
  message("Saved: ", out_path)
  invisible(final_grob)
}




#############################################
# 畫圖模擬                                  #
#############################################

plot_interval_tsne_grid_2x4 <- function(
    vm_tsne_2d, qm_tsne_2d, lU2_tsne_2d, cr1_tsne_2d,
    vm_tsne_noise, qm_tsne_noise, lU2_tsne_noise, cr1_tsne_noise,
    panel_in = 2.4,
    pad_tsne = 0.05,
    legend_height_in = 0.7,
    print_plot = TRUE,
    legend_ncol = NULL
) {
  suppressPackageStartupMessages({
    library(ggplot2)
    library(grid)
    library(gtable)
    library(gridExtra)
  })
  
  ## ---- themes
  panel_theme <- theme_classic(base_size = 11) +
    theme(
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      plot.margin  = margin(2, 2, 2, 2),
      panel.background = element_blank(),
      panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.6),
      plot.title = element_text(hjust = 0.5, face = "bold", size = 12),
      legend.title = element_text(size = 12, face = "bold"),
      legend.text  = element_text(size = 10)
    )
  
  ## ---- helpers
  get_bounds <- function(df){
    if (all(c("Dim1_Lower","Dim1_Upper","Dim2_Lower","Dim2_Upper") %in% names(df))) {
      c(
        xmin = min(df$Dim1_Lower, na.rm = TRUE),
        xmax = max(df$Dim1_Upper, na.rm = TRUE),
        ymin = min(df$Dim2_Lower, na.rm = TRUE),
        ymax = max(df$Dim2_Upper, na.rm = TRUE)
      )
    } else if (all(c("Lower_X","Upper_X","Lower_Y","Upper_Y") %in% names(df))) {
      c(
        xmin = min(df$Lower_X, na.rm = TRUE),
        xmax = max(df$Upper_X, na.rm = TRUE),
        ymin = min(df$Lower_Y, na.rm = TRUE),
        ymax = max(df$Upper_Y, na.rm = TRUE)
      )
    } else stop("Missing interval columns: need Dim*_Lower/Upper or Lower_*/Upper_*")
  }
  
  square_lims <- function(bounds, pad = 0.05) {
    xr <- bounds["xmax"] - bounds["xmin"]
    yr <- bounds["ymax"] - bounds["ymin"]
    r  <- max(xr, yr) / 2
    cx <- (bounds["xmin"] + bounds["xmax"]) / 2
    cy <- (bounds["ymin"] + bounds["ymax"]) / 2
    r  <- r * (1 + pad)
    list(x = c(cx - r, cx + r), y = c(cy - r, cy + r))
  }
  
  ## ---- global colors
  all_groups <- levels(factor(c(
    vm_tsne_2d$Group,
    qm_tsne_2d$Group,
    lU2_tsne_2d$Group,
    cr1_tsne_2d$Group,
    vm_tsne_noise$Group,
    qm_tsne_noise$Group,
    lU2_tsne_noise$Group,
    cr1_tsne_noise$Group
  )))
  cols <- setNames(grDevices::hcl.colors(length(all_groups), "Dark 3"), all_groups)
  
  ## ---- plot builder
  panel_plot_interval <- function(df, title, lims = NULL) {
    df$Group <- factor(df$Group, levels = all_groups)
    
    if (all(c("Dim1_Lower","Dim1_Upper","Dim2_Lower","Dim2_Upper") %in% names(df))) {
      p <- ggplot(df, aes(colour = Group)) +
        geom_rect(aes(
          xmin = Dim1_Lower, xmax = Dim1_Upper,
          ymin = Dim2_Lower, ymax = Dim2_Upper
        ), fill = NA, linewidth = 0.7)
    } else if (all(c("Lower_X","Upper_X","Lower_Y","Upper_Y") %in% names(df))) {
      p <- ggplot(df, aes(colour = Group)) +
        geom_rect(aes(
          xmin = Lower_X, xmax = Upper_X,
          ymin = Lower_Y, ymax = Upper_Y
        ), fill = NA, linewidth = 0.7)
    } else stop("Missing interval columns for plotting.")
    
    if (!is.null(lims)) {
      p <- p +
        coord_fixed(1) +
        scale_x_continuous(limits = lims$x, expand = expansion(mult = 0)) +
        scale_y_continuous(limits = lims$y, expand = expansion(mult = 0))
    } else {
      p <- p + coord_fixed(1)
    }
    
    p +
      scale_colour_manual(values = cols, drop = FALSE, name = "Group") +
      labs(title = title) +
      panel_theme +
      guides(colour = "none") +
      theme(legend.position = "none")
  }
  
  ## ---- first row: 2D
  p1 <- panel_plot_interval(
    vm_tsne_2d, "I-tSNE(VM)",
    square_lims(get_bounds(vm_tsne_2d), pad = pad_tsne)
  )
  p2 <- panel_plot_interval(
    qm_tsne_2d, "I-tSNE(QM)",
    square_lims(get_bounds(qm_tsne_2d), pad = pad_tsne)
  )
  p3 <- panel_plot_interval(
    lU2_tsne_2d, "I-tSNE(MM)",
    square_lims(get_bounds(lU2_tsne_2d), pad = pad_tsne)
  )
  p4 <- panel_plot_interval(
    cr1_tsne_2d, "I-tSNE(CR)",
    square_lims(get_bounds(cr1_tsne_2d), pad = pad_tsne)
  )
  
  ## ---- second row: Noise
  p5 <- panel_plot_interval(
    vm_tsne_noise, "I-tSNE(VM)",
    square_lims(get_bounds(vm_tsne_noise), pad = pad_tsne)
  )
  p6 <- panel_plot_interval(
    qm_tsne_noise, "I-tSNE(QM)",
    square_lims(get_bounds(qm_tsne_noise), pad = pad_tsne)
  )
  p7 <- panel_plot_interval(
    lU2_tsne_noise, "I-tSNE(MM)",
    square_lims(get_bounds(lU2_tsne_noise), pad = pad_tsne)
  )
  p8 <- panel_plot_interval(
    cr1_tsne_noise, "I-tSNE(CR)",
    square_lims(get_bounds(cr1_tsne_noise), pad = pad_tsne)
  )
  
  ## ---- 直接排成 2x4，不加左側標題
  row1 <- arrangeGrob(p1, p2, p3, p4, ncol = 4)
  row2 <- arrangeGrob(p5, p6, p7, p8, ncol = 4)
  
  ## ---- legend grob
  carrier_df <- cr1_tsne_noise
  carrier_df$Group <- factor(carrier_df$Group, levels = all_groups)
  
  carrier <- ggplot(carrier_df, aes(colour = Group)) +
    {
      if (all(c("Dim1_Lower","Dim1_Upper","Dim2_Lower","Dim2_Upper") %in% names(carrier_df))) {
        geom_rect(aes(
          xmin = Dim1_Lower, xmax = Dim1_Upper,
          ymin = Dim2_Lower, ymax = Dim2_Upper
        ), fill = NA, linewidth = 0.7)
      } else {
        geom_rect(aes(
          xmin = Lower_X, xmax = Upper_X,
          ymin = Lower_Y, ymax = Upper_Y
        ), fill = NA, linewidth = 0.7)
      }
    } +
    scale_colour_manual(values = cols, drop = FALSE, name = "Group") +
    theme_void() +
    theme(
      legend.position = "bottom",
      legend.title = element_text(size = 12, face = "bold"),
      legend.text  = element_text(size = 10)
    ) +
    guides(colour = guide_legend(
      nrow = 1,
      ncol = legend_ncol,
      byrow = TRUE,
      override.aes = list(linewidth = 1)
    ))
  
  carrier_grob <- ggplotGrob(carrier)
  leg <- gtable_filter(carrier_grob, "guide-box")
  
  ## ---- final layout
  final_grob <- arrangeGrob(
    row1,
    row2,
    leg,
    ncol = 1,
    heights = unit.c(
      unit(1, "null"),
      unit(1, "null"),
      unit(legend_height_in, "in")
    )
  )
  
  if (isTRUE(print_plot)) {
    grid.newpage()
    grid.draw(final_grob)
  }
  
  invisible(final_grob)
}


#############################################
# vm不能用的圖                              #
#############################################
plot_interval_grid_6 <- function(
    qm_pca, cr_pca, imds,
    qm_tsne, lU2_tsne_embedding, cr1_tsne_embedding,
    out_file = "interval_pca_tsne_comparison_6panels.pdf",
    out_dir  = ".",
    make_dir = TRUE,
    panel_in = 2.8,
    pad_pca = 0.02,
    pad_tsne = 0.05,
    legend_height_in = 0.7,
    use_cairo = TRUE,
    print_plot = TRUE,
    legend_ncol = NULL   # NULL = auto
) {
  stopifnot(is.character(out_file), length(out_file) == 1)
  stopifnot(is.character(out_dir),  length(out_dir)  == 1)
  
  suppressPackageStartupMessages({
    library(ggplot2)
    library(grid)
    library(gtable)
    library(gridExtra)
  })
  
  ## ---- output path
  if (!dir.exists(out_dir)) {
    if (isTRUE(make_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    else stop("out_dir does not exist: ", out_dir)
  }
  out_path <- file.path(out_dir, basename(out_file))
  
  ## ---- themes
  panel_theme <- theme_classic(base_size = 11) +
    theme(
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      plot.margin  = margin(2, 2, 2, 2),
      panel.background = element_blank(),
      panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.6),
      plot.title = element_text(hjust = 0.5, face = "bold", size = 12),
      legend.title = element_text(size = 12, face = "bold"),
      legend.text  = element_text(size = 10)
    )
  
  ## ---- helpers
  get_bounds <- function(df){
    if (all(c("Dim1_Lower","Dim1_Upper","Dim2_Lower","Dim2_Upper") %in% names(df))) {
      c(
        xmin = min(df$Dim1_Lower, na.rm = TRUE),
        xmax = max(df$Dim1_Upper, na.rm = TRUE),
        ymin = min(df$Dim2_Lower, na.rm = TRUE),
        ymax = max(df$Dim2_Upper, na.rm = TRUE)
      )
    } else if (all(c("Lower_X","Upper_X","Lower_Y","Upper_Y") %in% names(df))) {
      c(
        xmin = min(df$Lower_X, na.rm = TRUE),
        xmax = max(df$Upper_X, na.rm = TRUE),
        ymin = min(df$Lower_Y, na.rm = TRUE),
        ymax = max(df$Upper_Y, na.rm = TRUE)
      )
    } else {
      stop("Missing interval columns: need Dim*_Lower/Upper or Lower_*/Upper_*")
    }
  }
  
  square_lims <- function(bounds, pad = 0.02) {
    xr <- bounds["xmax"] - bounds["xmin"]
    yr <- bounds["ymax"] - bounds["ymin"]
    r  <- max(xr, yr) / 2
    cx <- (bounds["xmin"] + bounds["xmax"]) / 2
    cy <- (bounds["ymin"] + bounds["ymax"]) / 2
    r  <- r * (1 + pad)
    list(x = c(cx - r, cx + r), y = c(cy - r, cy + r))
  }
  
  ## ---- global colors
  all_groups <- levels(factor(c(
    qm_pca$Group, cr_pca$Group, imds$Group,
    qm_tsne$Group, lU2_tsne_embedding$Group, cr1_tsne_embedding$Group
  )))
  cols <- setNames(grDevices::hcl.colors(length(all_groups), "Dark 3"), all_groups)
  
  ## ---- plot builder
  panel_plot_interval <- function(df, title, lims = NULL) {
    df$Group <- factor(df$Group, levels = all_groups)
    
    if (all(c("Dim1_Lower","Dim1_Upper","Dim2_Lower","Dim2_Upper") %in% names(df))) {
      p <- ggplot(df, aes(colour = Group)) +
        geom_rect(aes(
          xmin = Dim1_Lower, xmax = Dim1_Upper,
          ymin = Dim2_Lower, ymax = Dim2_Upper
        ), fill = NA, linewidth = 0.7)
    } else if (all(c("Lower_X","Upper_X","Lower_Y","Upper_Y") %in% names(df))) {
      p <- ggplot(df, aes(colour = Group)) +
        geom_rect(aes(
          xmin = Lower_X, xmax = Upper_X,
          ymin = Lower_Y, ymax = Upper_Y
        ), fill = NA, linewidth = 0.7)
    } else {
      stop("Missing interval columns for plotting.")
    }
    
    if (!is.null(lims)) {
      p <- p +
        coord_fixed(1) +
        scale_x_continuous(limits = lims$x, expand = expansion(mult = 0)) +
        scale_y_continuous(limits = lims$y, expand = expansion(mult = 0))
    } else {
      p <- p + coord_fixed(1)
    }
    
    p +
      scale_colour_manual(values = cols, drop = FALSE, name = "Group") +
      labs(title = title) +
      panel_theme +
      guides(colour = "none") +
      theme(legend.position = "none")
  }
  
  ## ---- PCA shared limits (3 張共用)
  b_qm <- get_bounds(qm_pca)
  b_cr <- get_bounds(cr_pca)
  b_im <- get_bounds(imds)
  
  bounds_pca_all <- c(
    xmin = min(b_qm["xmin"], b_cr["xmin"], b_im["xmin"]),
    xmax = max(b_qm["xmax"], b_cr["xmax"], b_im["xmax"]),
    ymin = min(b_qm["ymin"], b_cr["ymin"], b_im["ymin"]),
    ymax = max(b_qm["ymax"], b_cr["ymax"], b_im["ymax"])
  )
  lims_pca_square <- square_lims(bounds_pca_all, pad = pad_pca)
  
  ## ---- build 6 panels
  p1_pca <- panel_plot_interval(qm_pca, "IPCA(QM)", lims_pca_square)
  p2_pca <- panel_plot_interval(cr_pca, "IPCA(CR)", lims_pca_square)
  p3_pca <- panel_plot_interval(imds,   "IMDS",     lims_pca_square)
  
  p1_tsne <- panel_plot_interval(qm_tsne,            "I-tSNE(QM)",  square_lims(get_bounds(qm_tsne),            pad_tsne))
  p2_tsne <- panel_plot_interval(cr1_tsne_embedding, "I-tSNE(CR)",  square_lims(get_bounds(cr1_tsne_embedding), pad_tsne))
  p3_tsne <- panel_plot_interval(lU2_tsne_embedding, "I-tSNE(MM)",  square_lims(get_bounds(lU2_tsne_embedding), pad_tsne))
  
  ## ---- legend grob
  carrier <- ggplot(qm_pca, aes(colour = factor(Group, levels = all_groups))) +
    geom_rect(aes(
      xmin = if ("Dim1_Lower" %in% names(qm_pca)) Dim1_Lower else Lower_X,
      xmax = if ("Dim1_Upper" %in% names(qm_pca)) Dim1_Upper else Upper_X,
      ymin = if ("Dim2_Lower" %in% names(qm_pca)) Dim2_Lower else Lower_Y,
      ymax = if ("Dim2_Upper" %in% names(qm_pca)) Dim2_Upper else Upper_Y
    ), fill = NA, linewidth = 0.7) +
    scale_colour_manual(values = cols, drop = FALSE, name = "Group") +
    theme_void() +
    theme(
      legend.position = "bottom",
      legend.title = element_text(size = 12, face = "bold"),
      legend.text  = element_text(size = 10)
    ) +
    guides(colour = guide_legend(
      nrow = 1,
      ncol = legend_ncol,
      byrow = TRUE,
      override.aes = list(linewidth = 1)
    ))
  
  carrier_grob <- ggplotGrob(carrier)
  leg <- gtable_filter(carrier_grob, "guide-box")
  
  ## ---- layout: 3 x 2
  row1 <- arrangeGrob(p1_pca, p2_pca, p3_pca, ncol = 3)
  row2 <- arrangeGrob(p1_tsne, p2_tsne, p3_tsne, ncol = 3)
  body <- arrangeGrob(row1, row2, ncol = 1, heights = unit.c(unit(1, "null"), unit(1, "null")))
  
  final_grob <- arrangeGrob(
    body,
    leg,
    ncol = 1,
    heights = unit.c(unit(1, "null"), unit(legend_height_in, "in"))
  )
  
  if (isTRUE(print_plot)) {
    grid.newpage()
    grid.draw(final_grob)
  }
  
  ## ---- export
  pdf_w <- 3 * panel_in
  pdf_h <- 2 * panel_in + legend_height_in
  
  if (isTRUE(use_cairo) && exists("cairo_pdf", mode = "function")) {
    cairo_pdf(filename = out_path, width = pdf_w, height = pdf_h)
    grid.draw(final_grob)
    dev.off()
  } else {
    pdf(file = out_path, width = pdf_w, height = pdf_h)
    grid.draw(final_grob)
    dev.off()
  }
  
  message("Saved: ", out_path)
  invisible(final_grob)
}