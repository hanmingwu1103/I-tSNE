read_example_data <- function(filename) {
  installed_path <- system.file("extdata", filename, package = "itsne")
  source_path <- file.path("inst", "extdata", filename)
  data_path <- if (nzchar(installed_path)) installed_path else source_path
  readr::read_csv(data_path, show_col_types = FALSE)
}

#' Run the manuscript real-data analyses
#'
#' Reproduce the manuscript real-data workflows for the Face and Digits datasets,
#' including preprocessing, model fitting, modified-LCMC evaluation, and export of
#' figures and summary tables.
#'
#' @param output_dir Directory where results and figures will be saved.
#'
#' @return A named list with two components, `face` and `digits`, each containing
#'   the prepared data, fitted method outputs, and modified-LCMC tables.
#'
#' @examples
#' if (interactive()) {
#'   out <- run_real_data_analysis(output_dir = tempdir())
#'   names(out)
#' }
#' @export
run_real_data_analysis <- function(output_dir = "results/real-data") {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  face_raw <- read_example_data("facedata.csv")
  face_data <- face_raw[, c("species", grep("_(Lower|Upper)$", names(face_raw), value = TRUE))]
  face_data <- standardize_interval_data(face_data, id_col = "species", method = 1)

  face_settings <- list(
    ipca_vm = list(dims = 2),
    ipca_qm = list(m = 5, dims = 2),
    ipca_cr = list(dims = 2),
    itsne_vm = list(perplexity = 50, eta = 200),
    itsne_qm = list(m = 5, perplexity = 50, eta = 200),
    itsne_mm = list(perplexity = 5, alpha = 0.5, penalty_lambda = 0.07, learning_rate = 5, initial_P_gain = 4),
    itsne_cr = list(perplexity = 5, lambda = 1, eta = 5, EE = 1)
  )

  face_methods <- run_method_suite(face_data, id_col = "species", settings = face_settings, include_baselines = TRUE)
  face_lcmc <- compute_lcmc_tables(face_data, face_methods, k_range = 1:26)

  save_interval_plot(
    plot_interval_method_grid(
      face_methods,
      titles = unname(default_method_name_map()[names(face_methods)]),
      ncol = 4
    ),
    file.path(output_dir, "face_methods.pdf"),
    width = 12,
    height = 6
  )
  save_interval_plot(
    plot_lcmc_grid(face_lcmc, method_name_map = default_method_name_map()),
    file.path(output_dir, "face_lcmc.pdf"),
    width = 10,
    height = 4
  )

  digits_main_raw <- read_example_data("digits_interval.csv")
  digits_vm_raw <- read_example_data("digits_interval_pca.csv")
  digits_main <- digits_main_raw[, c("label", grep("_(Lower|Upper)$", names(digits_main_raw), value = TRUE))]
  digits_vm <- digits_vm_raw[, c("label", grep("_(Lower|Upper)$", names(digits_vm_raw), value = TRUE))]
  digits_main <- standardize_interval_data(digits_main, id_col = "label", method = 3)
  digits_vm <- standardize_interval_data(digits_vm, id_col = "label", method = 3)

  digits_settings <- list(
    ipca_vm = list(dims = 2),
    ipca_qm = list(m = 5, dims = 2),
    ipca_cr = list(dims = 2),
    itsne_vm = list(perplexity = 400, eta = 200),
    itsne_qm = list(m = 5, perplexity = 100, eta = 200),
    itsne_mm = list(perplexity = 30, alpha = 0.5, penalty_lambda = 1.8, learning_rate = 5, initial_P_gain = 4),
    itsne_cr = list(perplexity = 50, lambda = 1, eta = 5, EE = 12)
  )

  digits_methods <- run_method_suite(
    digits_main,
    id_col = "label",
    settings = digits_settings,
    include_baselines = TRUE,
    vm_data = digits_vm
  )
  digits_lcmc <- compute_lcmc_tables(digits_main, digits_methods, k_range = 1:26)

  save_interval_plot(
    plot_interval_method_grid(
      digits_methods,
      titles = unname(default_method_name_map()[names(digits_methods)]),
      ncol = 4
    ),
    file.path(output_dir, "digits_methods.pdf"),
    width = 12,
    height = 6
  )
  save_interval_plot(
    plot_lcmc_grid(digits_lcmc, method_name_map = default_method_name_map()),
    file.path(output_dir, "digits_lcmc.pdf"),
    width = 10,
    height = 4
  )

  write_lcmc_tables(face_lcmc, output_dir, "face")
  write_lcmc_tables(digits_lcmc, output_dir, "digits")
  utils::write.csv(face_data, file.path(output_dir, "face_data.csv"), row.names = FALSE)
  utils::write.csv(digits_main, file.path(output_dir, "digits_data.csv"), row.names = FALSE)

  list(
    face = list(data = face_data, methods = face_methods, lcmc = face_lcmc),
    digits = list(data = digits_main, vm_data = digits_vm, methods = digits_methods, lcmc = digits_lcmc)
  )
}
