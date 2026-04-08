run_method_suite <- function(data, id_col, settings, include_baselines = FALSE, vm_data = NULL) {
  vm_input <- if (is.null(vm_data)) data else vm_data
  results <- list(
    itsne_vm = do.call(itsne_vm, c(list(interval_data = vm_input, id_col = id_col), settings$itsne_vm)),
    itsne_qm = do.call(itsne_qm, c(list(interval_data = data, id_col = id_col), settings$itsne_qm)),
    itsne_mm = do.call(itsne_mm, c(list(interval_data = data, id_col = id_col), settings$itsne_mm))$embedding,
    itsne_cr = do.call(itsne_cr, c(list(interval_data = data, id_col = id_col), settings$itsne_cr))$embedding
  )

  if (include_baselines) {
    results <- c(
      list(
        ipca_vm = do.call(ipca_vm, c(list(interval_data = vm_input, id_col = id_col), settings$ipca_vm)),
        ipca_qm = do.call(ipca_qm, c(list(interval_data = data, id_col = id_col), settings$ipca_qm)),
        ipca_cr = do.call(ipca_cr, c(list(interval_data = data, id_col = id_col), settings$ipca_cr))
      ),
      if (isTRUE(settings$include_imds)) {
        list(imds = do.call(imds_box, c(list(interval_data = data, id_col = id_col), settings$imds)))
      } else {
        list()
      },
      results
    )
  }

  results
}

write_lcmc_tables <- function(lcmc_tables, output_dir, prefix) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  for (name in names(lcmc_tables)) {
    utils::write.csv(
      lcmc_tables[[name]],
      file = file.path(output_dir, paste0(prefix, "_", tolower(name), ".csv")),
      row.names = FALSE
    )
  }
}

save_interval_plot <- function(plot_object, path, width = 10, height = 6) {
  ggplot2::ggsave(path, plot = plot_object, width = width, height = height, units = "in")
}

simulation_settings <- function() {
  list(
    point_high_signal = list(
      generator = function() simulate_point_aggregation(matrix(c(6, 0, -6, 0, 0, 6), ncol = 2, byrow = TRUE), noise_dims = 0),
      params = list(
        itsne_vm = list(perplexity = 30, eta = 200),
        itsne_qm = list(m = 4, perplexity = 30, eta = 5),
        itsne_mm = list(perplexity = 5, alpha = 0.5, penalty_lambda = 0.2, learning_rate = 200, initial_P_gain = 12),
        itsne_cr = list(perplexity = 5, lambda = 1, eta = 5, EE = 12, momentum = 0.1, final_momentum = 0.1)
      )
    ),
    point_high_noise = list(
      generator = function() simulate_point_aggregation(matrix(c(6, 0, -6, 0, 0, 6), ncol = 2, byrow = TRUE), noise_dims = 3),
      params = list(
        itsne_vm = list(perplexity = 300, eta = 200),
        itsne_qm = list(m = 5, perplexity = 30, eta = 5),
        itsne_mm = list(perplexity = 5, alpha = 0.5, penalty_lambda = 0.1, learning_rate = 200, initial_P_gain = 12),
        itsne_cr = list(perplexity = 5, lambda = 1, eta = 5, EE = 12, momentum = 0.1, final_momentum = 0.1)
      )
    ),
    point_low_signal = list(
      generator = function() simulate_point_aggregation(matrix(c(0, 1, -6, 0, 0, 6), ncol = 2, byrow = TRUE), noise_dims = 0),
      params = list(
        itsne_vm = list(perplexity = 30, eta = 200),
        itsne_qm = list(m = 4, perplexity = 30, eta = 5),
        itsne_mm = list(perplexity = 5, alpha = 0.5, penalty_lambda = 0.2, learning_rate = 200, initial_P_gain = 12),
        itsne_cr = list(perplexity = 5, lambda = 1, eta = 5, EE = 12, momentum = 0.1, final_momentum = 0.1)
      )
    ),
    point_low_noise = list(
      generator = function() simulate_point_aggregation(matrix(c(0, 1, -6, 0, 0, 6), ncol = 2, byrow = TRUE), noise_dims = 3),
      params = list(
        itsne_vm = list(perplexity = 200, eta = 200),
        itsne_qm = list(m = 9, perplexity = 20, eta = 5),
        itsne_mm = list(perplexity = 5, alpha = 0.5, penalty_lambda = 0.2, learning_rate = 200, initial_P_gain = 12),
        itsne_cr = list(perplexity = 5, lambda = 1, eta = 5, EE = 4, momentum = 0.1, final_momentum = 0.3)
      )
    ),
    direct_high_signal = list(
      generator = function() simulate_direct_intervals(matrix(c(-2, 0, 0, 2, 2, -2), ncol = 2, byrow = TRUE), noise_dims = 0),
      params = list(
        itsne_vm = list(perplexity = 30, eta = 200),
        itsne_qm = list(m = 5, perplexity = 40, eta = 200),
        itsne_mm = list(perplexity = 5, alpha = 0.5, penalty_lambda = 0.2, learning_rate = 200, initial_P_gain = 12),
        itsne_cr = list(perplexity = 5, lambda = 1, eta = 5, EE = 12)
      )
    ),
    direct_high_noise = list(
      generator = function() simulate_direct_intervals(matrix(c(-2, 0, 0, 2, 2, -2), ncol = 2, byrow = TRUE), noise_dims = 3),
      params = list(
        itsne_vm = list(perplexity = 50, eta = 200),
        itsne_qm = list(m = 5, perplexity = 30, eta = 200),
        itsne_mm = list(perplexity = 5, alpha = 0.5, penalty_lambda = 0.1, learning_rate = 200, initial_P_gain = 12),
        itsne_cr = list(perplexity = 5, lambda = 1, eta = 5, EE = 12)
      )
    ),
    direct_low_signal = list(
      generator = function() simulate_direct_intervals(matrix(c(-1.5, 0, 0, 1.5, 0, -1.5), ncol = 2, byrow = TRUE), noise_dims = 0),
      params = list(
        itsne_vm = list(perplexity = 30, eta = 200),
        itsne_qm = list(m = 5, perplexity = 40, eta = 200),
        itsne_mm = list(perplexity = 5, alpha = 0.5, penalty_lambda = 0.15, learning_rate = 200, initial_P_gain = 12),
        itsne_cr = list(perplexity = 5, lambda = 1, eta = 5, EE = 12, momentum = 0.1, final_momentum = 0.1)
      )
    ),
    direct_low_noise = list(
      generator = function() simulate_direct_intervals(matrix(c(-1.5, 0, 0, 1.5, 0, -1.5), ncol = 2, byrow = TRUE), noise_dims = 3),
      params = list(
        itsne_vm = list(perplexity = 300, eta = 2000),
        itsne_qm = list(m = 7, perplexity = 50, eta = 200),
        itsne_mm = list(perplexity = 5, alpha = 0.5, penalty_lambda = 0.06, learning_rate = 100, initial_P_gain = 12),
        itsne_cr = list(perplexity = 10, lambda = 1, eta = 5, EE = 12, momentum = 0.5, final_momentum = 0.5)
      )
    )
  )
}

#' Run the manuscript simulation study
#'
#' Reproduce the manuscript simulation workflows, including interval-data
#' generation, model fitting, modified-LCMC evaluation, and figure and table
#' export.
#'
#' @param output_dir Directory where results and figures will be saved.
#' @param seed Random seed. Included for interface consistency; the scenario
#'   generators currently manage their own seeds.
#'
#' @return A named list containing the simulated datasets, fitted method outputs,
#'   and modified-LCMC tables for each simulation setting.
#'
#' @examples
#' if (interactive()) {
#'   out <- run_simulation_studies(output_dir = tempdir())
#'   names(out)
#' }
#' @export
run_simulation_studies <- function(output_dir = "results/simulations", seed = 12345) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  settings <- simulation_settings()
  results <- list()

  for (scenario in names(settings)) {
    spec <- settings[[scenario]]
    dataset <- spec$generator()
    methods <- run_method_suite(dataset, id_col = "Group", settings = spec$params)
    lcmc <- compute_lcmc_tables(dataset, methods, k_range = 1:26)

    save_interval_plot(
      plot_interval_method_grid(
        methods,
        titles = c("I-tSNE(VM)", "I-tSNE(QM)", "I-tSNE(MM)", "I-tSNE(CR)"),
        ncol = 4
      ),
      file.path(output_dir, paste0(scenario, "_methods.pdf")),
      width = 10,
      height = 4
    )

    save_interval_plot(
      plot_lcmc_grid(lcmc, method_name_map = default_method_name_map()),
      file.path(output_dir, paste0(scenario, "_lcmc.pdf")),
      width = 10,
      height = 4
    )

    utils::write.csv(dataset, file.path(output_dir, paste0(scenario, "_data.csv")), row.names = FALSE)
    write_lcmc_tables(lcmc, output_dir, scenario)

    results[[scenario]] <- list(
      data = dataset,
      methods = methods,
      lcmc = lcmc
    )
  }

  results
}