ergodic_slicing <- function(df,
                            subset_df = NULL,
                            variable = "dna_norm",
                            bins = 10,
                            bin_add_name = "",
                            rate = 24,
                            rate_description = "doubling time (hours)",
                            icell = FALSE,
                            icell_number = 3,
                            icell_fusion = FALSE,
                            pca = FALSE,
                            pca_groups = FALSE,
                            hist_save = "figures/ergodicHist.pdf",
                            ergo_save = "figures/ergodicRates.pdf",
                            ergo_file_save = "data/ergo.csv",
                            icellate_fun = NULL,
                            pca_fun = NULL) {
  if (!is.data.frame(df)) {
    stop("`df` must be a data frame.", call. = FALSE)
  }
  if (is.null(subset_df)) {
    if (!"edu" %in% names(df)) {
      stop("`subset_df` is required when `df` has no `edu` column.", call. = FALSE)
    }
    subset_df <- df[df[["edu"]] == "Positive", , drop = FALSE]
  }
  if (!is.data.frame(subset_df) || nrow(subset_df) == 0L) {
    stop("`subset_df` must be a non-empty data frame.", call. = FALSE)
  }
  if (!variable %in% names(subset_df)) {
    stop("`variable` was not found in `subset_df`.", call. = FALSE)
  }
  if (length(bins) != 1L || !is.numeric(bins) || bins < 1 || bins %% 1 != 0) {
    stop("`bins` must be a positive integer.", call. = FALSE)
  }

  bins <- as.integer(bins)
  variable_values <- subset_df[[variable]]
  if (!is.numeric(variable_values) || any(!is.finite(variable_values))) {
    stop("The slicing variable must contain only finite numeric values.", call. = FALSE)
  }
  if (length(unique(variable_values)) < 2L) {
    stop("The slicing variable must contain at least two distinct values.", call. = FALSE)
  }

  bin_factor <- cut(variable_values, breaks = bins, include.lowest = TRUE)
  list_df <- split(subset_df, bin_factor, drop = FALSE)

  if (isTRUE(icell)) {
    if (!is.function(icellate_fun)) {
      stop("`icellate_fun` must be a function when `icell = TRUE`.", call. = FALSE)
    }
    for (bin_index in seq_along(list_df)) {
      bin_set <- list_df[[bin_index]]
      if (nrow(bin_set) > 0L) {
        icellate_fun(
          target_cells = bin_set,
          folder_name = paste0(bin_add_name, "bin ", bin_index),
          verify_size = FALSE,
          fuse = icell_fusion,
          randomize = TRUE,
          sampling_number = min(nrow(bin_set), icell_number),
          line_analyses = FALSE
        )
      }
    }
  }

  histogram_plot <- ggplot2::ggplot(
    subset_df,
    ggplot2::aes(x = rlang::.data[[variable]])
  ) +
    ggplot2::geom_density() +
    ggplot2::xlab(variable) +
    ggplot2::theme_classic()

  if (isTRUE(pca) && !is.function(pca_fun)) {
    stop("Supply `pca_fun` when `pca = TRUE`.", call. = FALSE)
  }

  ergo_df <- data.frame(
    bin_id = seq_len(bins),
    bin_rate = NA_real_,
    volume = integer(bins),
    mean_stat = NA_real_
  )
  subset_fraction <- nrow(subset_df) / nrow(df)

  for (bin_index in seq_along(list_df)) {
    bin_set <- list_df[[bin_index]]
    if (isTRUE(pca) && nrow(bin_set) > 0L) {
      pca_path <- file.path("figures", "pcas", paste0(bin_add_name, "bin_", bin_index, ".png"))
      dir.create(dirname(pca_path), recursive = TRUE, showWarnings = FALSE)
      pca_fun(df = bin_set, groupz = pca_groups, save_file = pca_path)
    }

    ergo_df[["volume"]][bin_index] <- nrow(bin_set)
    if (nrow(bin_set) > 0L) {
      ergo_df[["bin_rate"]][bin_index] <- rate *
        ((2 - subset_fraction) / (nrow(bin_set) / nrow(subset_df)))
      ergo_df[["mean_stat"]][bin_index] <- stats::median(bin_set[[variable]])
    }
  }

  finite_rates <- is.finite(ergo_df[["bin_rate"]])
  rate_mean <- mean(ergo_df[["bin_rate"]][finite_rates], na.rm = TRUE)
  if (is.finite(rate_mean) && rate_mean != 0) {
    ergo_df[["bin_rate"]] <- ergo_df[["bin_rate"]] / rate_mean
  }

  rate_plot_data <- ergo_df[is.finite(ergo_df[["bin_rate"]]), , drop = FALSE]
  rate_plot <- ggplot2::ggplot(
    rate_plot_data,
    ggplot2::aes(
      x = rlang::.data[["bin_id"]],
      y = rlang::.data[["bin_rate"]],
      fill = rlang::.data[["volume"]]
    )
  ) +
    ggplot2::geom_col(position = ggplot2::position_dodge(width = 0.7)) +
    ggplot2::ylab(paste0("Mean bin progression rate (", rate_description, ")")) +
    ggplot2::xlab(paste0("Median size of ", variable)) +
    ggplot2::theme_classic() +
    viridis::scale_fill_viridis() +
    ggplot2::theme(legend.position = "top")

  if (!identical(hist_save, FALSE) && is.character(hist_save) && length(hist_save) == 1L) {
    dir.create(dirname(hist_save), recursive = TRUE, showWarnings = FALSE)
    ggplot2::ggsave(filename = hist_save, plot = histogram_plot)
  }
  if (!identical(ergo_save, FALSE) && is.character(ergo_save) && length(ergo_save) == 1L) {
    dir.create(dirname(ergo_save), recursive = TRUE, showWarnings = FALSE)
    ggplot2::ggsave(filename = ergo_save, plot = rate_plot)
  }
  if (
    !identical(ergo_file_save, FALSE) &&
      is.character(ergo_file_save) &&
      length(ergo_file_save) == 1L
  ) {
    dir.create(dirname(ergo_file_save), recursive = TRUE, showWarnings = FALSE)
    utils::write.csv(ergo_df, file = ergo_file_save, row.names = FALSE)
  }

  invisible(list(
    data = ergo_df,
    bins = list_df,
    histogram = histogram_plot,
    rates_plot = rate_plot
  ))
}
