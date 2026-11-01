smoke_stack <- function(df,
                        norm_by = "zscore",
                        set_x = FALSE,
                        set_y = FALSE,
                        bins = 20,
                        bin_size = FALSE,
                        min_bin_number = 100,
                        lims = FALSE,
                        facet = FALSE,
                        facet_within = TRUE,
                        fixed_coord = TRUE,
                        mid_set = FALSE,
                        save_file = "figures/coloredBar.svg",
                        save_data = FALSE,
                        xdim = 5,
                        ydim = 0.5,
                        viridis_pallet = "inferno",
                        call_icellate = FALSE,
                        fuse_icellate = TRUE,
                        number_icellate = 10,
                        icellate_fun = NULL) {
  if (!is.data.frame(df)) {
    stop("`df` must be a data frame.", call. = FALSE)
  }

  x_name <- if (identical(set_x, FALSE)) {
    readline(prompt = "Which variable should be mapped on the x-axis? ")
  } else {
    as.character(set_x)[[1L]]
  }
  metric_names <- if (identical(set_y, FALSE)) {
    metric_count <- as.integer(readline(
      prompt = paste0("How many variables should be mapped across ", x_name, "? ")
    ))
    if (is.na(metric_count) || metric_count < 1L) {
      stop("At least one y-axis variable is required.", call. = FALSE)
    }
    vapply(
      seq_len(metric_count),
      function(index) readline(prompt = paste0("What is variable number ", index, "? ")),
      character(1)
    )
  } else {
    as.character(set_y)
  }

  facet_name <- if (identical(facet, FALSE)) NULL else as.character(facet)[[1L]]
  required_columns <- unique(c(x_name, metric_names, facet_name))
  missing_columns <- setdiff(required_columns, names(df))
  if (length(missing_columns) > 0L) {
    stop("Missing columns: ", paste(missing_columns, collapse = ", "), call. = FALSE)
  }
  if (!is.numeric(df[[x_name]])) {
    stop("The x-axis column must be numeric.", call. = FALSE)
  }
  if (!norm_by %in% c("zscore", "lq", "uq", "mean", "iqr", "none")) {
    stop("Unsupported `norm_by` value.", call. = FALSE)
  }

  interim <- df
  interim[["Name"]] <- interim[[x_name]]
  if (is.numeric(lims) && length(lims) == 2L) {
    keep <- interim[["Name"]] >= min(lims) & interim[["Name"]] <= max(lims)
    interim <- interim[keep, , drop = FALSE]
  }
  if (nrow(interim) == 0L) {
    stop("No rows remain after applying `lims`.", call. = FALSE)
  }

  x_range <- range(interim[["Name"]], na.rm = TRUE)
  if (!all(is.finite(x_range)) || x_range[[1L]] == x_range[[2L]]) {
    stop("The x-axis column must contain at least two distinct finite values.", call. = FALSE)
  }
  if (is.numeric(bin_size) && length(bin_size) == 1L && is.finite(bin_size) && bin_size > 0) {
    breaks <- seq(x_range[[1L]], x_range[[2L]], by = bin_size)
    if (breaks[[length(breaks)]] < x_range[[2L]]) {
      breaks <- c(breaks, x_range[[2L]])
    }
  } else {
    if (length(bins) != 1L || !is.numeric(bins) || bins < 1L) {
      stop("`bins` must be a positive integer.", call. = FALSE)
    }
    breaks <- seq(x_range[[1L]], x_range[[2L]], length.out = as.integer(bins) + 1L)
  }
  bin_count <- length(breaks) - 1L
  interim[["bin"]] <- as.integer(cut(
    interim[["Name"]],
    breaks = breaks,
    include.lowest = TRUE,
    labels = FALSE
  ))

  legend_name <- switch(
    norm_by,
    zscore = "Z-Score",
    lq = "Lower quartile",
    uq = "Upper quartile",
    mean = "% of Mean",
    iqr = "Interquartile range",
    none = "Intensity"
  )

  summarize_values <- function(values, normalized_values = NULL) {
    if (length(values) == 0L || all(is.na(values))) {
      return(NA_real_)
    }
    switch(
      norm_by,
      zscore = stats::median(normalized_values, na.rm = TRUE),
      lq = as.numeric(stats::quantile(values, 0.25, na.rm = TRUE, names = FALSE)),
      uq = as.numeric(stats::quantile(values, 0.75, na.rm = TRUE, names = FALSE)),
      mean = stats::median(normalized_values, na.rm = TRUE),
      iqr = stats::IQR(values, na.rm = TRUE),
      none = stats::median(values, na.rm = TRUE)
    )
  }

  facet_values <- if (is.null(facet_name)) "all" else unique(as.character(interim[[facet_name]]))
  result_parts <- vector("list", length(facet_values))

  for (facet_index in seq_along(facet_values)) {
    facet_value <- facet_values[[facet_index]]
    facet_rows <- if (is.null(facet_name)) {
      rep(TRUE, nrow(interim))
    } else {
      as.character(interim[[facet_name]]) == facet_value
    }
    facet_data <- interim[facet_rows, , drop = FALSE]
    facet_result <- data.frame(Bin = seq_len(bin_count))

    for (metric_name in metric_names) {
      metric_values <- as.numeric(facet_data[[metric_name]])
      reference_values <- if (isTRUE(facet_within) || is.null(facet_name)) {
        metric_values
      } else {
        as.numeric(interim[[metric_name]])
      }
      normalized <- NULL
      if (norm_by == "zscore") {
        center <- mean(reference_values, na.rm = TRUE)
        spread <- stats::sd(reference_values, na.rm = TRUE)
        normalized <- if (is.finite(spread) && spread > 0) {
          (metric_values - center) / spread
        } else {
          rep(NA_real_, length(metric_values))
        }
      } else if (norm_by == "mean") {
        center <- mean(reference_values, na.rm = TRUE)
        normalized <- if (is.finite(center) && center != 0) {
          100 * metric_values / center
        } else {
          rep(NA_real_, length(metric_values))
        }
      }

      summaries <- vapply(seq_len(bin_count), function(bin_id) {
        in_bin <- facet_data[["bin"]] == bin_id
        if (sum(in_bin, na.rm = TRUE) < min_bin_number) {
          return(NA_real_)
        }
        summarize_values(metric_values[in_bin], normalized[in_bin])
      }, numeric(1))
      facet_result[[metric_name]] <- summaries
    }

    if (!is.null(facet_name)) {
      facet_result[[facet_name]] <- facet_value
    }
    result_parts[[facet_index]] <- facet_result

    if (isTRUE(call_icellate)) {
      if (!is.function(icellate_fun)) {
        stop("`icellate_fun` must be a function when `call_icellate = TRUE`.", call. = FALSE)
      }
      for (bin_id in seq_len(bin_count)) {
        target_cells <- facet_data[facet_data[["bin"]] == bin_id, , drop = FALSE]
        if (nrow(target_cells) == 0L) {
          next
        }
        folder_prefix <- if (is.null(facet_name)) {
          ""
        } else {
          paste0(facet_name, "-", facet_value, "_")
        }
        icellate_fun(
          target_cells = target_cells,
          folder_name = paste0(folder_prefix, "Bin_", bin_id),
          fuse = fuse_icellate,
          sampling_number = min(number_icellate, nrow(target_cells)),
          line_analyses = FALSE,
          heat_intensity = TRUE
        )
      }
    }
  }

  wide_data <- do.call(rbind, result_parts)
  long_data <- tidyr::pivot_longer(
    wide_data,
    cols = tidyselect::all_of(metric_names),
    names_to = "Target",
    values_to = "Value"
  )
  if (!is.null(facet_name)) {
    long_data[["Target"]] <- paste0(
      facet_name,
      ": ",
      long_data[[facet_name]],
      " - ",
      long_data[["Target"]]
    )
  }

  bin_centers <- breaks[-length(breaks)] + diff(breaks) / 2
  label_indices <- seq.int(1L, bin_count, by = 2L)
  draft <- ggplot2::ggplot(
    long_data,
    ggplot2::aes(
      x = rlang::.data[["Bin"]],
      y = rlang::.data[["Target"]],
      fill = rlang::.data[["Value"]]
    )
  ) +
    ggplot2::geom_tile() +
    ggpubr::theme_classic2() +
    ggplot2::xlab(x_name) +
    ggplot2::ylab("") +
    ggplot2::guides(fill = ggplot2::guide_legend(title = legend_name)) +
    ggplot2::theme(legend.position = "bottom") +
    ggplot2::scale_x_continuous(
      breaks = label_indices,
      labels = round(bin_centers[label_indices], 1)
    )

  if (is.numeric(mid_set) && length(mid_set) == 1L) {
    draft <- draft + ggplot2::scale_fill_gradient2(
      low = "blue",
      midpoint = mid_set,
      mid = "black",
      high = "red"
    )
  } else {
    draft <- draft + viridis::scale_fill_viridis(option = viridis_pallet)
  }
  if (isTRUE(fixed_coord)) {
    draft <- draft + ggplot2::coord_fixed()
  }

  if (!identical(save_file, FALSE) && is.character(save_file) && length(save_file) == 1L) {
    dir.create(dirname(save_file), recursive = TRUE, showWarnings = FALSE)
    ggplot2::ggsave(filename = save_file, plot = draft, width = xdim, height = ydim)
  }
  if (!identical(save_data, FALSE) && is.character(save_data) && length(save_data) == 1L) {
    dir.create(dirname(save_data), recursive = TRUE, showWarnings = FALSE)
    utils::write.csv(long_data, save_data, row.names = FALSE)
  }

  invisible(list(
    data = long_data,
    plot = draft,
    binned_data = interim,
    breaks = breaks
  ))
}
