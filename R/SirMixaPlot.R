# Package-safe plotting and analysis helpers for MicroCyte.

.assert_data_columns <- function(data, columns) {
  if (!is.data.frame(data)) {
    stop("`data` must be a data frame.", call. = FALSE)
  }
  missing_columns <- setdiff(columns, names(data))
  if (length(missing_columns) > 0L) {
    stop(
      "Missing required columns: ",
      paste(missing_columns, collapse = ", "),
      call. = FALSE
    )
  }
  invisible(data)
}

.assert_numeric_columns <- function(data, columns) {
  .assert_data_columns(data, columns)
  numeric_columns <- vapply(data[columns], is.numeric, logical(1L))
  if (!all(numeric_columns)) {
    stop(
      "Columns must be numeric: ",
      paste(columns[!numeric_columns], collapse = ", "),
      call. = FALSE
    )
  }
  invisible(data)
}

.column_values <- function(data, column) {
  data[[column]]
}

#' Load and prepare a MicroCyte data file
#'
#' @param filo Path to a CSV file.
#' @param run Whether to create a scatter plot.
#' @param adjust Whether to calculate EdU and DNA normalization columns when possible.
#' @param x,y Column names used when `run` is `TRUE`.
#' @param x_name,y_name Optional axis labels.
#' @param dna_source Column to use for `log2_dna` when it is not already present.
#' @param output_file Optional output path. By default, a new `_cells.csv` file is
#'   created for an input that does not already have that suffix.
#' @param write_file Whether to write prepared data to `output_file`.
#' @param normalizer_args Named arguments passed to [normalizer()].
#' @return A list containing prepared `data`, an optional `plot`, normalization results,
#'   and the `output_file` path.
sirmixaplot <- function(filo, run = TRUE, adjust = TRUE, x = NULL, y = NULL,
                        x_name = NULL, y_name = NULL,
                        dna_source = "IntDen_ANC_anchor", output_file = NULL,
                        write_file = TRUE, normalizer_args = list()) {
  if (!is.character(filo) || length(filo) != 1L || !nzchar(filo)) {
    stop("`filo` must be one non-empty file path.", call. = FALSE)
  }
  if (!is.list(normalizer_args)) {
    stop("`normalizer_args` must be a list.", call. = FALSE)
  }

  data <- readr::read_csv(filo, show_col_types = FALSE)
  is_cells_file <- grepl("_cells[.]csv$", filo)
  if (is.null(output_file)) {
    output_file <- if (is_cells_file) filo else sub("[.]csv$", "_cells.csv", filo)
  }

  if (!"log2_dna" %in% names(data) && dna_source %in% names(data)) {
    data$log2_dna <- log2(data[[dna_source]])
  }
  mean_columns <- grep("Mean_", names(data), value = TRUE)
  for (column in mean_columns) {
    log_column <- paste0("L", column)
    if (!log_column %in% names(data)) {
      data[[log_column]] <- log10(data[[column]])
    }
  }

  normalization <- NULL
  can_normalize <- all(c("LMean_ANC_edu", "log2_dna") %in% names(data))
  if (isTRUE(adjust) && !"dna_norm" %in% names(data) && can_normalize) {
    normalization <- do.call(normalizer, c(list(data = data), normalizer_args))
    data <- normalization$data
  }

  plot <- NULL
  if (isTRUE(run)) {
    if (is.null(x) || is.null(y)) {
      stop("`x` and `y` are required when `run` is TRUE.", call. = FALSE)
    }
    plot <- grapho(data, x, y, x_name = x_name, y_name = y_name)
  }
  if (isTRUE(write_file)) {
    readr::write_csv(data, output_file)
  }

  list(
    data = data,
    plot = plot,
    normalization = normalization,
    output_file = output_file
  )
}

#' Create a scatter plot from explicitly selected columns
#'
#' @return The plot returned by [grapho()].
joiner <- function(df, x, y, x_name = NULL, y_name = NULL) {
  grapho(df, x, y, x_name = x_name, y_name = y_name)
}

#' Build a MicroCyte scatter plot
#'
#' @return A `ggplot` object.
grapho <- function(df, x, y, x_name = NULL, y_name = NULL) {
  .assert_numeric_columns(df, c(x, y))
  if (is.null(x_name)) {
    x_name <- x
  }
  if (is.null(y_name)) {
    y_name <- y
  }

  ggplot2::ggplot(
    df,
    ggplot2::aes(x = rlang::.data[[x]], y = rlang::.data[[y]])
  ) +
    ggplot2::geom_point(size = 2) +
    ggplot2::scale_x_continuous(name = x_name) +
    ggplot2::scale_y_continuous(name = y_name) +
    ggplot2::theme_linedraw() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = 40, face = "bold"),
      axis.title = ggplot2::element_text(size = 30, face = "bold"),
      axis.text = ggplot2::element_text(face = "bold", size = 16),
      axis.ticks = ggplot2::element_blank()
    )
}

#' Add quadrant summaries to a scatter plot
#'
#' @param df Data frame containing the plotted values.
#' @param x,y Numeric column names.
#' @param xi,yi Quadrant intercepts.
#' @param plot Optional base plot. A plot is created with [grapho()] when omitted.
#' @param save_file Optional path at which to save the quadrant plot.
#' @param graphit Whether to include density contours.
#' @return A list with `counts`, `percentages`, and `plot`.
modthequad <- function(df, x, y, xi = 1.5, yi = 1.5, plot = NULL,
                       save_file = NULL, graphit = TRUE) {
  .assert_numeric_columns(df, c(x, y))
  xi <- as.numeric(xi)
  yi <- as.numeric(yi)
  if (length(xi) != 1L || length(yi) != 1L || !is.finite(xi) || !is.finite(yi)) {
    stop("`xi` and `yi` must be finite numeric scalars.", call. = FALSE)
  }

  x_values <- df[[x]]
  y_values <- df[[y]]
  valid <- is.finite(x_values) & is.finite(y_values)
  counts <- c(
    Q1 = sum(valid & x_values < xi & y_values > yi),
    Q2 = sum(valid & x_values > xi & y_values > yi),
    Q3 = sum(valid & x_values < xi & y_values < yi),
    Q4 = sum(valid & x_values > xi & y_values < yi)
  )
  total <- sum(valid)
  percentages <- if (total == 0L) rep(NA_real_, 4L) else counts / total * 100
  summary <- data.frame(
    quadrant = names(counts),
    count = unname(counts),
    percent = unname(percentages)
  )

  if (is.null(plot)) {
    plot <- grapho(df, x, y)
  }
  ranges <- list(
    x = range(x_values[valid], na.rm = TRUE),
    y = range(y_values[valid], na.rm = TRUE)
  )
  if (total == 0L) {
    stop("No finite x/y pairs are available for quadrant analysis.", call. = FALSE)
  }
  labels <- paste0(names(counts), ": ", format(round(percentages, 2), nsmall = 2), "%")
  quadrant_plot <- plot +
    ggplot2::geom_hline(
      yintercept = yi,
      color = "red",
      linetype = "dashed",
      linewidth = 2
    ) +
    ggplot2::geom_vline(
      xintercept = xi,
      color = "red",
      linetype = "dashed",
      linewidth = 2
    ) +
    ggplot2::annotate("text", x = ranges$x[1L], y = ranges$y[2L],
                      label = labels[1L], hjust = 0, vjust = 1) +
    ggplot2::annotate("text", x = ranges$x[2L], y = ranges$y[2L],
                      label = labels[2L], hjust = 1, vjust = 1) +
    ggplot2::annotate("text", x = ranges$x[1L], y = ranges$y[1L],
                      label = labels[3L], hjust = 0, vjust = 0) +
    ggplot2::annotate("text", x = ranges$x[2L], y = ranges$y[1L],
                      label = labels[4L], hjust = 1, vjust = 0)
  if (isTRUE(graphit)) {
    quadrant_plot <- quadrant_plot + ggplot2::geom_density_2d()
  }
  if (!is.null(save_file)) {
    ggplot2::ggsave(save_file, plot = quadrant_plot)
  }

  list(counts = counts, percentages = summary, plot = quadrant_plot)
}

#' Remove observations on one side of a threshold
#'
#' @param df Input data frame.
#' @param parameter Numeric column to evaluate.
#' @param direction One of `"b"`, `"s"`, or `"="`. `"b"` removes values bigger
#'   than `value`; `"s"` removes values smaller than `value`.
#' @param value Numeric threshold.
#' @return The filtered data frame.
cull <- function(df, parameter, direction, value) {
  .assert_numeric_columns(df, parameter)
  value <- as.numeric(value)
  if (length(value) != 1L || !is.finite(value)) {
    stop("`value` must be one finite number.", call. = FALSE)
  }
  values <- df[[parameter]]
  keep <- switch(
    direction,
    b = values < value,
    s = values > value,
    `=` = values != value,
    stop("`direction` must be one of 'b', 's', or '='.", call. = FALSE)
  )
  df[!is.na(keep) & keep, , drop = FALSE]
}

#' Estimate point density in two dimensions
#'
#' @return A numeric density value for each input pair.
get_density <- function(px, py, ...) {
  if (!is.numeric(px) || !is.numeric(py) || length(px) != length(py)) {
    stop("`px` and `py` must be numeric vectors of equal length.", call. = FALSE)
  }
  finite <- is.finite(px) & is.finite(py)
  result <- rep(NA_real_, length(px))
  if (sum(finite) < 2L) {
    return(result)
  }
  density <- MASS::kde2d(px[finite], py[finite], ...)
  x_index <- findInterval(px[finite], density$x, all.inside = TRUE)
  y_index <- findInterval(py[finite], density$y, all.inside = TRUE)
  result[finite] <- density$z[cbind(x_index, y_index)]
  result
}

#' Color a scatter plot by local point density
#'
#' @return A list containing data with a `density` column and a `ggplot` object.
color_density <- function(df, x, y, x_name = NULL, y_name = NULL, n = 100L, ...) {
  .assert_numeric_columns(df, c(x, y))
  if (is.null(x_name)) {
    x_name <- x
  }
  if (is.null(y_name)) {
    y_name <- y
  }
  result <- df
  result$density <- get_density(result[[x]], result[[y]], n = n, ...)
  plot <- ggplot2::ggplot(
    result,
    ggplot2::aes(
      x = rlang::.data[[x]],
      y = rlang::.data[[y]],
      color = rlang::.data$density
    )
  ) +
    ggplot2::geom_point() +
    viridis::scale_color_viridis() +
    ggplot2::labs(x = x_name, y = y_name, color = "Density") +
    ggplot2::theme_bw(base_size = 16)
  list(data = result, plot = plot)
}

#' Normalize EdU and DNA measurements
#'
#' @param data Input data frame.
#' @param edu_column,dna_column Numeric input columns.
#' @param edu_cutoff Optional EdU-negative cutoff. Estimated from density when `NULL`.
#' @param diploid_peak Optional 2N DNA peak. Estimated from EdU-negative cells when `NULL`.
#' @param ploidy_cutoffs Two cutoffs on normalized DNA values for 4N and >4N calls.
#' @param remove_na Whether to remove rows with missing or infinite input values.
#' @return A list containing normalized `data`, diagnostic `plots`, and `thresholds`.
normalizer <- function(data, edu_column = "LMean_ANC_edu", dna_column = "log2_dna",
                       edu_cutoff = NULL, diploid_peak = NULL,
                       ploidy_cutoffs = c(1.6, 2.6), remove_na = TRUE) {
  .assert_numeric_columns(data, c(edu_column, dna_column))
  if (length(ploidy_cutoffs) != 2L || any(!is.finite(ploidy_cutoffs)) ||
        ploidy_cutoffs[1L] >= ploidy_cutoffs[2L]) {
    stop("`ploidy_cutoffs` must be two increasing finite numbers.", call. = FALSE)
  }

  result <- data
  valid <- is.finite(result[[edu_column]]) & is.finite(result[[dna_column]])
  if (isTRUE(remove_na)) {
    result <- result[valid, , drop = FALSE]
  } else if (!all(valid)) {
    stop("Input columns contain missing or infinite values.", call. = FALSE)
  }
  if (nrow(result) < 2L) {
    stop("At least two finite observations are required.", call. = FALSE)
  }

  edu_values <- result[[edu_column]]
  if (is.null(edu_cutoff)) {
    edu_density <- stats::density(edu_values)
    edu_cutoff <- edu_density$x[which.max(edu_density$y)] + 0.1
  }
  edu_cutoff <- as.numeric(edu_cutoff)
  if (length(edu_cutoff) != 1L || !is.finite(edu_cutoff)) {
    stop("`edu_cutoff` must be one finite number.", call. = FALSE)
  }

  edu_negative <- edu_values < edu_cutoff
  if (sum(edu_negative) < 2L) {
    stop("The EdU cutoff must identify at least two negative cells.", call. = FALSE)
  }
  negative_edu_values <- edu_values[edu_negative]
  if (is.null(diploid_peak)) {
    dna_density <- stats::density(result[[dna_column]][edu_negative])
    diploid_peak <- dna_density$x[which.max(dna_density$y)]
  }
  diploid_peak <- as.numeric(diploid_peak)
  if (length(diploid_peak) != 1L || !is.finite(diploid_peak)) {
    stop("`diploid_peak` must be one finite number.", call. = FALSE)
  }

  result$edu <- ifelse(edu_values > edu_cutoff, "Positive", "Negative")
  result$edu_norm <- (edu_values + 1) - mean(negative_edu_values)
  result$dna_norm <- (result[[dna_column]] + 1) - diploid_peak
  result$ploidy <- "2N"
  result$ploidy[result$dna_norm > ploidy_cutoffs[1L]] <- "4N"
  result$ploidy[result$dna_norm > ploidy_cutoffs[2L]] <- ">4N"

  edu_plot <- ggplot2::ggplot(
    result,
    ggplot2::aes(x = rlang::.data[[edu_column]])
  ) +
    ggplot2::geom_density() +
    ggplot2::geom_vline(xintercept = edu_cutoff) +
    ggplot2::labs(x = "Log EdU")
  dna_plot <- ggplot2::ggplot(
    result[edu_negative, , drop = FALSE],
    ggplot2::aes(x = rlang::.data[[dna_column]])
  ) +
    ggplot2::geom_density() +
    ggplot2::geom_vline(xintercept = diploid_peak, color = "green", linewidth = 1.5) +
    ggplot2::labs(x = "DNA content")

  list(
    data = result,
    plots = list(edu = edu_plot, dna = dna_plot),
    thresholds = list(
      edu_cutoff = edu_cutoff,
      diploid_peak = diploid_peak,
      ploidy_cutoffs = ploidy_cutoffs
    )
  )
}

#' Compensate one measurement for linear spectral overlap
#'
#' @param df Input data frame.
#' @param x,y Numeric columns. Compensation is applied to `y`.
#' @param point_one,point_two Numeric length-two vectors defining the overlap line.
#' @param recalculate Whether to create integrated-density and log columns for `y`.
#' @param area_column Area column used when `recalculate` is `TRUE`.
#' @return A list containing compensated `data`, fitted coefficients, and before/after plots.
compensate <- function(df, x, y, point_one, point_two, recalculate = FALSE,
                       area_column = "Area") {
  .assert_numeric_columns(df, c(x, y))
  if (length(point_one) != 2L || length(point_two) != 2L ||
        any(!is.finite(c(point_one, point_two)))) {
    stop("`point_one` and `point_two` must be finite x/y pairs.", call. = FALSE)
  }
  run <- point_two[1L] - point_one[1L]
  if (run == 0) {
    stop("Compensation points must have different x values.", call. = FALSE)
  }
  slope <- (point_two[2L] - point_one[2L]) / run
  intercept <- point_one[2L] - point_one[1L] * slope
  points <- data.frame(
    x = c(point_one[1L], point_two[1L]),
    y = c(point_one[2L], point_two[2L])
  )

  before_plot <- grapho(df, x, y) +
    ggplot2::geom_point(
      data = points,
      ggplot2::aes(x = rlang::.data$x, y = rlang::.data$y),
      color = "red",
      size = 4,
      inherit.aes = FALSE
    ) +
    ggplot2::geom_abline(slope = slope, intercept = intercept, linewidth = 2,
                         color = "blue")
  result <- df
  result[[y]] <- result[[y]] - (slope * result[[x]] + intercept)
  if (isTRUE(recalculate)) {
    .assert_numeric_columns(result, area_column)
    result[[paste0("I", y)]] <- result[[y]] * result[[area_column]]
    result[[paste0("L", y)]] <- log10(result[[y]] + 1)
  }
  after_plot <- grapho(result, x, y)

  list(
    data = result,
    coefficients = c(slope = slope, intercept = intercept),
    points = points,
    plots = list(before = before_plot, after = after_plot)
  )
}

#' Reconstruct cell locations as a point map
#'
#' @param df Input data frame.
#' @param highlight Either cell identifiers to highlight, a logical vector, a column name,
#'   or `FALSE`.
#' @param y,x,size Exact column names or prefixes matching exactly one column each.
#' @param id_column Identifier column used for numeric or character highlights.
#' @param colz Colors used for highlighted groups.
#' @return A `ggplot` object.
re_map <- function(df, highlight = FALSE, y = "X_ANC_", x = "Y_ANC_",
                   size = "Area_ANC_", x_name = "X position",
                   y_name = "Y position", id_column = "Number",
                   colz = c("#000000", "#E69F00", "#56B4E9", "#009E73",
                            "#F0E442", "#0072B2", "#D55E00", "#CC79A7")) {
  resolve_column <- function(prefix) {
    matches <- names(df)[startsWith(names(df), prefix)]
    if (length(matches) != 1L) {
      stop("`", prefix, "` must identify exactly one column.", call. = FALSE)
    }
    matches
  }
  x_column <- resolve_column(x)
  y_column <- resolve_column(y)
  size_column <- resolve_column(size)
  .assert_numeric_columns(df, c(x_column, y_column, size_column))

  plot_data <- df
  mapping <- ggplot2::aes(
    x = rlang::.data[[y_column]],
    y = -rlang::.data[[x_column]],
    size = rlang::.data[[size_column]]
  )
  use_color <- !identical(highlight, FALSE)
  if (use_color) {
    if (is.character(highlight) && length(highlight) == 1L && highlight %in% names(df)) {
      plot_data$.highlight <- as.factor(df[[highlight]])
    } else if (is.logical(highlight) && length(highlight) == nrow(df)) {
      plot_data$.highlight <- as.factor(highlight)
    } else {
      .assert_data_columns(df, id_column)
      plot_data$.highlight <- as.factor(df[[id_column]] %in% highlight)
    }
    mapping$colour <- rlang::expr(rlang::.data$.highlight)
  }

  plot <- ggplot2::ggplot(plot_data, mapping) +
    ggplot2::geom_point() +
    ggplot2::labs(x = x_name, y = y_name, color = "Highlight") +
    ggplot2::theme_classic() +
    ggplot2::theme(
      axis.text = ggplot2::element_text(face = "bold", size = 16),
      axis.ticks = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(size = 40, face = "bold"),
      axis.title = ggplot2::element_text(size = 30, face = "bold")
    )
  if (use_color) {
    plot <- plot + ggplot2::scale_color_manual(values = colz)
  }
  plot
}

#' Run principal component analysis or t-SNE
#'
#' @param df Input data frame.
#' @param plot_type Either `"pca"` or `"tsne"`.
#' @param save_file Optional image path.
#' @param subz Optional regular expression selecting numeric columns.
#' @param group_by Optional grouping column.
#' @param perplexity t-SNE perplexity.
#' @return A list with the model, coordinates, plot, and represented percentage.
run_ca <- function(df, plot_type = c("tsne", "pca"), save_file = NULL,
                   subz = NULL, group_by = NULL, perplexity = 30) {
  if (!is.data.frame(df)) {
    stop("`df` must be a data frame.", call. = FALSE)
  }
  plot_type <- match.arg(plot_type)
  if (!is.null(group_by)) {
    .assert_data_columns(df, group_by)
  }

  numeric_columns <- names(df)[vapply(df, is.numeric, logical(1L))]
  if (!is.null(subz)) {
    numeric_columns <- numeric_columns[grepl(subz, numeric_columns)]
  }
  if (length(numeric_columns) < 2L) {
    stop("At least two numeric analysis columns are required.", call. = FALSE)
  }
  analysis_data <- df[numeric_columns]
  analysis_data[] <- lapply(analysis_data, function(values) {
    values[!is.finite(values)] <- 0
    values
  })
  constant_columns <- foo(analysis_data)
  if (length(constant_columns) > 0L) {
    analysis_data <- analysis_data[, -constant_columns, drop = FALSE]
  }
  if (ncol(analysis_data) < 2L) {
    stop("At least two non-constant numeric columns are required.", call. = FALSE)
  }

  if (identical(plot_type, "pca")) {
    model <- stats::prcomp(analysis_data, center = TRUE, scale. = TRUE)
    coordinates <- data.frame(x = model$x[, 1L], y = model$x[, 2L])
    variance <- model$sdev^2 / sum(model$sdev^2) * 100
    axis_labels <- c(
      paste0("PC1 (", round(variance[1L], 2), " %)"),
      paste0("PC2 (", round(variance[2L], 2), " %)")
    )
  } else {
    if (nrow(analysis_data) < 2L) {
      stop("t-SNE requires at least two observations.", call. = FALSE)
    }
    max_perplexity <- (nrow(analysis_data) - 1) / 3
    if (!is.numeric(perplexity) || length(perplexity) != 1L || perplexity <= 0 ||
          perplexity >= max_perplexity) {
      stop("`perplexity` must be positive and less than (nrow(df) - 1) / 3.",
           call. = FALSE)
    }
    model <- Rtsne::Rtsne(as.matrix(analysis_data), perplexity = perplexity)
    coordinates <- data.frame(x = model$Y[, 1L], y = model$Y[, 2L])
    axis_labels <- c("t-SNE 1", "t-SNE 2")
  }

  if (!is.null(group_by)) {
    coordinates$group <- as.factor(df[[group_by]])
  }
  plot <- ggplot2::ggplot(
    coordinates,
    ggplot2::aes(x = rlang::.data$x, y = rlang::.data$y)
  )
  if (is.null(group_by)) {
    plot <- plot +
      ggplot2::geom_point(size = 3, alpha = 0.5) +
      ggplot2::geom_density_2d()
  } else {
    plot <- plot +
      ggplot2::geom_point(
        ggplot2::aes(color = rlang::.data$group),
        size = 3,
        alpha = 0.5
      ) +
      ggplot2::geom_density_2d(ggplot2::aes(color = rlang::.data$group))
  }
  plot <- plot +
    ggplot2::labs(x = axis_labels[1L], y = axis_labels[2L], color = group_by) +
    ggplot2::theme_classic()
  if (!is.null(save_file)) {
    ggplot2::ggsave(save_file, plot = plot, dpi = 300, height = 4, width = 4)
  }

  represented <- nrow(analysis_data) * ncol(analysis_data) / (nrow(df) * ncol(df)) * 100
  list(
    model = model,
    coordinates = coordinates,
    plot = plot,
    represented_percent = round(represented, 1),
    analysis_columns = names(analysis_data)
  )
}

foo <- function(dat) {
  unique_counts <- vapply(dat, function(values) length(unique(values)), integer(1L))
  which(unique_counts <= 1L)
}

#' Find each home cell's nearest target cell
#'
#' @param cells Data frame to augment.
#' @param home,target Data frames identifying home and candidate target cells.
#' @param label Label appended to output column names.
#' @param id_column,x_column,y_column Identifier and coordinate columns.
#' @return `cells` with nearest-cell identifier and distance columns.
nearest_n <- function(cells, home, target, label = "CellToCell",
                      id_column = "Number", x_column = "X_ANC_dna",
                      y_column = "Y_ANC_dna") {
  required <- c(id_column, x_column, y_column)
  .assert_data_columns(cells, required)
  .assert_data_columns(home, id_column)
  .assert_data_columns(target, required)
  .assert_numeric_columns(cells, c(x_column, y_column))
  .assert_numeric_columns(target, c(x_column, y_column))

  result <- cells
  nearest_column <- paste0("nearest_", label)
  distance_column <- paste0(nearest_column, "_distance")
  result[[nearest_column]] <- NA
  result[[distance_column]] <- NA_real_
  home_ids <- unique(home[[id_column]])

  for (index in seq_len(nrow(result))) {
    if (!result[[id_column]][index] %in% home_ids) {
      next
    }
    candidates <- target[target[[id_column]] != result[[id_column]][index], , drop = FALSE]
    if (nrow(candidates) == 0L) {
      next
    }
    distances <- sqrt(
      (candidates[[x_column]] - result[[x_column]][index])^2 +
        (candidates[[y_column]] - result[[y_column]][index])^2
    )
    nearest_index <- which.min(distances)
    result[[nearest_column]][index] <- candidates[[id_column]][nearest_index]
    result[[distance_column]][index] <- distances[nearest_index]
  }
  result
}

#' Characterize spatial clusters bounded by edge cells
#'
#' @param df Input data frame.
#' @param clust_pop Population-classification column.
#' @param center Value identifying potential cluster centers.
#' @param edge Value that stops nearest-neighbor expansion.
#' @param cluster_cutoff Minimum Jaccard overlap used to join center neighborhoods.
#' @param id_column,x_column,y_column Identifier and coordinate columns.
#' @return `df` with `cluster_size`, `cluster_mates`, and `cluster_id` columns.
cluster_count <- function(df, clust_pop = "tag", center = "Positive",
                          edge = "Negative", cluster_cutoff = 0.5,
                          id_column = "Number", x_column = "X_ANC_dna",
                          y_column = "Y_ANC_dna") {
  required <- c(clust_pop, id_column, x_column, y_column)
  .assert_data_columns(df, required)
  .assert_numeric_columns(df, c(x_column, y_column))
  if (!is.numeric(cluster_cutoff) || length(cluster_cutoff) != 1L ||
        cluster_cutoff < 0 || cluster_cutoff > 1) {
    stop("`cluster_cutoff` must be between zero and one.", call. = FALSE)
  }

  result <- df
  result$cluster_size <- 0L
  result$cluster_mates <- NA_character_
  result$cluster_id <- NA_integer_
  center_rows <- which(result[[clust_pop]] == center)
  neighborhoods <- vector("list", length(center_rows))

  for (position in seq_along(center_rows)) {
    index <- center_rows[position]
    candidates <- result[-index, , drop = FALSE]
    distances <- sqrt(
      (candidates[[x_column]] - result[[x_column]][index])^2 +
        (candidates[[y_column]] - result[[y_column]][index])^2
    )
    candidates <- candidates[order(distances), , drop = FALSE]
    edge_position <- which(candidates[[clust_pop]] == edge)[1L]
    if (is.na(edge_position)) {
      neighborhood <- candidates[[id_column]]
    } else if (edge_position == 1L) {
      neighborhood <- candidates[[id_column]][0L]
    } else {
      neighborhood <- candidates[[id_column]][seq_len(edge_position - 1L)]
    }
    neighborhood <- sort(unique(c(result[[id_column]][index], neighborhood)))
    neighborhoods[[position]] <- neighborhood
    result$cluster_size[index] <- length(neighborhood)
    result$cluster_mates[index] <- paste(neighborhood, collapse = "_")
  }

  if (length(center_rows) == 0L) {
    return(result)
  }
  parent <- seq_along(center_rows)
  root <- function(index) {
    while (parent[index] != index) {
      index <- parent[index]
    }
    index
  }
  for (left in seq_along(neighborhoods)) {
    if (left == length(neighborhoods)) {
      next
    }
    for (right in seq.int(left + 1L, length(neighborhoods))) {
      union_size <- length(union(neighborhoods[[left]], neighborhoods[[right]]))
      overlap <- length(intersect(neighborhoods[[left]], neighborhoods[[right]]))
      score <- if (union_size == 0L) 0 else overlap / union_size
      if (score > cluster_cutoff) {
        left_root <- root(left)
        right_root <- root(right)
        parent[right_root] <- left_root
      }
    }
  }
  roots <- vapply(seq_along(parent), root, integer(1L))
  result$cluster_id[center_rows] <- match(roots, unique(roots))
  result
}

color_blind <- c(
  "#999999", "#E69F00", "#56B4E9", "#009E73",
  "#F0E442", "#0072B2", "#D55E00", "#CC79A7"
)

cb_black <- c(
  "#000000", "#E69F00", "#56B4E9", "#009E73",
  "#F0E442", "#0072B2", "#D55E00", "#CC79A7"
)

cell_cycle_colors <- c(
  "#D4D4D4", "#98C84C", "#23B8CC", "#F16B1A", "#E5001C"
)
