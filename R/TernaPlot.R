# Generate and customize ternary plots without modifying package-global state.

ternaplot <- function(data,
                      customize_it = TRUE,
                      save_it = TRUE,
                      output_path = NULL) {
  if (is.data.frame(data)) {
    source_data <- data
  } else if (is.character(data) && length(data) == 1L) {
    source_data <- utils::read.csv(file = data, header = TRUE)
  } else {
    stop("`data` must be a data frame or a single CSV file path.", call. = FALSE)
  }

  print("Data categories detected:")
  print(names(source_data))
  samples <- readline(prompt = "Which category is the sample ID: ")
  cat("\n")
  lname <- readline(prompt = "Which category is left point: ")
  llabel <- readline(prompt = "What should the left label be: ")
  cat("\n")
  tname <- readline(prompt = "Which category is top point: ")
  tlabel <- readline(prompt = "What should the top label be: ")
  cat("\n")
  rname <- readline(prompt = "Which category is right point: ")
  rlabel <- readline(prompt = "What should the right label be: ")
  cat("\n")

  requested_columns <- c(samples, lname, tname, rname)
  missing_columns <- setdiff(requested_columns, names(source_data))
  if (length(missing_columns) > 0L) {
    stop(
      "Unknown column(s): ",
      paste(missing_columns, collapse = ", "),
      call. = FALSE
    )
  }

  ternary_data <- data.frame(
    sampleName = source_data[[samples]],
    lpoint = source_data[[lname]],
    tpoint = source_data[[tname]],
    rpoint = source_data[[rname]],
    stringsAsFactors = FALSE
  )
  print("Dataframe generated.")
  cat("\n")

  labels <- list(left = llabel, top = tlabel, right = rlabel)
  plot <- ggtern::ggtern(
    data = ternary_data,
    ggplot2::aes(
      x = rlang::.data$lpoint,
      y = rlang::.data$tpoint,
      z = rlang::.data$rpoint
    )
  ) +
    ggplot2::geom_point(size = 3) +
    ggplot2::theme_bw() +
    ggtern::Tlab(labels$top) +
    ggtern::Rlab(labels$right) +
    ggtern::Llab(labels$left)

  print(plot)
  customization <- NULL
  if (customize_it == TRUE) {
    customization <- customize(
      plot = plot,
      data = source_data,
      ternary_data = ternary_data
    )
    ternary_data <- customization$data
    plot <- customization$plot
  } else {
    cat("\n")
    print("Call customize() with a plot and its data to make things fancy.")
  }

  saved <- NULL
  if (save_it == TRUE) {
    if (is.null(output_path)) {
      stop("`output_path` must be supplied when `save_it` is TRUE.", call. = FALSE)
    }
    saved <- save_png(plot = plot, path = output_path)
  } else {
    print("Call save_png() with a plot and path to save a high-resolution PNG.")
    cat("\n")
  }

  structure(
    list(
      source_data = source_data,
      ternary_data = ternary_data,
      plot = plot,
      customization = customization,
      saved = saved
    ),
    class = c("microcyte_ternaplot", "list")
  )
}

save_png <- function(plot, path) {
  if (!inherits(plot, "ggplot")) {
    stop("`plot` must be a ggplot object.", call. = FALSE)
  }
  if (!is.character(path) || length(path) != 1L || !nzchar(path)) {
    stop("`path` must be a single non-empty file path.", call. = FALSE)
  }

  if (!grepl("\\.png$", path, ignore.case = TRUE)) {
    path <- paste0(path, ".png")
  }
  print(plot)
  ggplot2::ggsave(
    filename = path,
    plot = plot,
    dpi = 300,
    units = "in",
    width = 11,
    height = 8.5
  )
  print(paste("PNG saved to", path))

  list(plot = plot, path = path)
}

customize <- function(plot, data, ternary_data) {
  if (!inherits(plot, "ggplot")) {
    stop("`plot` must be a ggplot object.", call. = FALSE)
  }
  if (!is.data.frame(data) || !is.data.frame(ternary_data)) {
    stop("`data` and `ternary_data` must be data frames.", call. = FALSE)
  }

  print("Data categories detected:")
  print(names(data))

  colz <- readline(prompt = "Would you like to color the dots based on a variable (y/n): ")
  coln <- NULL
  if (colz == "y") {
    colv <- readline(prompt = "Which category should color be based on: ")
    if (!colv %in% names(data)) {
      stop("Unknown color column: ", colv, call. = FALSE)
    }
    coln <- readline(prompt = "What is the name of this variable: ")
    ternary_data[["color_value"]] <- data[[colv]]
  }
  cat("\n")

  shapez <- readline(
    prompt = "Would you like to alter shape based on a categorical variable (y/n): "
  )
  shapen <- NULL
  if (shapez == "y") {
    shapev <- readline(prompt = "Which category should shape be based on: ")
    if (!shapev %in% names(data)) {
      stop("Unknown shape column: ", shapev, call. = FALSE)
    }
    shapen <- readline(prompt = "What is the name of this variable: ")
    ternary_data[["shape_value"]] <- data[[shapev]]
  }
  cat("\n")

  conz <- readline(prompt = "Would you like to add a confidence interval (y/n): ")
  confidence <- NULL
  if (conz == "y") {
    confidence <- as.integer(readline(
      prompt = "What percentage of confidence should be used: "
    )) / 100
  }

  plot$data <- ternary_data
  if (colz == "y" && shapez == "y") {
    custom_plot <- plot +
      ggplot2::aes(
        color = rlang::.data$color_value,
        shape = rlang::.data$shape_value
      ) +
      ggplot2::guides(
        shape = ggplot2::guide_legend(title = shapen),
        color = ggplot2::guide_legend(title = coln)
      )
  } else if (colz == "y" && shapez != "y") {
    custom_plot <- plot +
      ggplot2::aes(color = rlang::.data$color_value) +
      ggplot2::guides(color = ggplot2::guide_legend(title = coln))
  } else if (colz != "y" && shapez == "y") {
    custom_plot <- plot +
      ggplot2::aes(shape = rlang::.data$shape_value) +
      ggplot2::guides(shape = ggplot2::guide_legend(title = shapen))
  } else {
    custom_plot <- plot
  }

  if (conz == "y") {
    custom_plot <- custom_plot + ggtern::geom_confidence_tern(breaks = confidence)
  }

  cell_cycle_colors <- c(
    "#D4D4D4", "#98C84C", "#23B8CC", "#F16B1A", "#E5001C"
  )
  custom_plot <- custom_plot + ggtern::theme_custom(
    base_size = 20,
    col.T = cell_cycle_colors[[2L]],
    col.L = cell_cycle_colors[[1L]],
    col.R = cell_cycle_colors[[3L]],
    tern.panel.background = "white",
    tern.plot.background = NULL
  )

  print(custom_plot)
  cat("\n")

  list(
    plot = custom_plot,
    data = ternary_data,
    settings = list(
      color = coln,
      shape = shapen,
      confidence = confidence
    )
  )
}
