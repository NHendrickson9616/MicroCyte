.require_icellate_pkgs <- function(caller) {
  if (!requireNamespace("imager", quietly = TRUE)) {
    stop(
      "`", caller, "()` requires the optional package 'imager'. ",
      "Install it with `install.packages(\"imager\")`.",
      call. = FALSE
    )
  }
  if (!requireNamespace("magick", quietly = TRUE)) {
    stop(
      "`", caller, "()` requires the optional package 'magick'. ",
      "Install it with `install.packages(\"magick\")`.",
      call. = FALSE
    )
  }
}

.microcyte_cell_geometry <- function(target_cell, d_multiplier) {
  area_name <- grep("^Area_ANC", names(target_cell), value = TRUE)[1L]
  x_name <- grep("^X_ANC", names(target_cell), value = TRUE)[1L]
  y_name <- grep("^Y_ANC", names(target_cell), value = TRUE)[1L]
  if (anyNA(c(area_name, x_name, y_name))) {
    stop(
      "`target_cells` must include columns beginning with Area_ANC, X_ANC, and Y_ANC.",
      call. = FALSE
    )
  }

  area <- as.numeric(target_cell[[area_name]][1L])
  x_pos <- as.integer(round(as.numeric(target_cell[[x_name]][1L])))
  y_pos <- as.integer(round(as.numeric(target_cell[[y_name]][1L])))
  if (anyNA(c(area, x_pos, y_pos)) || area < 0) {
    stop("Cell area and centroid coordinates must be finite numeric values.", call. = FALSE)
  }
  cell_length <- as.integer(round(d_multiplier * round(sqrt(area / pi))))
  list(
    cell_length = max(1L, cell_length),
    x_pos = x_pos,
    y_pos = y_pos
  )
}

.microcyte_crop_pixels <- function(pixels, geometry) {
  keep <- pixels$x >= geometry$x_pos - geometry$cell_length &
    pixels$x <= geometry$x_pos + geometry$cell_length &
    pixels$y >= geometry$y_pos - geometry$cell_length &
    pixels$y <= geometry$y_pos + geometry$cell_length
  cropped <- pixels[keep, , drop = FALSE]
  if (!nrow(cropped)) {
    stop("The requested cell crop does not overlap the source image.", call. = FALSE)
  }
  original_min_x <- min(cropped$x)
  original_min_y <- min(cropped$y)
  new_x <- 1 + geometry$x_pos - original_min_x
  new_y <- 1 + geometry$y_pos - original_min_y
  cropped$x <- 1 + cropped$x - original_min_x
  cropped$y <- 1 + cropped$y - original_min_y
  list(pixels = cropped, mid_x = new_x, mid_y = new_y)
}

.microcyte_add_fusion <- function(accumulator, pixels, mid_x, mid_y) {
  if (is.null(accumulator)) {
    pixels$number <- 1L
    return(list(data = pixels, mid_x = mid_x, mid_y = mid_y))
  }

  incoming <- pixels
  incoming$x <- incoming$x + accumulator$mid_x - mid_x
  incoming$y <- incoming$y + accumulator$mid_y - mid_y
  incoming <- incoming[incoming$x > 0 & incoming$y > 0, , drop = FALSE]
  old <- accumulator$data
  count <- max(old$number) + 1L
  keys <- intersect(c("x", "y", "cc"), intersect(names(old), names(incoming)))
  merged <- merge(
    old[, c(keys, "value"), drop = FALSE],
    incoming[, c(keys, "value"), drop = FALSE],
    by = keys,
    all = TRUE,
    suffixes = c(".old", ".new")
  )
  merged$value.old[is.na(merged$value.old)] <- 0
  merged$value.new[is.na(merged$value.new)] <- 0
  merged$value <- merged$value.old + merged$value.new
  merged$number <- count
  merged <- merged[, c(keys, "value", "number"), drop = FALSE]
  list(data = merged, mid_x = accumulator$mid_x, mid_y = accumulator$mid_y)
}

icellate <- function(target_cells,
                     folder_name = "singleCells",
                     d_multiplier = 0.3,
                     verify_size = FALSE,
                     fuse = FALSE,
                     marked = FALSE,
                     heat_intensity = FALSE,
                     heat_color = "inferno",
                     verify_image = "overlay.png",
                     randomize = TRUE,
                     sampling_number = 5,
                     line_analyses = TRUE,
                     number_of_lines = 2,
                     line_graphs = TRUE,
                     initial_angle = 0,
                     random_sampler_denom = 10,
                     root = ".",
                     input_dir = file.path(root, "files"),
                     output_dir = file.path(root, "icellates"),
                     heat_map_fun = NULL) {
  .require_icellate_pkgs("icellate")
  if (!is.data.frame(target_cells) || !nrow(target_cells)) {
    stop("`target_cells` must be a non-empty data frame.", call. = FALSE)
  }
  required <- c("name_id", "image", "Number")
  missing_columns <- setdiff(required, names(target_cells))
  if (length(missing_columns)) {
    stop(
      "`target_cells` is missing required columns: ",
      paste(missing_columns, collapse = ", "),
      call. = FALSE
    )
  }
  if (!dir.exists(input_dir)) {
    stop("Input directory not found: ", input_dir, call. = FALSE)
  }
  if (heat_intensity && !is.function(heat_map_fun)) {
    stop("Supply `heat_map_fun` when `heat_intensity = TRUE`.", call. = FALSE)
  }
  if (!is.numeric(sampling_number) || length(sampling_number) != 1L ||
        is.na(sampling_number) || sampling_number < 1) {
    stop("`sampling_number` must be a positive number.", call. = FALSE)
  }

  folder_dir <- file.path(output_dir, folder_name)
  dir.create(folder_dir, recursive = TRUE, showWarnings = FALSE)
  if (randomize) {
    sample_size <- min(nrow(target_cells), as.integer(sampling_number))
    target_cells <- target_cells[sample.int(nrow(target_cells), sample_size), , drop = FALSE]
  }

  size_verified <- FALSE
  fusion <- list()
  cell_results <- list()

  for (cell_index in seq_len(nrow(target_cells))) {
    target_cell <- target_cells[cell_index, , drop = FALSE]
    target_name_id <- as.character(target_cell$name_id[1L])
    target_image <- as.character(target_cell$image[1L])
    cell_id <- paste0(target_name_id, "-", target_cell$Number[1L])
    cell_dir <- file.path(folder_dir, cell_id)
    source_dir <- file.path(input_dir, target_name_id, target_image, "PNGS")
    if (!dir.exists(source_dir)) {
      warning("Image directory not found; skipping cell: ", source_dir, call. = FALSE)
      next
    }

    dir.create(cell_dir, recursive = TRUE, showWarnings = FALSE)
    source_images <- list.files(
      source_dir,
      pattern = "\\.png$",
      full.names = TRUE,
      ignore.case = TRUE
    )
    if (!length(source_images)) {
      warning("No PNG images found; skipping cell: ", source_dir, call. = FALSE)
      next
    }
    copied <- file.copy(source_images, cell_dir, overwrite = TRUE)
    if (!all(copied)) {
      stop("Failed to copy one or more images into: ", cell_dir, call. = FALSE)
    }
    cell_images <- file.path(cell_dir, basename(source_images))
    geometry <- .microcyte_cell_geometry(target_cell, d_multiplier)

    if (verify_size && !size_verified) {
      verify_path <- file.path(cell_dir, verify_image)
      if (!file.exists(verify_path)) {
        message(
          "The size-verification image was not found at ", verify_path,
          "; using the supplied `d_multiplier`."
        )
      } else {
        repeat {
          geometry <- .microcyte_cell_geometry(target_cell, d_multiplier)
          preview <- .microcyte_crop_pixels(
            as.data.frame(imager::load.image(verify_path)),
            geometry
          )
          graphics::plot(suppressWarnings(imager::as.cimg(preview$pixels)))
          message("d_multiplier is currently set to: ", d_multiplier)
          answer <- readline("Does this look like the appropriate size (Y/n)? ")
          if (tolower(answer) != "n") {
            size_verified <- TRUE
            break
          }
          proposed <- suppressWarnings(as.numeric(readline("New d_multiplier: ")))
          if (!is.finite(proposed) || proposed <= 0) {
            message("Please enter a positive numeric multiplier.")
          } else {
            d_multiplier <- proposed
          }
        }
      }
    }

    generated <- character()
    last_crop <- NULL
    for (image_path in cell_images) {
      loaded <- imager::load.image(image_path)
      crop <- .microcyte_crop_pixels(as.data.frame(loaded), geometry)
      cropped_image <- suppressWarnings(imager::as.cimg(crop$pixels))
      cropped_magick <- imager::cimg2magick(cropped_image, rotate = TRUE)
      magick::image_write(cropped_magick, path = image_path, format = "png")
      generated <- c(generated, image_path)
      last_crop <- crop

      if (marked) {
        marked_pixels <- crop$pixels
        crosshair <- max(1, geometry$cell_length / 20)
        on_horizontal <- marked_pixels$x > crop$mid_x - crosshair &
          marked_pixels$x < crop$mid_x + crosshair & marked_pixels$y == crop$mid_y
        on_vertical <- marked_pixels$y > crop$mid_y - crosshair &
          marked_pixels$y < crop$mid_y + crosshair & marked_pixels$x == crop$mid_x
        marked_pixels$value[on_horizontal | on_vertical] <- 1
        marked_image <- suppressWarnings(imager::as.cimg(marked_pixels))
        marked_path <- file.path(
          cell_dir,
          paste0(tools::file_path_sans_ext(basename(image_path)), "_marked.png")
        )
        magick::image_write(
          imager::cimg2magick(marked_image, rotate = TRUE),
          path = marked_path,
          format = "png"
        )
        generated <- c(generated, marked_path)
      }

      if (heat_intensity) {
        heat_result <- heat_map_fun(
          image_location = image_path,
          heat_color = heat_color,
          intensity_binning = FALSE,
          output_dir = cell_dir
        )
        heat_path <- file.path(cell_dir, paste0("heatMap_", basename(image_path)))
        if (!file.rename(heat_result$output_path, heat_path)) {
          stop("Failed to create heat-map output: ", heat_path, call. = FALSE)
        }
        generated <- c(generated, heat_path)
      }

      if (fuse && max(crop$pixels$x) == max(crop$pixels$y)) {
        image_key <- tools::file_path_sans_ext(basename(image_path))
        fusion[[image_key]] <- .microcyte_add_fusion(
          fusion[[image_key]],
          crop$pixels,
          crop$mid_x,
          crop$mid_y
        )
      }
    }

    line_result <- NULL
    if (line_analyses && !is.null(last_crop)) {
      line_result <- line_analysis(
        mid_x = last_crop$mid_x,
        mid_y = last_crop$mid_y,
        image_name = last_crop$pixels,
        lines_number = number_of_lines,
        line_length = geometry$cell_length,
        total_corr_pix_fraction = random_sampler_denom,
        graphing = line_graphs,
        offset = initial_angle,
        input_dir = cell_dir,
        output_dir = file.path(cell_dir, "line_analyses")
      )
    }

    cell_results[[cell_id]] <- list(
      cell_dir = cell_dir,
      source_dir = source_dir,
      output_paths = generated,
      geometry = geometry,
      line_analysis = line_result
    )
  }

  fusion_paths <- character()
  fusion_metrics_path <- NULL
  fusion_metrics <- data.frame(image = character(), number = integer())
  if (fuse && length(fusion)) {
    metric_rows <- vector("list", length(fusion))
    fusion_names <- names(fusion)
    for (fusion_index in seq_along(fusion_names)) {
      image_key <- fusion_names[[fusion_index]]
      fused <- fusion[[image_key]]$data
      count <- max(fused$number)
      fused$value <- fused$value / count
      fused <- fused[, setdiff(names(fused), "number"), drop = FALSE]
      fused_path <- file.path(folder_dir, paste0(image_key, ".png"))
      fused_image <- suppressWarnings(imager::as.cimg(fused))
      magick::image_write(
        imager::cimg2magick(fused_image, rotate = TRUE),
        path = fused_path,
        format = "png"
      )
      fusion_paths <- c(fusion_paths, fused_path)
      metric_rows[[fusion_index]] <- data.frame(
        image = image_key,
        number = count,
        stringsAsFactors = FALSE
      )
    }
    fusion_metrics <- do.call(rbind, metric_rows)
    row.names(fusion_metrics) <- NULL
    fusion_metrics_path <- file.path(folder_dir, "ClassAveragingMetrics.csv")
    utils::write.csv(fusion_metrics, fusion_metrics_path, row.names = FALSE)
  }

  invisible(list(
    output_dir = folder_dir,
    cells = cell_results,
    fusion = fusion,
    fusion_paths = fusion_paths,
    fusion_metrics = fusion_metrics,
    fusion_metrics_path = fusion_metrics_path,
    d_multiplier = d_multiplier
  ))
}

line_analysis <- function(mid_x,
                          mid_y,
                          image_name,
                          lines_number = 2,
                          line_length,
                          total_corr_pix_fraction = 10,
                          graphing = TRUE,
                          offset = 0,
                          input_dir = ".",
                          output_dir = file.path(input_dir, "line_analyses")) {
  if (!requireNamespace("imager", quietly = TRUE)) {
    stop(
      "`line_analysis()` requires the optional package 'imager'. ",
      "Install it with `install.packages(\"imager\")`.",
      call. = FALSE
    )
  }
  if (!dir.exists(input_dir)) {
    stop("Input directory not found: ", input_dir, call. = FALSE)
  }
  if (!is.numeric(total_corr_pix_fraction) || total_corr_pix_fraction <= 0) {
    stop("`total_corr_pix_fraction` must be positive.", call. = FALSE)
  }
  if (!is.numeric(lines_number) || lines_number < 1) {
    stop("`lines_number` must be at least 1.", call. = FALSE)
  }

  csv_dir <- file.path(output_dir, "csv")
  dir.create(csv_dir, recursive = TRUE, showWarnings = FALSE)
  png_paths <- list.files(
    input_dir,
    pattern = "\\.png$",
    full.names = TRUE,
    ignore.case = TRUE
  )
  excluded <- grepl("_marked\\.png$", basename(png_paths), ignore.case = TRUE) |
    basename(png_paths) == "overlay.png" |
    grepl("^heatMap_", basename(png_paths))
  png_paths <- png_paths[!excluded]
  if (!length(png_paths)) {
    stop("No analyzable PNG images found in: ", input_dir, call. = FALSE)
  }

  channel_frames <- lapply(png_paths, function(image_path) {
    pixels <- as.data.frame(imager::grayscale(imager::load.image(image_path)))
    pixels$value <- pixels$value - min(pixels$value, na.rm = TRUE)
    channel <- tools::file_path_sans_ext(basename(image_path))
    names(pixels)[names(pixels) == "value"] <- channel
    pixels[, c("x", "y", channel), drop = FALSE]
  })
  line_data <- Reduce(
    function(left, right) merge(left, right, by = c("x", "y"), all = FALSE),
    channel_frames
  )
  channel_names <- setdiff(names(line_data), c("x", "y"))

  random_count <- min(
    nrow(line_data),
    max(0L, as.integer(floor(nrow(line_data) / total_corr_pix_fraction)))
  )
  random_points <- if (random_count) {
    line_data[sample.int(nrow(line_data), random_count), , drop = FALSE]
  } else {
    line_data[FALSE, , drop = FALSE]
  }
  random_path <- file.path(output_dir, "WholeCellRandom.csv")
  utils::write.csv(random_points, random_path, row.names = FALSE)

  line_paths <- character()
  if (lines_number < 3) {
    x_line <- line_data[
      line_data$y == mid_y &
        line_data$x <= mid_x + line_length & line_data$x >= mid_x - line_length,
      ,
      drop = FALSE
    ]
    y_line <- line_data[
      line_data$x == mid_x &
        line_data$y <= mid_y + line_length & line_data$y >= mid_y - line_length,
      ,
      drop = FALSE
    ]
    line_paths <- file.path(csv_dir, c("xLine.csv", "yLine.csv"))
    utils::write.csv(x_line, line_paths[1L], row.names = FALSE)
    utils::write.csv(y_line, line_paths[2L], row.names = FALSE)
  } else {
    angles <- offset + seq_len(as.integer(lines_number)) * (360 / lines_number)
    for (line_angle in angles) {
      angle_radians <- line_angle * pi / 180
      x_start <- round(sin(angle_radians) * line_length + mid_x)
      y_start <- round(cos(angle_radians) * line_length + mid_y)
      x_end <- round(mid_x - sin(angle_radians) * line_length)
      y_end <- round(mid_y - cos(angle_radians) * line_length)
      if (mid_x != x_start && mid_y != y_start) {
        line_slope <- -(mid_y - y_start) / (mid_x - x_start)
        y_intercept <- mid_y + line_slope * mid_x
        keep <- line_data$x <= max(x_start, x_end) &
          line_data$x >= min(x_start, x_end) &
          line_data$y == -round(line_slope * line_data$x - y_intercept)
      } else {
        keep <- line_data$x <= max(x_start, x_end) &
          line_data$x >= min(x_start, x_end) &
          line_data$y <= max(y_start, y_end) &
          line_data$y >= min(y_start, y_end)
      }
      new_line <- line_data[keep, , drop = FALSE]
      line_path <- file.path(csv_dir, paste0("lineAngle_", line_angle, ".csv"))
      utils::write.csv(new_line, line_path, row.names = FALSE)
      line_paths <- c(line_paths, line_path)
    }
  }

  graph_paths <- character()
  plots <- list()
  long_paths <- character()
  if (graphing) {
    for (line_path in line_paths) {
      graph_data <- utils::read.csv(line_path, check.names = FALSE)
      graph_data$pix <- seq_len(nrow(graph_data))
      available_channels <- intersect(channel_names, names(graph_data))
      if (!length(available_channels)) {
        next
      }
      graph_data <- tidyr::pivot_longer(
        graph_data,
        cols = tidyselect::all_of(available_channels),
        names_to = "color",
        values_to = "intensity"
      )
      long_path <- file.path(csv_dir, paste0("longer_", basename(line_path)))
      utils::write.csv(graph_data, long_path, row.names = FALSE)
      long_paths <- c(long_paths, long_path)

      plot_object <- ggplot2::ggplot(
        graph_data,
        ggplot2::aes(
          x = rlang::.data$pix,
          y = rlang::.data$intensity,
          color = rlang::.data$color
        )
      ) +
        ggplot2::geom_step(linewidth = 1.5) +
        ggplot2::xlab("Relative pixel position") +
        ggplot2::ylab("Relative pixel intensity") +
        ggplot2::theme_classic() +
        ggplot2::theme(
          axis.line = ggplot2::element_line(colour = "black", linewidth = 1.5),
          axis.ticks = ggplot2::element_line(colour = "black", linewidth = 1.5),
          axis.text = ggplot2::element_text(face = "bold", color = "black", size = 12),
          axis.title = ggplot2::element_text(face = "bold", color = "black", size = 20),
          legend.position = "top"
        )
      graph_path <- file.path(
        output_dir,
        paste0(tools::file_path_sans_ext(basename(line_path)), ".pdf")
      )
      ggplot2::ggsave(
        graph_path,
        plot = plot_object,
        height = 5,
        width = 5,
        units = "in"
      )
      plots[[basename(line_path)]] <- plot_object
      graph_paths <- c(graph_paths, graph_path)
    }
  }

  invisible(list(
    data = line_data,
    output_dir = output_dir,
    random_sample_path = random_path,
    line_paths = line_paths,
    long_data_paths = long_paths,
    graph_paths = graph_paths,
    plots = plots,
    image_name = image_name
  ))
}
