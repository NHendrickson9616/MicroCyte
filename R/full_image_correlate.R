full_image_correlate <- function(cor_images = c("dna", "edu"),
                                 sub_percent = 50,
                                 background_subtract = TRUE,
                                 threshold_intensity = 0.1,
                                 remove_overexposed = TRUE,
                                 graph_dist = TRUE,
                                 root = ".",
                                 input_dir = file.path(root, "files"),
                                 schema_path = file.path(root, "schema.csv"),
                                 output_dir = file.path(root, "data"),
                                 figure_dir = file.path(root, "figures")) {
  if (!requireNamespace("imager", quietly = TRUE)) {
    stop(
      "`full_image_correlate()` requires the optional package 'imager'. ",
      "Install it with `install.packages(\"imager\")`.",
      call. = FALSE
    )
  }
  if (!requireNamespace("magick", quietly = TRUE)) {
    stop(
      "`full_image_correlate()` requires the optional package 'magick'. ",
      "Install it with `install.packages(\"magick\")`.",
      call. = FALSE
    )
  }
  if (!file.exists(schema_path)) {
    stop("Schema file not found: ", schema_path, call. = FALSE)
  }
  if (!dir.exists(input_dir)) {
    stop("Input directory not found: ", input_dir, call. = FALSE)
  }
  if (!length(cor_images) || anyNA(cor_images) || any(!nzchar(cor_images))) {
    stop("`cor_images` must contain at least one non-empty image name.", call. = FALSE)
  }
  if (!is.numeric(sub_percent) || length(sub_percent) != 1L ||
        is.na(sub_percent) || sub_percent < 0 || sub_percent > 100) {
    stop("`sub_percent` must be a number between 0 and 100.", call. = FALSE)
  }

  scheme <- readr::read_csv(schema_path, show_col_types = FALSE)
  sample_dirs <- list.dirs(input_dir, recursive = FALSE, full.names = TRUE)
  experiment_parts <- list()
  image_outputs <- character()
  sample_outputs <- character()
  skipped_images <- character()

  for (sample_dir in sample_dirs) {
    sample_id <- basename(sample_dir)
    message("Working in file: ", sample_id)
    image_dirs <- list.dirs(sample_dir, recursive = FALSE, full.names = TRUE)
    sample_parts <- list()

    for (image_dir in image_dirs) {
      image_id <- basename(image_dir)
      png_dir <- file.path(image_dir, "PNGS")
      if (!dir.exists(png_dir)) {
        next
      }
      message("Working in image: ", image_id)
      channel_data <- list()

      for (channel in cor_images) {
        image_path <- file.path(png_dir, paste0(channel, ".png"))
        if (!file.exists(image_path)) {
          skipped_images <- c(skipped_images, image_path)
          message("Unable to find ", image_path, "; skipping this image.")
          next
        }

        image <- imager::load.image(image_path)
        pixels <- as.data.frame(image)
        if ("cc" %in% names(pixels) && 4L %in% unique(pixels$cc)) {
          message("Alpha channel detected in ", image_path, "; removing it.")
          pixels <- pixels[pixels$cc != 4L, , drop = FALSE]
          image <- suppressWarnings(imager::as.cimg(pixels))
          converted <- imager::cimg2magick(image, rotate = TRUE)
          converted <- magick::image_flop(converted)
          magick::image_write(converted, path = image_path, format = "png")
          image <- imager::load.image(image_path)
        }

        pixels <- as.data.frame(imager::grayscale(image))
        if (background_subtract && nrow(pixels)) {
          pixels$value <- pixels$value - min(pixels$value, na.rm = TRUE)
        }
        names(pixels)[names(pixels) == "value"] <- channel
        channel_data[[channel]] <- pixels[, c("x", "y", channel), drop = FALSE]
      }

      if (!all(cor_images %in% names(channel_data))) {
        next
      }

      image_df <- Reduce(
        function(left, right) merge(left, right, by = c("x", "y"), all = FALSE),
        channel_data[cor_images]
      )
      image_df$image <- image_id
      for (channel in cor_images) {
        keep <- !is.na(image_df[[channel]]) & image_df[[channel]] > threshold_intensity
        if (remove_overexposed) {
          keep <- keep & image_df[[channel]] < 1
        }
        image_df <- image_df[keep, , drop = FALSE]
      }

      row_keep <- min(nrow(image_df), round((sub_percent / 100) * nrow(image_df)))
      if (row_keep > 0L) {
        image_df <- image_df[sample.int(nrow(image_df), row_keep), , drop = FALSE]
      } else {
        image_df <- image_df[FALSE, , drop = FALSE]
      }

      image_output <- file.path(png_dir, "image_corr.csv")
      readr::write_csv(image_df, image_output)
      image_outputs <- c(image_outputs, image_output)
      sample_parts[[image_id]] <- image_df
    }

    if (!length(sample_parts)) {
      next
    }

    sample_data <- do.call(rbind, unname(sample_parts))
    row.names(sample_data) <- NULL
    metadata <- scheme[scheme$name_id == sample_id, , drop = FALSE]
    if (nrow(metadata) && ncol(metadata) >= 7L) {
      for (metadata_name in names(metadata)[seq.int(7L, ncol(metadata))]) {
        sample_data[[metadata_name]] <- as.character(metadata[[metadata_name]][1L])
      }
    }

    sample_output <- file.path(sample_dir, "all_image_corr.csv")
    readr::write_csv(sample_data, sample_output)
    sample_outputs <- c(sample_outputs, sample_output)
    experiment_parts[[sample_id]] <- sample_data
  }

  experiment_data <- if (length(experiment_parts)) {
    result <- do.call(rbind, unname(experiment_parts))
    row.names(result) <- NULL
    result
  } else {
    data.frame()
  }

  experiment_output <- NULL
  figure_outputs <- character()
  plots <- list()
  if (nrow(experiment_data)) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    experiment_output <- file.path(output_dir, "experiment_correlate.csv")
    readr::write_csv(experiment_data, experiment_output)

    if (graph_dist) {
      dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
      for (channel in cor_images) {
        draft <- ggplot2::ggplot(
          experiment_data,
          ggplot2::aes(x = rlang::.data[[channel]])
        ) +
          ggplot2::geom_density() +
          ggplot2::xlab(paste0("Intensity profile of ", channel)) +
          ggplot2::theme_bw()
        figure_path <- file.path(figure_dir, paste0("Corr_plot_", channel, ".png"))
        ggplot2::ggsave(figure_path, plot = draft, dpi = 300)
        plots[[channel]] <- draft
        figure_outputs <- c(figure_outputs, figure_path)
      }
    }
  }

  invisible(list(
    data = experiment_data,
    output_path = experiment_output,
    image_output_paths = image_outputs,
    sample_output_paths = sample_outputs,
    figure_paths = figure_outputs,
    plots = plots,
    skipped_images = unique(skipped_images)
  ))
}
