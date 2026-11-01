.require_image_pkgs <- function(caller) {
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

.microcyte_rgb_image <- function(pixels, colors) {
  rgb <- grDevices::col2rgb(colors) / 255
  color_pixels <- pixels
  if (!"cc" %in% names(color_pixels)) {
    color_pixels$cc <- 1L
    red <- color_pixels
    green <- color_pixels
    blue <- color_pixels
    green$cc <- 2L
    blue$cc <- 3L
    color_pixels <- rbind(red, green, blue)
  }
  for (channel in seq_len(3L)) {
    channel_rows <- color_pixels$cc == channel
    color_pixels$value[channel_rows] <- rgb[channel, ]
  }
  suppressWarnings(imager::as.cimg(color_pixels))
}

heat_map_image <- function(image_location,
                           heat_color = "inferno",
                           intensity_binning = FALSE,
                           output_dir = dirname(image_location)) {
  .require_image_pkgs("heat_map_image")
  if (!file.exists(image_location)) {
    stop("Image file not found: ", image_location, call. = FALSE)
  }
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  extension <- tolower(tools::file_ext(image_location))
  converted_input <- NULL
  working_image <- image_location
  if (extension %in% c("tif", "tiff")) {
    message("TIF detected. Creating PNG.")
    converted_input <- file.path(
      output_dir,
      paste0(tools::file_path_sans_ext(basename(image_location)), ".png")
    )
    magick::image_write(
      magick::image_read(image_location),
      path = converted_input,
      format = "png"
    )
    working_image <- converted_input
  } else if (extension != "png") {
    stop("`image_location` must point to a PNG, TIF, or TIFF image.", call. = FALSE)
  }

  old_image <- imager::load.image(working_image)
  pixels <- as.data.frame(old_image)
  if ("cc" %in% names(pixels) && 4L %in% unique(pixels$cc)) {
    message("Alpha channel detected; removing it from the heat-mapped output.")
    pixels <- pixels[pixels$cc != 4L, , drop = FALSE]
    old_image <- suppressWarnings(imager::as.cimg(pixels))
  }
  gray <- as.data.frame(imager::grayscale(old_image))

  palette_plot <- ggplot2::ggplot(
    gray,
    ggplot2::aes(
      x = rlang::.data$x,
      y = rlang::.data$y,
      fill = rlang::.data$value
    )
  ) + ggplot2::geom_point()
  if (intensity_binning) {
    palette_plot <- palette_plot + ggplot2::scale_fill_viridis_b(option = heat_color)
  } else {
    palette_plot <- palette_plot + ggplot2::scale_fill_viridis_c(option = heat_color)
  }
  built <- ggplot2::ggplot_build(palette_plot)
  gray$hex <- built$data[[1L]]$fill

  bins <- NULL
  bins_path <- NULL
  base_name <- tools::file_path_sans_ext(basename(working_image))
  if (intensity_binning) {
    split_gray <- split(gray, gray$hex, drop = TRUE)
    bins <- do.call(
      rbind,
      lapply(split_gray, function(bin) {
        data.frame(
          Image = working_image,
          HexCode = bin$hex[1L],
          Mean_X = stats::median(bin$x),
          Mean_Y = stats::median(bin$y),
          Mean = mean(bin$value),
          Median = stats::median(bin$value),
          Range = diff(range(bin$value)),
          Total = sum(bin$value),
          Area_px = nrow(bin),
          stringsAsFactors = FALSE
        )
      })
    )
    row.names(bins) <- NULL
    bins_path <- file.path(output_dir, paste0(base_name, "_heatBins.csv"))
    readr::write_csv(bins, bins_path)
  }

  mapped_image <- .microcyte_rgb_image(pixels, gray$hex)
  mapped_image <- imager::cimg2magick(mapped_image, rotate = TRUE)
  mapped_image <- magick::image_flop(mapped_image)
  suffix <- if (intensity_binning) {
    "_binnedHeatMapped.png"
  } else {
    paste0("_", heat_color, "_HeatMapped.png")
  }
  output_path <- file.path(output_dir, paste0(base_name, suffix))
  magick::image_write(mapped_image, path = output_path, format = "png")
  message("Complete: ", output_path)

  invisible(list(
    output_path = output_path,
    converted_input_path = converted_input,
    bins_path = bins_path,
    bins = bins,
    pixel_data = gray
  ))
}

heatmap_it_all <- function(heatmap_color = "inferno",
                           binned_colors = FALSE,
                           root = ".",
                           input_dir = file.path(root, "files"),
                           output_dir = NULL) {
  .require_image_pkgs("heatmap_it_all")
  if (!dir.exists(input_dir)) {
    stop("Input directory not found: ", input_dir, call. = FALSE)
  }

  sample_dirs <- list.dirs(input_dir, recursive = FALSE, full.names = TRUE)
  results <- list()
  for (sample_dir in sample_dirs) {
    message("Working in file: ", basename(sample_dir))
    image_dirs <- list.dirs(sample_dir, recursive = FALSE, full.names = TRUE)
    for (image_dir in image_dirs) {
      png_dir <- file.path(image_dir, "PNGS")
      if (!dir.exists(png_dir) || grepl("ratio", basename(image_dir), fixed = TRUE) ||
            grepl(heatmap_color, basename(image_dir), fixed = TRUE)) {
        next
      }
      message("Working in image: ", basename(image_dir))
      images <- list.files(
        png_dir,
        pattern = "\\.png$",
        full.names = TRUE,
        ignore.case = TRUE
      )
      images <- images[!grepl("HeatMapped", basename(images), fixed = TRUE)]
      for (image_path in images) {
        destination <- if (is.null(output_dir)) {
          png_dir
        } else {
          file.path(output_dir, basename(sample_dir), basename(image_dir), "PNGS")
        }
        result <- suppressWarnings(heat_map_image(
          image_location = image_path,
          heat_color = heatmap_color,
          intensity_binning = binned_colors,
          output_dir = destination
        ))
        results[[image_path]] <- result
      }
    }
  }

  invisible(list(
    results = results,
    output_paths = vapply(results, `[[`, character(1L), "output_path"),
    input_dir = input_dir,
    output_dir = output_dir
  ))
}
