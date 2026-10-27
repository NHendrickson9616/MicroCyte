#!/usr/bin/Rscript

.purgo_directory <- function(path, argument = "directory") {
  if (!is.character(path) || length(path) != 1L || !dir.exists(path)) {
    stop(argument, " must identify an existing directory.")
  }
  normalizePath(path, mustWork = TRUE)
}

.purgo_file <- function(path, base_dir) {
  is_absolute <- grepl("^(/|[A-Za-z]:[/\\\\]|[/\\\\]{2})", path)
  if (is_absolute) path else file.path(base_dir, path)
}

# Convert all TIFF images in one image directory to PNG files.
abra <- function(image_dir = ".") {
  image_dir <- .purgo_directory(image_dir, "image_dir")
  png_dir <- file.path(image_dir, "PNGS")
  originals_dir <- file.path(image_dir, "originals")

  if (!dir.exists(png_dir) && !dir.create(png_dir)) {
    stop("Could not create directory: ", png_dir)
  }

  f_all <- list.files(path = image_dir, pattern = ".tif")
  image_number <- basename(image_dir)
  if (length(f_all) == 0L) {
    print(paste0("No images found in: ", image_number))
  } else if (length(f_all) > 4L) {
    print(paste0("Excess (>4) images found in: ", image_number))
  }

  if (!dir.exists(originals_dir)) {
    if (!dir.create(originals_dir)) {
      stop("Could not create directory: ", originals_dir)
    }
    if (length(f_all) > 0L) {
      file.copy(file.path(image_dir, f_all), originals_dir)
    }
  }

  output_files <- character()
  for (file_name in f_all) {
    image_name <- substr(file_name, 0L, nchar(file_name) - 4L)
    output_file <- file.path(png_dir, paste0(image_name, ".png"))
    image_item <- magick::image_read(file.path(image_dir, file_name))
    magick::image_write(image_item, path = output_file, format = "png")
    output_files <- c(output_files, output_file)
  }

  invisible(output_files)
}

# Extract grayscale pixel information from the PNG files in one image directory.
kadabra <- function(sample_number = 10000, image_dir = ".") {
  if (length(sample_number) != 1L || !is.numeric(sample_number) || is.na(sample_number) ||
        sample_number < 1 || sample_number != as.integer(sample_number)) {
    stop("sample_number must be a positive whole number.")
  }

  image_dir <- .purgo_directory(image_dir, "image_dir")
  png_dir <- file.path(image_dir, "PNGS")
  if (!dir.exists(png_dir)) {
    stop("No PNGS directory found in image_dir: ", image_dir)
  }

  qc_dir <- file.path(image_dir, "QC")
  if (!dir.exists(qc_dir) && !dir.create(qc_dir)) {
    stop("Could not create directory: ", qc_dir)
  }

  f_all <- list.files(path = png_dir, pattern = ".png")
  if (length(f_all) == 0L) {
    stop("No PNG images found in: ", png_dir)
  }

  dato <- NULL
  for (file_name in f_all) {
    interim <- imager::load.image(file.path(png_dir, file_name))
    image_name <- substr(file_name, 0L, nchar(file_name) - 4L)
    interim <- as.data.frame(imager::grayscale(interim))
    names(interim)[3L] <- image_name
    if (is.null(dato)) {
      dato <- interim[, 1:3]
    } else {
      dato <- cbind(dato, interim[, 3L])
      names(dato)[ncol(dato)] <- image_name
    }
  }

  if (sample_number > nrow(dato)) {
    stop("sample_number exceeds the number of available pixels.")
  }
  dato <- dato[sample(seq_len(nrow(dato)), sample_number), ]
  write.csv(dato, file = file.path(qc_dir, "sampledPixels.csv"), row.names = FALSE)

  invisible(dato)
}

# Correct overlap between two channels in one image directory.
alakazam <- function(df = file.path("QC", "sampledPixels.csv"),
                     color = "green",
                     sample_number = 10000,
                     corr_cut = 1,
                     overlapper = "CH1_3",
                     overlappee = "CH1_4",
                     amplification = 2,
                     breck = 10,
                     pool = 4,
                     image_dir = ".") {
  image_dir <- .purgo_directory(image_dir, "image_dir")
  if (!is.character(df) || length(df) != 1L) {
    stop("df must identify one CSV file.")
  }
  df_path <- .purgo_file(df, image_dir)
  if (!file.exists(df_path)) {
    stop("Pixel data file not found: ", df_path)
  }
  if (!is.character(overlapper) || length(overlapper) != 1L || !nzchar(overlapper) ||
        !is.character(overlappee) || length(overlappee) != 1L || !nzchar(overlappee)) {
    stop("overlapper and overlappee must be non-empty channel names.")
  }

  dato <- read.csv(df_path)

  overlapper_matches <- list.files(path = image_dir, pattern = overlapper)
  overlappee_matches <- list.files(path = image_dir, pattern = overlappee)
  if (length(overlapper_matches) == 0L || length(overlappee_matches) == 0L) {
    stop("Could not find TIFF files matching overlapper and overlappee in: ", image_dir)
  }
  overlapper_file <- overlapper_matches[1L]
  overlappee_file <- overlappee_matches[1L]
  overlapper <- substr(overlapper_file, 1L, nchar(overlapper_file) - 4L)
  overlappee <- substr(overlappee_file, 1L, nchar(overlappee_file) - 4L)

  dato <- split(dato, cut(dato[, grepl(overlapper, names(dato))], breaks = breck))
  distro <- NULL
  for (i in dato[1:4]) {
    selected <- head(i[order(i[, grepl(overlappee, names(i))]), ], n = pool)
    if (is.null(distro)) {
      distro <- selected
    } else {
      distro <- rbind(distro, selected)
    }
  }

  dato <- read.csv(df_path)
  line_eq <- suppressMessages(lm(
    distro[, grepl(overlappee, names(distro))] ~
      distro[, grepl(overlapper, names(distro))],
    data = distro
  ))
  line_co <- as.numeric(line_eq$coefficients[2L])
  line_in <- as.numeric(line_eq$coefficients[1L])
  line_pv <- summary(line_eq[[4L]][8L])

  qc_dir <- file.path(image_dir, "QC")
  draft <- ggplot2::ggplot(
    data = dato,
    ggplot2::aes(
      x = dato[, grepl(overlapper, names(dato))],
      y = dato[, grepl(overlappee, names(dato))]
    )
  ) +
    ggplot2::geom_point() +
    ggplot2::theme_classic() +
    ggplot2::ylab(overlappee) +
    ggplot2::xlab(overlapper) +
    ggplot2::geom_abline(
      intercept = line_in,
      slope = line_co,
      color = "red",
      linetype = "dashed",
      size = 1.5
    ) +
    ggplot2::ggtitle(paste0(
      "Initial correlation: ", round(line_co, 2), "x + ", round(line_in, 2),
      ", p-val = ", line_pv
    ))
  suppressMessages(ggplot2::ggsave(
    filename = file.path(qc_dir, paste0("InitialCorr_", overlappee, "x", overlapper, ".png")),
    plot = draft
  ))

  png_dir <- file.path(image_dir, "PNGS")
  interim_ee <- as.data.frame(imager::grayscale(
    imager::load.image(file.path(png_dir, paste0(overlappee, ".png")))
  ))
  interim_er <- as.data.frame(imager::grayscale(
    imager::load.image(file.path(png_dir, paste0(overlapper, ".png")))
  ))
  names(interim_ee)[3L] <- overlappee
  interim <- cbind(interim_ee, interim_er[, 3L])
  names(interim)[4L] <- overlapper

  corrected_name <- paste0(overlappee, "_corrected")
  interim[corrected_name] <- amplification *
    (interim[overlappee] - (interim[overlapper] * line_co))
  interim[corrected_name] <- interim[corrected_name] + abs(min(interim[corrected_name]))

  if (nrow(interim) < 10000L) {
    stop("Fewer than 10000 pixels are available for corrected QC sampling.")
  }
  dato <- interim[sample(seq_len(nrow(interim)), 10000L), ]
  write.csv(dato, file = file.path(qc_dir, "corrPixels.csv"), row.names = FALSE)

  draft <- ggplot2::ggplot(
    data = dato,
    ggplot2::aes(
      x = dato[, grepl(overlapper, names(dato))],
      y = dato[, grepl(corrected_name, names(dato))]
    )
  ) +
    ggplot2::geom_point() +
    ggplot2::theme_classic() +
    ggplot2::ylab(paste0("corrected ", overlappee)) +
    ggplot2::xlab(overlapper) +
    ggplot2::geom_abline(
      intercept = line_in,
      slope = line_co,
      color = "red",
      linetype = "dashed",
      size = 1.5
    ) +
    ggplot2::ggtitle(paste0(
      "Initial correlation: ", round(line_co, 2), "x + ", round(line_in, 2),
      ", p-val = ", line_pv
    ))
  suppressMessages(ggplot2::ggsave(
    filename = file.path(qc_dir, paste0("AdjustedCorr_", overlappee, "x", overlapper, ".png")),
    plot = draft
  ))

  if (color == "green") {
    g <- interim[, 1:2]
    g$x <- abs(g$x - (max(g$x) + 1))
    g$cc <- 2
    g$value <- interim[, 5L]
    r <- g
    r$cc <- 1
    r$value <- 0
    b <- g
    b$cc <- 3
    b$value <- 0
  } else if (color == "red") {
    r <- interim[, 1:2]
    r$x <- abs(r$x - (max(r$x) + 1))
    r$cc <- 1
    r$value <- interim[, 5L]
    g <- r
    g$cc <- 2
    g$value <- 0
    b <- r
    b$cc <- 3
    b$value <- 0
  } else {
    b <- interim[, 1:2]
    b$x <- abs(b$x - (max(b$x) + 1))
    b$cc <- 3
    b$value <- interim[, 5L]
    r <- b
    r$cc <- 1
    r$value <- 0
    g <- b
    g$cc <- 2
    g$value <- 0
  }

  new_dat <- rbind(r, g, b)
  new_plot <- suppressWarnings(imager::as.cimg(new_dat))
  new_plot <- imager::cimg2magick(new_plot, rotate = TRUE)
  png_output <- file.path(png_dir, paste0(overlappee, ".png"))
  tif_output <- file.path(image_dir, paste0(overlappee, ".tif"))
  magick::image_write(new_plot, path = png_output, format = "png")
  magick::image_write(new_plot, path = tif_output, format = "tif")

  invisible(list(
    overlapper = overlapper,
    overlappee = overlappee,
    coefficient = line_co,
    intercept = line_in,
    png = png_output,
    tif = tif_output
  ))
}

# Run purgo processing across every image directory in a project.
purgo <- function(run_type = "none",
                  gen_qc = "y",
                  project_dir = ".",
                  overlapper = NULL,
                  overlappee = NULL) {
  project_dir <- .purgo_directory(project_dir, "project_dir")
  run_type <- tolower(run_type)
  if (length(run_type) != 1L || !run_type %in% c("auto", "manual", "full", "none")) {
    stop("run_type must be one of 'auto', 'manual', 'full', or 'none'.")
  }
  if (run_type == "manual" &&
        (is.null(overlapper) || is.null(overlappee) ||
           length(overlapper) != 1L || length(overlappee) != 1L ||
           !nzchar(overlapper) || !nzchar(overlappee))) {
    stop("Manual mode requires non-empty overlapper and overlappee values.")
  }

  generate_qc <- if (is.logical(gen_qc) && length(gen_qc) == 1L && !is.na(gen_qc)) {
    gen_qc
  } else if (is.character(gen_qc) && length(gen_qc) == 1L && !is.na(gen_qc)) {
    tolower(gen_qc) != "n"
  } else {
    stop("gen_qc must be a single logical value or character response.")
  }

  schema_path <- file.path(project_dir, "schema.csv")
  files_dir <- file.path(project_dir, "files")
  if (!file.exists(schema_path)) {
    stop("No schema.csv found in project_dir: ", project_dir)
  }
  if (!dir.exists(files_dir)) {
    stop("No files directory found in project_dir: ", project_dir)
  }

  scheme <- names(read.csv(schema_path))[2:5]
  if (run_type == "full" && (length(scheme) < 4L || anyNA(scheme))) {
    stop("Full mode requires at least four target columns in schema.csv.")
  }

  processed <- character()
  condition_names <- list.files(path = files_dir)
  for (condition_name in condition_names) {
    condition_dir <- file.path(files_dir, condition_name)
    if (grepl("ijm", condition_name) || grepl(".csv", condition_name) ||
          !dir.exists(condition_dir)) {
      next
    }

    image_names <- list.files(path = condition_dir)
    for (image_name in image_names) {
      image_dir <- file.path(condition_dir, image_name)
      if (image_name == "images" || image_name == "Thumbs.db" ||
            grepl(".csv", image_name) || !dir.exists(image_dir)) {
        next
      }
      if (dir.exists(file.path(image_dir, "originals"))) {
        next
      }

      print(paste0(
        "Running purgo on file ", condition_name, " picture ", image_name, "."
      ))
      abra(image_dir = image_dir)

      if (run_type == "auto") {
        kadabra(image_dir = image_dir)
        alakazam(image_dir = image_dir)
      } else if (run_type == "manual") {
        kadabra(image_dir = image_dir)
        alakazam(
          overlapper = overlapper,
          overlappee = overlappee,
          image_dir = image_dir
        )
      } else if (run_type == "full") {
        kadabra(image_dir = image_dir)
        alakazam(
          overlapper = scheme[4L], overlappee = scheme[3L], color = "red",
          image_dir = image_dir
        )
        alakazam(
          overlapper = scheme[3L], overlappee = scheme[2L], color = "blue",
          image_dir = image_dir
        )
        alakazam(
          overlapper = scheme[2L], overlappee = scheme[1L], color = "green",
          image_dir = image_dir
        )
      }

      if (generate_qc) {
        sampling <- kadabra(image_dir = image_dir)
        qc_dir <- file.path(image_dir, "QC")
        channel_names <- names(sampling)[3:ncol(sampling)]
        for (i in channel_names) {
          for (j in channel_names) {
            existing_plot <- paste0(j, "by", i, ".png")
            if (!existing_plot %in% list.files(path = qc_dir) && i != j) {
              draft <- ggpubr::ggscatter(data = sampling, x = i, y = j)
              ggplot2::ggsave(
                filename = file.path(qc_dir, paste0(i, "by", j, ".png")),
                plot = draft,
                dpi = 300,
                width = 10,
                height = 10,
                units = "cm"
              )
            }
          }
        }
      }

      processed <- c(processed, image_dir)
    }
  }

  invisible(processed)
}

run_purgo_cli <- function(project_dir = ".", prompt = readline) {
  if (!is.function(prompt)) {
    stop("prompt must be a function.")
  }

  ask <- function(message) prompt(message)
  run_type <- tolower(trimws(ask(
    "What should the runtype be (auto/manual/full/none): "
  )))
  if (!nzchar(run_type)) {
    run_type <- "none"
  }

  gen_qc <- trimws(ask("should QC graphs be generated (y/n): "))
  if (!nzchar(gen_qc)) {
    gen_qc <- "y"
  }

  overlapper <- NULL
  overlappee <- NULL
  if (run_type == "manual") {
    project_dir_path <- .purgo_directory(project_dir, "project_dir")
    schema_path <- file.path(project_dir_path, "schema.csv")
    if (!file.exists(schema_path)) {
      stop("No schema.csv found in project_dir: ", project_dir_path)
    }
    scheme <- names(read.csv(schema_path))[2:5]
    print(scheme)
    overlapper <- trimws(ask("Which is the overlapper: "))
    overlappee <- trimws(ask("Which is the overlappee: "))
  }

  invisible(purgo(
    run_type = run_type,
    gen_qc = gen_qc,
    project_dir = project_dir,
    overlapper = overlapper,
    overlappee = overlappee
  ))
}

if (sys.nframe() == 0L) {
  run_purgo_cli()
}
