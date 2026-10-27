#!/usr/bin/Rscript

.concato_path <- function(path, project_dir = ".") {
  if (length(path) != 1L || is.na(path) || !nzchar(path)) {
    stop("A single, non-empty path is required.", call. = FALSE)
  }

  is_absolute <- grepl("^(/|[A-Za-z]:[/\\\\]|\\\\\\\\)", path)
  if (!is_absolute) {
    path <- file.path(project_dir, path)
  }
  normalizePath(path, winslash = "/", mustWork = FALSE)
}

concato <- function(dirz = "files",
                    count = FALSE,
                    normalization = FALSE,
                    project_dir = ".") {
  files_dir <- .concato_path(dirz, project_dir)
  if (!dir.exists(files_dir)) {
    stop(
      "Input directory does not exist: ", files_dir,
      ". Set 'dirz' or 'project_dir' to the MicroCyte project location.",
      call. = FALSE
    )
  }

  x_list <- list.dirs(files_dir, full.names = FALSE, recursive = FALSE)
  x_list <- x_list[!grepl("ijm", x_list) & x_list != "Thumbs.db"]
  if (length(x_list) == 0L) {
    stop("No condition directories were found in: ", files_dir, call. = FALSE)
  }

  min_num <- NULL
  output_paths <- character()

  for (xL in x_list) {
    print(paste0("Concatenating condition file ", xL))
    condition_dir <- file.path(files_dir, xL)
    y_list <- list.dirs(condition_dir, full.names = FALSE, recursive = FALSE)
    y_list <- y_list[y_list != "Thumbs.db"]
    if (length(y_list) == 0L) {
      stop("No image directories were found in condition directory: ", condition_dir, call. = FALSE)
    }

    cells <- NULL
    rois <- NULL
    peris <- NULL
    roi_go <- FALSE
    peri_go <- FALSE

    for (yL in y_list) {
      image_dir <- file.path(condition_dir, yL)
      z_path <- file.path(image_dir, paste0(yL, "_WN_all.csv"))
      r_path <- file.path(image_dir, paste0(yL, "_ROI_all.csv"))
      p_path <- file.path(image_dir, paste0(yL, "_PERI_all.csv"))

      if (!file.exists(z_path)) {
        stop(
          "Missing combined image input: ", z_path,
          ". Run ima_gen() for image '", yL, "' before concato().",
          call. = FALSE
        )
      }

      if (is.null(cells)) {
        cells <- read.csv(z_path)
      } else {
        interim <- read.csv(z_path)
        for (i in names(cells)[!names(cells) %in% names(interim)]) {
          interim[i] <- 0
        }
        for (i in names(interim)[!names(interim) %in% names(cells)]) {
          cells[i] <- 0
        }
        cells <- rbind(cells, interim)
      }

      if (file.exists(r_path)) {
        roi_go <- TRUE
        if (is.null(rois)) {
          rois <- read.csv(r_path)
        } else {
          rinterim <- read.csv(r_path)
          for (i in names(rois)[!names(rois) %in% names(rinterim)]) {
            rinterim[i] <- 0
          }
          for (i in names(rinterim)[!names(rinterim) %in% names(rois)]) {
            rois[i] <- 0
          }
          rois <- rbind(rois, rinterim)
        }
      }

      if (file.exists(p_path)) {
        peri_go <- TRUE
        if (is.null(peris)) {
          peris <- read.csv(p_path)
        } else {
          pinterim <- read.csv(p_path)
          for (i in names(peris)[!names(peris) %in% names(pinterim)]) {
            pinterim[i] <- 0
          }
          for (i in names(pinterim)[!names(pinterim) %in% names(peris)]) {
            peris[i] <- 0
          }
          peris <- rbind(peris, pinterim)
        }
      }
    }

    if (normalization == FALSE) {
      norm_targets <- names(cells)[grepl("Mean_", names(cells)) | grepl("IntDen", names(cells))]
      for (imageNum in unique(cells$image)) {
        for (target_name in norm_targets) {
          cells[target_name][cells$image == imageNum, ] <-
            cells[target_name][cells$image == imageNum, ] /
            mean(cells[target_name][cells$image == imageNum, ])
          cells[target_name][cells$image == imageNum, ] <-
            cells[target_name][cells$image == imageNum, ] +
            min(cells[target_name][cells$image == imageNum, ])
          quantos <- quantile(unlist(cells[target_name][
            cells$image == imageNum &
              !is.na(cells[target_name][cells$image == imageNum, ]),
          ]))[2]
          cells[target_name][cells$image == imageNum, ] <-
            cells[target_name][cells$image == imageNum, ] / quantos
        }
      }
    }

    cells_path <- file.path(condition_dir, paste0(xL, "_all.csv"))
    write.csv(cells, file = cells_path, row.names = FALSE)
    output_paths <- c(output_paths, cells_path)

    if (roi_go == FALSE) {
      roi_path <- file.path(condition_dir, paste0(xL, "_roi_all.csv"))
      write.csv(rois, file = roi_path, row.names = FALSE)
      output_paths <- c(output_paths, roi_path)
    }
    if (peri_go == FALSE) {
      peri_path <- file.path(condition_dir, paste0(xL, "_peri_all.csv"))
      write.csv(peris, file = peri_path, row.names = FALSE)
      output_paths <- c(output_paths, peri_path)
    }

    if (is.null(min_num) || nrow(cells) < min_num) {
      min_num <- nrow(cells)
    }
  }

  if (count == FALSE) {
    cat(paste0("The minimum cell count recorded is: ", min_num))
    cat("\n")
  }

  invisible(output_paths)
}
