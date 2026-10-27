#!/usr/bin/Rscript

.unite_path <- function(path, project_dir = ".") {
  if (length(path) != 1L || is.na(path) || !nzchar(path)) {
    stop("A single, non-empty path is required.", call. = FALSE)
  }

  is_absolute <- grepl("^(/|[A-Za-z]:[/\\\\]|\\\\\\\\)", path)
  if (!is_absolute) {
    path <- file.path(project_dir, path)
  }
  normalizePath(path, winslash = "/", mustWork = FALSE)
}

.unite_output_name <- function(filename) {
  if (length(filename) != 1L || is.na(filename) || !nzchar(filename)) {
    stop("The output filename must be a single non-empty value.", call. = FALSE)
  }
  if (filename %in% c(".", "..") || grepl("[<>:\"/\\\\|?*]", filename)) {
    stop(
      "The grouping value cannot be used as a cross-platform filename: ", filename,
      call. = FALSE
    )
  }
  filename
}

.unite_sample <- function(data, sample_size, condition_name) {
  if (nrow(data) < sample_size) {
    stop(
      "Condition '", condition_name, "' contains ", nrow(data),
      " rows, fewer than the requested sample size of ", sample_size, ".",
      call. = FALSE
    )
  }
  data[sample.int(nrow(data), sample_size), , drop = FALSE]
}

unite <- function(dirz = "files",
                  full = FALSE,
                  pattern = FALSE,
                  filename = "Experiment",
                  sample_size = 1000,
                  project_dir = ".",
                  output_dir = NULL,
                  fixed = FALSE) {
  if (length(full) != 1L || !is.logical(full) || is.na(full)) {
    stop("full must be TRUE or FALSE.", call. = FALSE)
  }
  if (length(fixed) != 1L || !is.logical(fixed) || is.na(fixed)) {
    stop("fixed must be TRUE or FALSE.", call. = FALSE)
  }
  if (!full && (length(sample_size) != 1L || !is.numeric(sample_size) ||
                  is.na(sample_size) || sample_size < 1 ||
                  sample_size != as.integer(sample_size))) {
    stop("sample_size must be a positive whole number.", call. = FALSE)
  }

  files_dir <- .unite_path(dirz, project_dir)
  if (!dir.exists(files_dir)) {
    stop(
      "Input directory does not exist: ", files_dir,
      ". Set 'dirz' or 'project_dir' to the MicroCyte project location.",
      call. = FALSE
    )
  }

  if (is.null(output_dir)) {
    output_dir <- file.path(project_dir, "data")
  } else {
    output_dir <- .unite_path(output_dir, project_dir)
  }
  output_dir <- normalizePath(output_dir, winslash = "/", mustWork = FALSE)
  if (!dir.exists(output_dir)) {
    stop(
      "Output directory does not exist: ", output_dir,
      ". Create it or set 'output_dir' explicitly.",
      call. = FALSE
    )
  }

  y_list <- list.dirs(files_dir, full.names = FALSE, recursive = FALSE)
  if (!identical(pattern, FALSE)) {
    if (length(pattern) != 1L || !is.character(pattern) || is.na(pattern) || !nzchar(pattern)) {
      stop("pattern must be FALSE or one non-empty character value.", call. = FALSE)
    }
    y_list <- y_list[grepl(pattern, y_list, fixed = fixed)]
    filename <- .unite_output_name(pattern)
    print(paste0("Uniting sample datasets containing the phrase: ", pattern))
  } else {
    filename <- .unite_output_name(filename)
    print(paste0("Uniting all sample datasets found in '", files_dir, "'."))
  }
  if (length(y_list) == 0L) {
    stop("No matching condition directories were found in: ", files_dir, call. = FALSE)
  }

  cells <- NULL
  peris <- NULL
  rois <- NULL

  for (yL in y_list) {
    condition_dir <- file.path(files_dir, yL)
    all_cells_path <- file.path(condition_dir, paste0(yL, "_all_cells.csv"))
    all_path <- file.path(condition_dir, paste0(yL, "_all.csv"))
    if (file.exists(all_cells_path)) {
      z_path <- all_cells_path
    } else if (file.exists(all_path)) {
      z_path <- all_path
    } else {
      stop(
        "Missing concatenated condition input in ", condition_dir,
        ". Expected '", basename(all_cells_path), "' or '", basename(all_path), "'.",
        call. = FALSE
      )
    }

    if (is.null(cells)) {
      cells <- read.csv(z_path)
      if (!full) {
        cells <- .unite_sample(cells, sample_size, yL)
      }
    } else {
      interim <- read.csv(z_path)
      if (!full) {
        interim <- .unite_sample(interim, sample_size, yL)
      }

      if (ncol(interim) != ncol(cells)) {
        for (i in names(cells)[!names(cells) %in% names(interim)]) {
          interim[i] <- 0
        }
        for (i in names(interim)[!names(interim) %in% names(cells)]) {
          cells[i] <- 0
        }
      }
      cells <- rbind(cells, interim)
    }

    peri_list <- list.files(condition_dir, pattern = "peri_all", full.names = TRUE)
    roi_list <- list.files(condition_dir, pattern = "roi_all", full.names = TRUE)

    if (length(peri_list) != 0L) {
      if (is.null(peris)) {
        peris <- read.csv(peri_list[1])
      } else {
        peri_interim <- read.csv(peri_list[1])
        if (ncol(peri_interim) != ncol(peris)) {
          for (i in names(peris)[!names(peris) %in% names(peri_interim)]) {
            peri_interim[i] <- 0
          }
          for (i in names(peri_interim)[!names(peri_interim) %in% names(peris)]) {
            peris[i] <- 0
          }
        }
        peris <- rbind(peris, peri_interim)
      }
    }

    if (length(roi_list) != 0L) {
      if (is.null(rois)) {
        rois <- read.csv(roi_list[1])
      } else {
        roi_interim <- read.csv(roi_list[1])
        if (ncol(roi_interim) != ncol(rois)) {
          for (i in names(rois)[!names(rois) %in% names(roi_interim)]) {
            roi_interim[i] <- 0
          }
          for (i in names(roi_interim)[!names(roi_interim) %in% names(rois)]) {
            rois[i] <- 0
          }
        }
        rois <- rbind(rois, roi_interim)
      }
    }
  }

  if ("log2_dna" %in% names(cells)) {
    cells_path <- file.path(output_dir, paste0(filename, "_all_cells.csv"))
  } else {
    cells_path <- file.path(output_dir, paste0(filename, "_all.csv"))
  }
  write.csv(cells, file = cells_path, row.names = FALSE)
  output_paths <- cells_path

  if (!is.null(peris)) {
    peri_path <- file.path(output_dir, paste0(filename, "_peri_all.csv"))
    write.csv(peris, file = peri_path, row.names = FALSE)
    output_paths <- c(output_paths, peri_path)
  }

  if (!is.null(rois)) {
    roi_path <- file.path(output_dir, paste0(filename, "_roi_all.csv"))
    write.csv(rois, file = roi_path, row.names = FALSE)
    output_paths <- c(output_paths, roi_path)
  }

  invisible(output_paths)
}
