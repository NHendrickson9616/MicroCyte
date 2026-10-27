#!/usr/bin/Rscript

.rename_project_dir <- function(project_dir) {
  if (!is.character(project_dir) || length(project_dir) != 1L || !dir.exists(project_dir)) {
    stop("project_dir must identify an existing directory.")
  }
  normalizePath(project_dir, mustWork = TRUE)
}

.rename_file <- function(from, to) {
  if (identical(from, to)) {
    return(FALSE)
  }
  if (file.exists(to)) {
    stop("Refusing to overwrite existing file while renaming '", from, "' to '", to, "'.")
  }
  if (!file.rename(from, to)) {
    stop("Could not rename '", from, "' to '", to, "'.")
  }
  TRUE
}

.rename_target <- function(stem, target_names) {
  matches <- target_names[vapply(
    target_names,
    function(target_name) grepl(target_name, stem, fixed = TRUE),
    logical(1)
  )]
  if (length(matches) == 0L) {
    return(NULL)
  }

  # Prefer the most specific channel name, e.g. CH10 rather than CH1.
  matches[which.max(nchar(matches))]
}

re_name <- function(name_number = 6, project_dir = ".") {
  if (length(name_number) != 1L || !is.numeric(name_number) || is.na(name_number)) {
    stop("name_number must be one numeric value.")
  }

  project_dir <- .rename_project_dir(project_dir)
  schema_path <- file.path(project_dir, "schema.csv")
  files_dir <- file.path(project_dir, "files")
  if (!file.exists(schema_path)) {
    stop("No schema.csv found in project_dir: ", project_dir)
  }
  if (!dir.exists(files_dir)) {
    stop("No files directory found in project_dir: ", project_dir)
  }

  schema <- read.csv(schema_path)
  if (!"notebook_id" %in% names(schema)) {
    stop("schema.csv must contain a notebook_id column.")
  }

  last_target_column <- min(as.integer(name_number), ncol(schema))
  if (last_target_column < 2L) {
    stop("schema.csv must contain at least one target column after notebook_id.")
  }
  name_targets <- names(schema)[seq.int(2L, last_target_column)]
  renamed <- character()
  condition_names <- list.files(path = files_dir)
  for (condition_name in condition_names) {
    condition_dir <- file.path(files_dir, condition_name)
    if (grepl(".ijm", x = condition_name) || !dir.exists(condition_dir)) {
      next
    }

    if (condition_name %in% unique(schema$notebook_id)) {
      hit <- schema[schema$notebook_id == condition_name, , drop = FALSE][1L, ]
    } else if ("name_id" %in% names(schema) &&
                 condition_name %in% unique(schema$name_id)) {
      hit <- schema[schema$name_id == condition_name, , drop = FALSE][1L, ]
    } else {
      warning("No schema row found for directory: ", condition_name)
      next
    }

    image_names <- list.files(path = condition_dir)
    for (image_name in image_names) {
      image_dir <- file.path(condition_dir, image_name)
      if (!dir.exists(image_dir)) {
        next
      }

      tif_names <- list.files(path = image_dir, pattern = "[.]tif$", ignore.case = TRUE)
      for (file_name in tif_names) {
        stem <- tools::file_path_sans_ext(file_name)
        target_name <- .rename_target(stem, name_targets)
        source_path <- file.path(image_dir, file_name)

        if (grepl("Overlay[.]tif$", file_name, ignore.case = TRUE)) {
          destination_path <- file.path(image_dir, "overlay.tif")
          if (.rename_file(source_path, destination_path)) {
            renamed <- c(renamed, destination_path)
          }
        } else if (!is.null(target_name)) {
          new_name <- paste0(as.character(hit[[target_name]][1L]), ".tif")
          destination_path <- file.path(image_dir, new_name)
          if (.rename_file(source_path, destination_path)) {
            renamed <- c(renamed, destination_path)
          }
        }
      }

      png_dir <- file.path(image_dir, "PNGS")
      if (dir.exists(png_dir)) {
        png_names <- list.files(path = png_dir, pattern = "[.]png$", ignore.case = TRUE)
        for (file_name in png_names) {
          stem <- tools::file_path_sans_ext(file_name)
          target_name <- .rename_target(stem, name_targets)
          source_path <- file.path(png_dir, file_name)

          if (grepl("Overlay[.]png$", file_name, ignore.case = TRUE)) {
            destination_path <- file.path(png_dir, "overlay.png")
            if (.rename_file(source_path, destination_path)) {
              renamed <- c(renamed, destination_path)
            }
          } else if (!is.null(target_name)) {
            new_name <- paste0(as.character(hit[[target_name]][1L]), ".png")
            destination_path <- file.path(png_dir, new_name)
            if (.rename_file(source_path, destination_path)) {
              renamed <- c(renamed, destination_path)
            }
          }
        }
      }
    }
  }

  invisible(renamed)
}

# Restore original TIFF files and remove generated output from each image directory.
reset_names <- function(project_dir = ".") {
  project_dir <- .rename_project_dir(project_dir)
  files_dir <- file.path(project_dir, "files")
  if (!dir.exists(files_dir)) {
    stop("No files directory found in project_dir: ", project_dir)
  }

  reset_dirs <- character()
  condition_names <- list.files(path = files_dir)
  for (condition_name in condition_names) {
    condition_dir <- file.path(files_dir, condition_name)
    if (grepl(".ijm", x = condition_name) || grepl(".csv", x = condition_name) ||
          !dir.exists(condition_dir)) {
      next
    }

    image_names <- list.files(path = condition_dir)
    for (image_name in image_names) {
      image_dir <- file.path(condition_dir, image_name)
      if (grepl(".ijm", x = image_name) || grepl(".csv", x = image_name) ||
            !dir.exists(image_dir)) {
        next
      }

      originals_dir <- file.path(image_dir, "originals")
      if (!dir.exists(originals_dir)) {
        next
      }

      original_names <- list.files(path = originals_dir)
      if (length(original_names) == 0L) {
        next
      }

      tif_names <- list.files(path = image_dir, pattern = "tif")
      if (length(tif_names) > 0L) {
        file.remove(file.path(image_dir, tif_names))
      }

      generated_names <- list.files(path = image_dir)
      generated_names <- generated_names[generated_names != "originals"]
      for (generatedName in generated_names) {
        unlink(file.path(image_dir, generatedName), recursive = TRUE)
      }

      file.copy(
        from = file.path(originals_dir, original_names),
        to = image_dir
      )
      unlink(originals_dir, recursive = TRUE)
      reset_dirs <- c(reset_dirs, image_dir)
    }
  }

  invisible(reset_dirs)
}
