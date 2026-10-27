#!/usr/bin/Rscript

# The ima_gen() script takes any number of csv files that are within your working directory
# and combines them. Take caution, however, as this script is designed for a particular imagej
# csv outpt.
# It is important that the jetData uses the Set Measurements with these settings active:
#   Area, Standard Deviation, Min & Max gray value, Center of Mass, Mean gray value,
#   Perimenter, and Display label
# Additionally, for the best results, decimal places should be set to 9

# If you have used the YggData macro, this will all be automatically set for you

.ima_gen_path <- function(path, base_dir = ".") {
  if (length(path) != 1L || is.na(path) || !nzchar(path)) {
    stop("A single, non-empty path is required.", call. = FALSE)
  }

  is_absolute <- grepl("^(/|[A-Za-z]:[/\\\\]|\\\\\\\\)", path)
  if (!is_absolute) {
    path <- file.path(base_dir, path)
  }
  normalizePath(path, winslash = "/", mustWork = FALSE)
}

# Assign nuclear number to each ROI that may or may not be localized to the nucleus.
join_r <- function(y, x = "dna.csv", directory = ".", project_dir = ".") {
  base_dir <- .ima_gen_path(directory, project_dir)
  roi_path <- .ima_gen_path(y, base_dir)
  anchor_path <- .ima_gen_path(x, base_dir)
  if (!file.exists(roi_path)) {
    stop("ROI input does not exist: ", roi_path, call. = FALSE)
  }
  if (!file.exists(anchor_path)) {
    stop("Anchor input does not exist: ", anchor_path, call. = FALSE)
  }

  # Files are opened
  tagman <- read.csv(anchor_path)
  fraction <- read.csv(roi_path)
  # The number column is named, if not present
  if (!"Number" %in% names(tagman)) {
    names(tagman)[1] <- "Number"
    write.csv(tagman, file = anchor_path, row.names = FALSE)
  }
  if (!"Number" %in% names(fraction)) {
    names(fraction)[1] <- "Number"
    write.csv(fraction, file = roi_path, row.names = FALSE)
  }

  # Now we skip fractions with empty lists
  if (nrow(fraction) > 0) {
    # Default cellid is given
    fraction$cell <- "unknown"
    # Default distance is given
    tagman$distance <- 0
    fraction$NucDist <- 0

    # The closest nucleus is calculated for each ROI
    for (i in seq_len(nrow(fraction))) {
      tagman$distance <- sqrt((tagman$x - fraction$x[i])^2 + (tagman$y - fraction$y[i])^2)
      closest <- min(tagman$distance)
      fraction$NucDist[i] <- closest
      fraction$cell[i] <- unique(tagman$Number[tagman$distance == closest])
      # If two or more nuclei are picked, it lets you know
      if (length(closest) > 1) {
        print("PING! More than one nucleus at minimum distance. Something is wrong (probably)")
      }
    }
    write.csv(fraction, file = roi_path, row.names = FALSE)
    print(paste0("Saved as ", roi_path))
  } else {
    print(paste0("No data found for file: ", roi_path))
  }

  invisible(roi_path)
}

ima_gen <- function(directory = "./",
                    colorz = FALSE,
                    label_ids = FALSE,
                    peri = FALSE,
                    wc = FALSE,
                    project_dir = ".",
                    condition_id = NULL,
                    image_id = NULL,
                    schema_path = NULL,
                    prompt = readline) {
  if (!is.function(prompt)) {
    stop("prompt must be a function.", call. = FALSE)
  }

  project_dir <- .ima_gen_path(project_dir)
  directory <- .ima_gen_path(directory, project_dir)
  if (!dir.exists(directory)) {
    stop("Image PNGS directory does not exist: ", directory, call. = FALSE)
  }

  # The directory is set and a list of CSV's is generated
  directory_n <- file.path(directory, "Anchor_extraction")
  directory_c <- file.path(directory, "NonAnchor_extraction")
  directory_p <- file.path(directory, "PeriAnchor_extraction")
  image_directory <- dirname(directory)

  if (is.null(image_id)) {
    image_id <- basename(image_directory)
  }
  if (is.null(condition_id)) {
    condition_id <- basename(dirname(image_directory))
  }
  if (is.null(schema_path)) {
    schema_path <- file.path(project_dir, "schema.csv")
    if (!file.exists(schema_path)) {
      legacy_schema_path <- normalizePath(
        file.path(directory, "..", "..", "..", "..", "schema.csv"),
        winslash = "/",
        mustWork = FALSE
      )
      if (file.exists(legacy_schema_path)) {
        schema_path <- legacy_schema_path
      }
    }
  } else {
    schema_path <- .ima_gen_path(schema_path, project_dir)
  }

  if (!dir.exists(directory_n)) {
    stop("Anchor extraction directory does not exist: ", directory_n, call. = FALSE)
  }
  # Generates the list of data files
  filez <- list.files(directory_n, pattern = ".csv")
  filez <- filez[!grepl("all.csv", filez)]
  if (length(filez) == 0L) {
    stop("No anchor CSV inputs were found in: ", directory_n, call. = FALSE)
  }

  # Use an explicit anchor file; otherwise, expect dna or ask the user to define it.
  if ("anchor.csv" %in% filez) {
    anchor_name <- "anchor"
  } else if ("dna.csv" %in% filez) {
    anchor_name <- "dna"
  } else {
    print(filez)
    anchor_name <- prompt("Anchor file not detected. Which file should be the anchor? ")
  }

  # The first file is opened to serve as a template
  # Any file works because they share the same anchor-dependent extraction order.
  cells <- read.csv(file.path(directory_n, filez[1]))

  # If the csv name is not simply the target name (colorz == F),
  # then you will be asked to name a target for each image

  # REMEMBER: One of the targets must be 'dna'
  if (colorz == FALSE) {
    cat("REMEMBER: One of these color's must be match the given anchor name.")
    cat("\n")
    cat(paste0("This file is ", filez[1]))
    cat("\n")
    colo <- prompt(paste0("What color is in ", filez[1], ": "))
  } else {
    # Otherwise, the target is taken from the csv name
    colo <- substr(filez[1], 1, nchar(filez[1]) - 4)
  }

  # The area is artificially increased to make small numbers bigger
  cells$Area <- cells$Area * 100

  # The template's columns are appended with the anchor marker and the 'color
  names(cells) <- paste0(names(cells), "_ANC_", colo)

  # The first column with designated as the cell number which is dictated by the nucleus
  # and is consistent across all downstream applications
  names(cells)[1] <- "Number"

  # The label is removed as this isn't required after 'colo' labeling
  cells <- cells[, -2]

  # The rest of the csv's in the Nuclear directory are opened and cbinded to the first set
  if (length(filez) > 1L) {
    for (i in filez[-1L]) {
      # Already bound csv's (_all tagged) are ignored
      if (!grepl("_all.csv", i)) {
        interim <- read.csv(file.path(directory_n, i))
        if (colorz == FALSE) {
          cat("REMEMBER: One of these color's must match the given anchor name")
          cat("/n")
          cat(paste0("This file is ", i))
          cat("/n")
          colo <- prompt(paste0("What color is in ", i, ": "))
        } else {
          colo <- substr(i, 1, nchar(i) - 4)
        }
        names(interim) <- paste0(names(interim), "_ANC_", colo)
        cells <- cbind(cells, interim[, 4:8], interim[, 11:12], interim[, 15:16])
      }
    }
  }
  # The anchor dependent file is saved
  filnam <- image_id
  xoo <- paste0(filnam, "_ANC_all.csv")
  anchor_output_path <- file.path(directory_n, xoo)
  write.csv(cells, file = anchor_output_path, row.names = FALSE)

  # Now that the nuclear data is collected, the ROI data will be parsed and added
  # First, you navigate to the NonAnchor file folder
  if (wc == FALSE) {
    cat("Parsing anchor-indepedent data...\n")
    if (!dir.exists(directory_c)) {
      stop("Non-anchor extraction directory does not exist: ", directory_c, call. = FALSE)
    }

    # Then you get the list of files
    filez <- list.files(directory_c)
    filez <- filez[!grepl("all.csv", filez)]
    if (length(filez) == 0L) {
      stop("No non-anchor inputs were found in: ", directory_c, call. = FALSE)
    }

    # Every ROI is designated an anchor number and the file is saved
    for (j in filez) {
      join_r(
        file.path(directory_c, j),
        file.path(directory_n, paste0(anchor_name, ".csv"))
      )
    }

    # The anchor-extraction file is opened
    cells <- read.csv(anchor_output_path)

    # Each ROI file is summarized and added to the anchor data
    for (i in filez) {
      # First let smake sure they have labels (they do)
      if (!grepl("_all.csv", i)) {
        if (colorz == FALSE) {
          cat("REMEMBER: One of these color's must match the given anchor'")
          cat("\n")
          cat(paste0("This file is ", i))
          cat("\n")
          colo <- prompt(paste0("What color is in ", i, ": "))
        } else {
          colo <- substr(i, 1, nchar(i) - 4)
        }
      }

      # Temporary dataset is made for easy column labeling and given summary columns
      interim <- cells
      interim$ROI_Num <- 0
      interim$ROI_Area <- 0
      interim$ROI_Mean <- 0
      interim$ROI_Stdev <- 0
      interim$ROI_Mode <- 0
      interim$ROI_Perimeter <- 0
      interim$ROI_IntDen <- 0
      interim$ROI_IntTotal <- 0
      interim$ROI_NucDist <- 0

      # The ROI files is opened as 'fraction'
      if (!grepl("_all.csv", i)) {
        fraction <- read.csv(file.path(directory_c, i))

        # if there are ROI values, it continues
        if (nrow(fraction) > 0L) {
          # Each nucleus number is analyzed for ROI's designated to it
          for (j in seq_len(nrow(interim))) {
            if (interim$Number[j] %in% unique(fraction$cell)) {
              # The ROIs that are assigned to a particular nucleus number are subsetted
              hitz <- fraction[fraction$cell == interim$Number[j], , drop = FALSE]
              if (colo != anchor_name) {
                cat(paste0(
                  "found ", nrow(hitz), " ", colo,
                  " ROIS for anchor number ", interim$Number[j], "\n"
                ))
              }

              # Individual ROI nuclear distances are calculated
              # The number of ROIS assigned to that cell
              interim$ROI_Num[j] <- nrow(hitz)

              # The average average of all the ROIs
              interim$ROI_Area[j] <- mean(hitz$Area) * 100

              # The average Mean intensity of all the ROIs
              interim$ROI_Mean[j] <- mean(hitz$Mean)

              # The average standard deviation of all the ROIs
              interim$ROI_Stdev[j] <- mean(hitz$StdDev)

              # The average mode of pixel intensities of all the ROIs
              interim$ROI_Mode[j] <- mean(hitz$Mode)

              # The average perimeter of all the ROIs
              interim$ROI_Perimeter[j] <- mean(hitz$Perim.)

              # The average integrated mean/density of the ROIs
              interim$ROI_IntDen[j] <- mean(hitz$RawIntDen)

              # The total of the integrated means/densities of all the ROIs
              interim$ROI_IntTotal[j] <- sum(hitz$RawIntDen)

              # The average distance from the nucleus of all the ROIs
              interim$ROI_NucDist[j] <- mean(hitz$NucDist)
            }
          }

          # A target tag is added for later identification
          names(interim) <- paste0(names(interim), "_", colo)
          # The new variables are added to their cell ID
          cells <- cbind(cells, interim[, (ncol(cells) + 1):ncol(interim)])
        }
      }
    }
    whole_cell_output_path <- file.path(directory_c, paste0(filnam, "_WC_all.csv"))
    write.csv(cells, file = whole_cell_output_path, row.names = FALSE)

    # Now begins generating a final ROI file for funnsies
    interim <- read.csv(file.path(directory_c, filez[1]))
    colo <- substr(filez[1], 1, nchar(i) - 2)
    interim$roi <- colo

    if (length(filez) > 1) {
      for (i in filez[2:length(filez)]) {
        if (!grepl("_all.csv", i)) {
          x <- read.csv(file.path(directory_c, i))
          if (nrow(x > 0)) {
            colo <- substr(i, 1, nchar(i) - 4)
            for (j in names(interim)) {
              if (!j %in% names(x)) {
                x[j] <- "NA"
              }
            }
            for (j in names(x)) {
              if (!j %in% names(interim)) {
                interim[j] <- "NA"
              }
            }
            x$roi <- colo
            interim <- rbind(interim, x)
          }
        }
      }
    }

    # Now things get saved
    roi_output_path <- file.path(image_directory, paste0(filnam, "_ROI_all.csv"))
    cells_output_path <- file.path(image_directory, paste0(filnam, "_WN_all.csv"))
    write.csv(interim, file = roi_output_path, row.names = FALSE)
    write.csv(cells, file = cells_output_path, row.names = FALSE)
  }
  # This ends wc collection------------------------------------

  # start of peri collection
  if (peri == FALSE) {
    if (!dir.exists(directory_p)) {
      stop("Peri-anchor extraction directory does not exist: ", directory_p, call. = FALSE)
    }
    # Then you get the list of files
    filez <- list.files(directory_p)
    filez <- filez[!grepl("_all.csv", filez)]
    if (length(filez) == 0L) {
      stop("No peri-anchor inputs were found in: ", directory_p, call. = FALSE)
    }

    # If the ROI-independent files have been called, that's what we want to pull
    if (wc) {
      nucleoid <- file.path(image_directory, paste0(filnam, "_WN_all.csv"))
    } else {
      nucleoid <- anchor_output_path
    }
    print(nucleoid)

    # Every ROI is designated a nucleus number
    for (j in filez) {
      join_r(
        file.path(directory_p, j),
        file.path(directory_n, paste0(anchor_name, ".csv"))
      )
    }
    print("Peris joined")

    # The tagger file is opened
    cells <- read.csv(nucleoid)

    # Each ROI file is summarized and added to the nucleus data
    for (i in filez) {
      if (!grepl("_all.csv", i)) {
        if (colorz == FALSE) {
          cat("REMEMBER: One of these color's must be labeled 'dna'")
          cat("\n")
          cat(paste0("This file is ", i))
          cat("\n")
          colo <- prompt(paste0("What color is in ", i, ": "))
        } else {
          colo <- substr(i, 1, nchar(i) - 4)
        }
      }
      # Temporary dataset is made for easy column labeling
      interim <- cells
      # The ROI files is opened
      if (!grepl("_all.csv", i)) {
        fraction <- read.csv(file.path(directory_p, i))

        # Each nucleus number is analyzed for ROI's designated to it
        for (j in seq_len(nrow(interim))) {
          interim$PERI_Area[j] <- 0
          interim$PERI_Mean[j] <- 0
          interim$PERI_Stdev[j] <- 0
          interim$PERI_Mode[j] <- 0
          interim$PERI_Perimeter[j] <- 0
          interim$PERI_IntDen[j] <- 0
          if (interim$Number[j] %in% unique(fraction$cell)) {
            # The PERI's that are assigned to a particular nucleus number are subsetted
            hitz <- fraction[fraction$cell == interim$Number[j], , drop = FALSE]

            # Individual ROI nuclear distances are calculated
            # The number of ROIS assigned to that cell

            # The average average of all the ROIs
            interim$PERI_Area[j] <- mean(hitz$Area) * 100

            # The average Mean intensity of all the ROIs
            interim$PERI_Mean[j] <- mean(hitz$Mean)

            # The average standard deviation of all the ROIs
            interim$PERI_Stdev[j] <- mean(hitz$StdDev)

            # The average mode of pixel intensities of all the ROIs
            interim$PERI_Mode[j] <- mean(hitz$Mode)

            # The average perimeter of all the ROIs
            interim$PERI_Perimeter[j] <- mean(hitz$Perim.)

            # The average integrated mean/density of the ROIs
            interim$PERI_IntDen[j] <- mean(hitz$RawIntDen)
          }
        }
        # A target tag is added for later identification
        names(interim) <- paste0(names(interim), "_", colo)
        # The new variables are added to their cell ID
        cells <- cbind(cells, interim[, (ncol(cells) + 1):ncol(interim)])
      }
    }
    # Now begins generating a final ROI file for funnies
    interim <- read.csv(file.path(directory_p, filez[1]))
    colo <- substr(filez[1], 1, nchar(filez[1]) - 4)
    interim$peri <- colo
    if (length(filez) > 1L) {
      for (i in filez[-1L]) {
        if (!grepl("_all.csv", i)) {
          x <- read.csv(file.path(directory_p, i))
          colo <- substr(i, 1, nchar(i) - 4)
          for (j in names(interim)) {
            if (!j %in% names(x)) {
              x[j] <- "NA"
            }
          }
          for (j in names(x)) {
            if (!j %in% names(interim)) {
              interim[j] <- "NA"
            }
          }
          x$peri <- colo
          interim <- rbind(interim, x)
        }
      }
    }

    # Now things get saved
    peri_output_path <- file.path(image_directory, paste0(filnam, "_PERI_all.csv"))
    peri_nuclear_output_path <- file.path(directory_p, paste0(filnam, "_PN_all.csv"))
    cells_output_path <- file.path(image_directory, paste0(filnam, "_WN_all.csv"))
    write.csv(interim, file = peri_output_path, row.names = FALSE)
    write.csv(cells, file = peri_nuclear_output_path, row.names = FALSE)
    write.csv(cells, file = cells_output_path, row.names = FALSE)
  }
  #-------------------------------------------------------------------------
  # This ends the PERI_all collection

  # #Now that we have a nuclear dataset and a WC dataset, we might as well put it all together...
  # if(wc==T & peri==T){
  # }
  # if(wc == T){
  # }
  # if(peri==T){
  # }
  cells_output_path <- file.path(image_directory, paste0(filnam, "_WN_all.csv"))
  if (wc || peri) {
    cells <- read.csv(cells_output_path)
  } else {
    cells <- read.csv(anchor_output_path)
  }

  if (peri == FALSE) {
    think_tank <- names(cells)[grepl("PERI_Area", names(cells))]
    for (z in think_tank) {
      zim <- strsplit(z, "Area")[[1]][2]
      tic <- names(cells)[grepl("Area_ANC", names(cells))][1]
      cells[paste0("PERI_SubArea", zim)] <- cells[paste0("PERI_Area", zim)] - cells[tic]
      tic <- paste0("IntDen_ANC", zim)
      cells[paste0("PERI_SubIntDen", zim)] <- cells[paste0("PERI_IntDen", zim)] - cells[tic]
      cells[paste0("PERI_SubMean", zim)] <-
        (100 * cells[paste0("PERI_SubIntDen", zim)]) /
        cells[paste0("PERI_SubArea", zim)]
    }
  }

  # Add metadata and rename NumberIDs because this is not assumed to be a stitched image.
  cells$image <- image_id
  cells$Number <- paste0(image_id, "_", cells$Number)
  if (!file.exists(schema_path)) {
    stop(
      "Schema file does not exist: ", schema_path,
      ". Set 'schema_path' or 'project_dir' explicitly.",
      call. = FALSE
    )
  }
  schema <- read.csv(schema_path)
  if (condition_id %in% unique(schema$notebook_id)) {
    hit <- schema[schema$notebook_id == condition_id, , drop = FALSE]
  } else if ("name_id" %in% names(schema)) {
    hit <- schema[schema$name_id == condition_id, , drop = FALSE]
  } else {
    hit <- schema[FALSE, , drop = FALSE]
  }
  if (nrow(hit) == 0L) {
    stop(
      "No schema row matched condition '", condition_id,
      "' by notebook_id or name_id in: ", schema_path,
      call. = FALSE
    )
  }
  if (nrow(hit) > 1) {
    print("More than one scheme match found.")
    print(paste0(
      "The first match will be used but it is recommended that your ammend ",
      "the schema file to avoid other issues"
    ))
    hit <- hit[1, ]
  }
  if (ncol(hit) >= 7L) {
    for (ab in seq.int(7L, ncol(hit))) {
      cells[names(hit)[ab]] <- hit[ab]
      if (peri == FALSE) {
        peris <- read.csv(peri_output_path)
        peris[names(hit)[ab]] <- hit[ab]
        write.csv(peris, file = peri_output_path, row.names = FALSE)
      }
      if (wc == FALSE) {
        rois <- read.csv(roi_output_path)
        rois[names(hit)[ab]] <- hit[ab]
        write.csv(rois, file = roi_output_path, row.names = FALSE)
      }
    }
  }
  # And finally save the whole thing
  write.csv(cells, file = cells_output_path, row.names = FALSE)

  output_paths <- c(anchor = anchor_output_path, cells = cells_output_path)
  if (wc == FALSE) {
    output_paths <- c(output_paths, roi = roi_output_path, wholeCell = whole_cell_output_path)
  }
  if (peri == FALSE) {
    output_paths <- c(output_paths, peri = peri_output_path, periNuclear = peri_nuclear_output_path)
  }
  invisible(output_paths)
}

run_ima_gen <- function(project_dir = ".", files_dir = NULL, schema_path = NULL, ...) {
  project_dir <- .ima_gen_path(project_dir)
  if (!dir.exists(project_dir)) {
    stop("Project directory does not exist: ", project_dir, call. = FALSE)
  }
  if (is.null(files_dir)) {
    files_dir <- file.path(project_dir, "files")
  } else {
    files_dir <- .ima_gen_path(files_dir, project_dir)
  }
  if (!dir.exists(files_dir)) {
    stop(
      "Input directory does not exist: ", files_dir,
      ". Set 'files_dir' or 'project_dir' to the MicroCyte project location.",
      call. = FALSE
    )
  }
  if (is.null(schema_path)) {
    schema_path <- file.path(project_dir, "schema.csv")
  } else {
    schema_path <- .ima_gen_path(schema_path, project_dir)
  }
  if (!file.exists(schema_path)) {
    stop("Schema file does not exist: ", schema_path, call. = FALSE)
  }

  extra_args <- list(...)
  reserved_args <- c(
    "directory", "project_dir", "condition_id", "image_id", "schema_path", "peri", "wc"
  )
  duplicate_args <- intersect(names(extra_args), reserved_args)
  if (length(duplicate_args) > 0L) {
    stop(
      "run_ima_gen() controls these arguments: ", paste(duplicate_args, collapse = ", "),
      ". Pass image-specific calls directly to ima_gen() instead.",
      call. = FALSE
    )
  }

  x_list <- list.dirs(files_dir, full.names = FALSE, recursive = FALSE)
  x_list <- x_list[!grepl("ijm", x_list) & x_list != "Thumbs.db"]
  if (length(x_list) == 0L) {
    stop("No condition directories were found in: ", files_dir, call. = FALSE)
  }

  results <- list()
  for (xL in x_list) {
    condition_dir <- file.path(files_dir, xL)
    print(paste0("Running ImaGen on the ", xL, " folder:"))
    y_list <- list.dirs(condition_dir, full.names = FALSE, recursive = FALSE)
    y_list <- y_list[y_list != "Thumbs.db"]
    if (length(y_list) == 0L) {
      stop("No image directories were found in condition directory: ", condition_dir, call. = FALSE)
    }

    for (yL in y_list) {
      image_png_dir <- file.path(condition_dir, yL, "PNGS")
      if (!dir.exists(image_png_dir)) {
        stop("Image PNGS directory does not exist: ", image_png_dir, call. = FALSE)
      }
      anchor_directory <- file.path(image_png_dir, "Anchor_extraction")
      if (!dir.exists(anchor_directory)) {
        stop("Anchor extraction directory does not exist: ", anchor_directory, call. = FALSE)
      }
      peri_go <- dir.exists(file.path(image_png_dir, "PeriAnchor_extraction"))
      wc_go <- dir.exists(file.path(image_png_dir, "NonAnchor_extraction"))
      combined_exists <- length(list.files(
        path = anchor_directory,
        pattern = "ANC_all.csv"
      )) > 0L

      ima_gen_args <- c(
        list(
          directory = image_png_dir,
          peri = peri_go,
          wc = wc_go,
          project_dir = project_dir,
          condition_id = xL,
          image_id = yL,
          schema_path = schema_path
        ),
        extra_args
      )
      results[[file.path(xL, yL)]] <- do.call(ima_gen, ima_gen_args)
      if (combined_exists) {
        print(paste0("Combined data already detected for: ", yL))
      } else {
        print(paste0("Combining data from image ", yL))
      }
    }
  }

  invisible(results)
}

if (sys.nframe() == 0L) {
  run_ima_gen()
}
