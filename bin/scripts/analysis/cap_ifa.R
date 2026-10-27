# These functions are used for Core Analysis of Plaque IFA (CAP-IFA). After gating the
# dataset as infected or not infected (`variable`), you may call `foci_caller` to generate
# plaque IDs, a list of infected nearest neighbors within the "wiggle" room, and its length.
# An optional third function collects plaque data independently of the individual cell data.

suppressPackageStartupMessages(library(hash))

# This first function generates the list of nearby infected cells for every cell in the data.
# It is imperative that `df` contains cells from a single image. If the data contain cells
# from multiple images, subset them by image ID and run the function on each subset separately.
# The second function, which condenses the lists based on overlap, is called automatically.
foci_caller <- function(df = cells,
                        variable = "n",
                        variable_subset = "Positive",
                        variable_intensity = "Mean_ANC_n",
                        midpointer_x = "XM_ANC_dna",
                        midpointer_y = "YM_ANC_dna",
                        wiggle = 2,
                        overlap_req = 2,
                        file_name = "NN_all_cells.csv") {
  cat("Running foci_caller...\n")

  # Collect only the infected cells
  foci_cells <- df[df[variable] == variable_subset, ]

  # Add some place-holder names
  df["focus"] <- 0
  df["NNs"] <- 0
  df["F_ids"] <- "NA"

  # Generate distances from each infected cell to every cell in the image.
  for (a in seq_len(nrow(foci_cells))) {
    # Get the cell ID
    cell_id_a <- foci_cells[a, ]$Number
    # create a dummy dataset to work the math on
    dist_list <- df
    # Get the X and Y coordinates of the infected cell
    x_midpoint <- as.numeric(foci_cells[a, ][midpointer_x])
    y_midpoint <- as.numeric(foci_cells[a, ][midpointer_y])
    # Calculate distances from the infected cell and store them in the dummy data set.
    dist_list["distance"] <- sqrt(
      (dist_list[midpointer_x] - x_midpoint)**2 +
        (dist_list[midpointer_y] - y_midpoint)**2
    )
    # Order the dummy dataset  by the distances from closest to furthest
    dist_list <- dist_list[order(dist_list$distance), ]
    # Generate a misses variable to use the wiggle variable
    misses <- 0
    focus <- dist_list[dist_list$Number == cell_id_a, ]

    # Now we go through the sorted dummy list to find a list of the closest, infected cells
    for (b in seq_len(nrow(dist_list[dist_list$Number != cell_id_a, ]))) {
      # Starting with itself, add the cell's number ID if it is infected.
      if (dist_list[b, ][variable] == variable_subset) {
        focus <- rbind(focus, dist_list[b, ])

        # Increase misses for each uninfected cell until it surpasses the wiggle value.
        # If this happens, the iteration through the sorted dummy dataset is ended and
        # the length and composition of the infected list is added to the original dataset
      } else {
        misses <- misses + 1
        if (misses > wiggle) {
          focus_size <- nrow(focus)
          df[df["Number"] == cell_id_a, ]$NNs <- focus_size

          f_id <- ""
          for (c in unique(focus$Number)) {
            f_id <- paste(f_id, c, sep = " ")
          }
          f_id <- paste0(f_id, " ")
          df[df["Number"] == cell_id_a, ]$F_ids <- f_id
          rm(focus)
          break
        }
      }
    }
  }
  cat("foci called\n")
  # The new dataset is saved
  write.csv(df, paste0("data/", file_name), row.names = FALSE)
  # The next function will condense the lists of infected neighbors based on overlapping number IDs
  foci_compact(
    file_name = file_name,
    variable = variable,
    variable_subset = variable_subset,
    variable_intensity = variable_intensity,
    overlap = overlap_req
  )
}

foci_compact <- function(file_name = "assigned.csv",
                         variable = "n",
                         variable_subset = "Positive",
                         variable_intensity = "Mean_ANC_n",
                         overlap = 2) {
  cat(paste0("Running foci_compact from ", file_name, ".\n"))
  cells_new <- read.csv(paste0("data/", file_name))
  # First, we get the cells that are infected
  f_list <- cells_new[cells_new[variable] == variable_subset, ]
  # Then we'll order them from great to lowest
  f_list <- f_list[order(f_list$NNs, decreasing = TRUE), ]
  # Then we create the hash to store the focal numbers and initialize the count
  focal_hash <- hash()
  foci_count <- 1
  cat("Hash generated...\n")

  for (a in seq_len(nrow(f_list))) {
    found <- FALSE
    # Starting with the highest NNs value, split NNs by the space delimiter.
    f_ids <- strsplit(
      as.character(f_list[order(f_list$NNs, decreasing = TRUE), ][a, ]$F_ids),
      split = " "
    )
    f_ids <- f_ids[1][[1]][2:length(f_ids[1][[1]])]
    # If no hash entries are detected, the first infected cell NN becomes the first focus entry
    if (length(focal_hash) < 1) {
      focal_hash[[as.character(foci_count)]] <- f_ids
      foci_count <- foci_count + 1
      found <- TRUE
      # Otherwise, the entries are compared the the previously entered hashes
    } else {
      hit <- 0
      for (b in f_ids) {
        for (c in seq_along(focal_hash)) {
          if (b %in% focal_hash[[as.character(c)]]) {
            hit <- hit + 1
            if (hit >= overlap) {
              focal_hash[[as.character(c)]] <- unique(append(focal_hash[[as.character(c)]], f_ids))
              found <- TRUE
              break
            }
          }
        }
      }
    }
    if (found == FALSE) {
      focal_hash[[as.character(foci_count)]] <- f_ids
      foci_count <- foci_count + 1
    }
  }

  cat("Finding overlapping hashes...\n")
  for (a in ls(focal_hash)) {
    found <- FALSE
    b <- focal_hash[[a]]
    for (c in ls(focal_hash)) {
      if (a != c) {
        if (!is.na(as.numeric(table(b %in% focal_hash[[c]])["TRUE"]) >= overlap) &&
              found == FALSE) {
          if (as.numeric(table(b %in% focal_hash[[c]])["TRUE"]) >= overlap) {
            print(paste(
              "matching", a, "with", c, "because they overlap by",
              as.numeric(table(b %in% focal_hash[[c]])["TRUE"])
            ))
            focal_hash[[c]] <- unique(append(focal_hash[[c]], b))
            print(length(focal_hash[[c]]))
            del(a, focal_hash)
            found <- TRUE
          }
        }
      }
      if (found == TRUE) {
        break
      }
    }
  }

  cat("Determining probabilities...\n")
  cells_new$focus <- 0
  cells_new$focal_size <- 0
  for (a in ls(focal_hash)) {
    c <- length(unique(focal_hash[[a]]))
    for (b in focal_hash[[a]]) {
      cells_new[cells_new$Number == b, ]$focus <- a
    }
    cells_new[cells_new["focus"] == a, ]["focal_size"] <- c
  }
  # Find the right name
  uni_x_name <- names(cells_new)[grepl("X_ANC", names(cells_new))][1]
  uni_y_name <- names(cells_new)[grepl("Y_ANC", names(cells_new))][1]


  cells_new$distance_to_center <- NA
  cells_new$intensity_rank <- 0
  cells_new$distance_rank <- 0
  cells_new$nn_rank <- 0

  for (a in unique(cells_new$focus)) {
    # Find the median x and y values and calculate distance
    uni_x <- exp(mean(log(cells_new[cells_new$focus == a, ][, uni_x_name])))
    uni_y <- exp(mean(log(cells_new[cells_new$focus == a, ][, uni_y_name])))

    ticket <- subset(cells_new, focus == a)
    ticket$distance_to_center <- sqrt(
      (ticket[, uni_x_name] - uni_x)^2 + (ticket[, uni_y_name] - uni_y)^2
    )
    ticket <- ticket[order(ticket$distance_to_center), ]
    ticket$distance_rank <- seq_len(nrow(ticket))

    ticket <- ticket[order(-ticket$NNs), ]
    ticket$nn_rank <- seq_len(nrow(ticket))
    if (max(ticket$nn_rank) > 1 && unique(ticket$focus) != 0) {
      for (rowrank in seq_len(nrow(ticket) - 1)) {
        if (ticket[rowrank, ]$NNs == ticket[rowrank + 1, ]$NNs) {
          ticket[rowrank + 1, ]$nn_rank <- ticket[rowrank, ]$nn_rank
        } else {
          ticket[rowrank + 1, ]$nn_rank <- ticket[rowrank, ]$nn_rank + 1
        }
      }
    }

    ticket <- ticket[order(-unlist(ticket[variable_intensity])), ]
    ticket$intensity_rank <- seq_len(nrow(ticket))

    nn_var <- var(ticket$NNs)
    distance_var <- var(ticket$distance_to_center)
    intensity_var <- var(unlist(ticket[variable_intensity]))
    total_var <- nn_var + distance_var + intensity_var

    nn_var <- 1 + nn_var / total_var
    distance_var <- 1 + distance_var / total_var
    intensity_var <- 1 + intensity_var / total_var

    ticket$intensity_rank <- ticket$intensity_rank * intensity_var
    ticket$distance_rank <- ticket$distance_rank * distance_var
    ticket$nn_rank <- ticket$nn_rank * nn_var
    ticket$cap_rank <- (
      ticket$nn_rank + ticket$distance_rank + ticket$intensity_rank
    ) / nrow(ticket)

    ticket <- ticket[order(ticket$cap_rank), ]
    ticket$core_rank <- seq_len(nrow(ticket))

    if (!exists("core_data")) {
      core_data <- ticket
    } else {
      core_data <- rbind(core_data, ticket)
    }
  }
  write.csv(core_data, paste0("data/capIFA_", file_name), row.names = FALSE)
  cat("Foci compacted.\n")
  foci <<- core_data
}

# Run this function from the `files` directory to collect basic data on previously called
# plaques: image ID, infection status, plaque number, cell count, lowest and highest core
# probabilities, and the plaque's average X/Y coordinate.
plaque_data_collect <- function() {
  file_list <- list.files(path = "data/", pattern = "capIFA")
  if (exists("plaque_info")) {
    rm(plaque_info)
  }

  if (exists("all_cells")) {
    rm(all_cells)
  }

  for (a in file_list) {
    print(paste("Opening file:", a))
    cells <- read.csv(paste0("data/", a))
    cells$top_core <- "Negative"
    for (b in unique(cells$focus)) {
      cells[
        cells$focus == b & cells$core_prob == max(cells[cells$focus == b, ]$core_prob),
      ]$top_core <- "Positive"
    }
    if (!exists("all_cells")) {
      all_cells <- cells
    } else {
      if (FALSE %in% (names(cells) == names(all_cells))) {
        for (i in names(all_cells)[!names(all_cells) %in% names(cells)]) {
          cells[i] <- 0
        }
        for (i in names(cells)[!names(cells) %in% names(all_cells)]) {
          all_cells[i] <- 0
        }
      }
      all_cells <- rbind(all_cells, cells)
    }
    for (c in unique(cells$focus)) {
      tittering <- subset(cells, focus == c)
      if (!exists("plaque_info")) {
        plaque_info <- data.frame(
          "image" = unique(tittering$image),
          "time" = unique(tittering$time),
          "virus" = unique(tittering$virus),
          "infection" = unique(tittering$infection),
          "infected" = !unique(tittering$focus) == 0,
          "plaque_number" = c,
          "plaque_size" = max(tittering$focal_size),
          "percent_G1" = 100 * nrow(subset(
            tittering, ploidy == "2N" & edu == "Negative"
          )) / nrow(tittering),
          "percent_G2" = 100 * nrow(subset(
            tittering, ploidy == "4N" & edu == "Negative"
          )) / nrow(tittering),
          "percent_S" = 100 * nrow(subset(tittering, edu != "Negative")) / nrow(tittering),
          "low_prob" = min(tittering$core_prob),
          "high_prob" = max(tittering$core_prob),
          "Mean_x_position" = mean(tittering$X_ANC_dna),
          "Mean_y_position" = mean(tittering$Y_ANC_dna)
        )
      } else {
        interim <- data.frame(
          "image" = unique(tittering$image),
          "time" = unique(tittering$time),
          "virus" = unique(tittering$virus),
          "infection" = unique(tittering$infection),
          "infected" = !unique(tittering$focus) == 0,
          "plaque_number" = c,
          "plaque_size" = max(tittering$focal_size),
          "percent_G1" = 100 * nrow(subset(
            tittering, ploidy == "2N" & edu == "Negative"
          )) / nrow(tittering),
          "percent_G2" = 100 * nrow(subset(
            tittering, ploidy == "4N" & edu == "Negative"
          )) / nrow(tittering),
          "percent_S" = 100 * nrow(subset(tittering, edu != "Negative")) / nrow(tittering),
          "low_prob" = min(tittering$core_prob),
          "high_prob" = max(tittering$core_prob),
          "Mean_x_position" = mean(tittering$X_ANC_dna),
          "Mean_y_position" = mean(tittering$Y_ANC_dna)
        )
        plaque_info <- rbind(plaque_info, interim)
      }
    }
  }
  write.csv(plaque_info, "data/plaque_data.csv", row.names = FALSE)
  write.csv(all_cells, "data/experiment_data_w_coreProb_all_cells.csv", row.names = FALSE)
}
