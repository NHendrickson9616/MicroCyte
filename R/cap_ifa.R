# Core Analysis of Plaque IFA (CAP-IFA).

foci_caller <- function(df,
                        variable = "n",
                        variable_subset = "Positive",
                        variable_intensity = "Mean_ANC_n",
                        midpointer_x = "XM_ANC_dna",
                        midpointer_y = "YM_ANC_dna",
                        wiggle = 2,
                        overlap_req = 2,
                        output_path = NULL,
                        compact_output_path = NULL) {
  if (!is.data.frame(df)) {
    stop("`df` must be a data frame.", call. = FALSE)
  }
  required_columns <- c("Number", variable, variable_intensity, midpointer_x, midpointer_y)
  missing_columns <- setdiff(required_columns, names(df))
  if (length(missing_columns) > 0L) {
    stop(
      "Missing required column(s): ",
      paste(missing_columns, collapse = ", "),
      call. = FALSE
    )
  }

  cat("Running foci_caller...\n")
  foci_cells <- df[df[[variable]] == variable_subset, , drop = FALSE]

  df[["focus"]] <- 0
  df[["NNs"]] <- 0
  df[["F_ids"]] <- NA_character_

  for (a in seq_len(nrow(foci_cells))) {
    cell_id_a <- foci_cells[["Number"]][[a]]
    x_midpoint <- as.numeric(foci_cells[[midpointer_x]][[a]])
    y_midpoint <- as.numeric(foci_cells[[midpointer_y]][[a]])

    dist_list <- df
    dist_list[["distance"]] <- sqrt(
      (as.numeric(dist_list[[midpointer_x]]) - x_midpoint)^2 +
        (as.numeric(dist_list[[midpointer_y]]) - y_midpoint)^2
    )
    dist_list <- dist_list[order(dist_list[["distance"]]), , drop = FALSE]

    focus <- dist_list[dist_list[["Number"]] == cell_id_a, , drop = FALSE]
    misses <- 0L
    for (b in seq_len(nrow(dist_list))) {
      if (dist_list[["Number"]][[b]] == cell_id_a) {
        next
      }
      if (dist_list[[variable]][[b]] == variable_subset) {
        focus <- rbind(focus, dist_list[b, , drop = FALSE])
      } else {
        misses <- misses + 1L
        if (misses > wiggle) {
          break
        }
      }
    }

    target_row <- df[["Number"]] == cell_id_a
    df[target_row, "NNs"] <- nrow(focus)
    focus_ids <- unique(focus[["Number"]])
    df[target_row, "F_ids"] <- paste0(" ", paste(focus_ids, collapse = " "), " ")
  }

  cat("foci called\n")
  if (!is.null(output_path)) {
    utils::write.csv(df, output_path, row.names = FALSE)
  }

  foci <- foci_compact(
    data = df,
    midpointer_x = midpointer_x,
    midpointer_y = midpointer_y,
    variable = variable,
    variable_subset = variable_subset,
    variable_intensity = variable_intensity,
    overlap = overlap_req,
    output_path = compact_output_path
  )

  list(cells = df, foci = foci, output_path = output_path)
}

foci_compact <- function(data,
                         midpointer_x,
                         midpointer_y,
                         variable = "n",
                         variable_subset = "Positive",
                         variable_intensity = "Mean_ANC_n",
                         overlap = 2,
                         output_path = NULL) {
  if (is.data.frame(data)) {
    cells_new <- data
    source <- NULL
  } else if (is.character(data) && length(data) == 1L) {
    cells_new <- utils::read.csv(data)
    source <- data
  } else {
    stop("`data` must be a data frame or a single CSV file path.", call. = FALSE)
  }

  required_columns <- c(
    "Number", "NNs", "F_ids", variable, variable_intensity,
    midpointer_x, midpointer_y
  )
  missing_columns <- setdiff(required_columns, names(cells_new))
  if (length(missing_columns) > 0L) {
    stop(
      "Missing required column(s): ",
      paste(missing_columns, collapse = ", "),
      call. = FALSE
    )
  }

  if (is.null(source)) {
    cat("Running foci_compact from supplied data.\n")
  } else {
    cat("Running foci_compact from ", source, ".\n", sep = "")
  }

  f_list <- cells_new[cells_new[[variable]] == variable_subset, , drop = FALSE]
  f_list <- f_list[order(f_list[["NNs"]], decreasing = TRUE), , drop = FALSE]
  focal_groups <- list()
  foci_count <- 1L
  cat("Focus groups initialized...\n")

  for (a in seq_len(nrow(f_list))) {
    found <- FALSE
    raw_ids <- trimws(as.character(f_list[["F_ids"]][[a]]))
    f_ids <- if (nzchar(raw_ids)) strsplit(raw_ids, "\\s+")[[1L]] else character()

    if (length(focal_groups) < 1L) {
      focal_groups[[as.character(foci_count)]] <- f_ids
      foci_count <- foci_count + 1L
      found <- TRUE
    } else {
      for (group_name in names(focal_groups)) {
        hits <- sum(f_ids %in% focal_groups[[group_name]])
        if (hits >= overlap) {
          focal_groups[[group_name]] <- unique(c(focal_groups[[group_name]], f_ids))
          found <- TRUE
          break
        }
      }
    }

    if (found == FALSE) {
      focal_groups[[as.character(foci_count)]] <- f_ids
      foci_count <- foci_count + 1L
    }
  }

  cat("Finding overlapping focus groups...\n")
  for (group_name in names(focal_groups)) {
    if (is.null(focal_groups[[group_name]])) {
      next
    }
    found <- FALSE
    group_ids <- focal_groups[[group_name]]
    for (candidate_name in names(focal_groups)) {
      if (group_name != candidate_name && !is.null(focal_groups[[candidate_name]])) {
        hits <- sum(group_ids %in% focal_groups[[candidate_name]])
        if (hits >= overlap && found == FALSE) {
          print(paste(
            "matching", group_name, "with", candidate_name,
            "because they overlap by", hits
          ))
          focal_groups[[candidate_name]] <- unique(c(
            focal_groups[[candidate_name]], group_ids
          ))
          focal_groups[[group_name]] <- NULL
          found <- TRUE
        }
      }
      if (found == TRUE) {
        break
      }
    }
  }

  cat("Determining probabilities...\n")
  cells_new[["focus"]] <- 0L
  cells_new[["focal_size"]] <- 0L
  for (group_name in names(focal_groups)) {
    group_ids <- unique(focal_groups[[group_name]])
    group_rows <- as.character(cells_new[["Number"]]) %in% group_ids
    cells_new[group_rows, "focus"] <- as.integer(group_name)
    cells_new[group_rows, "focal_size"] <- length(group_ids)
  }

  cells_new[["distance_to_center"]] <- NA_real_
  cells_new[["intensity_rank"]] <- 0
  cells_new[["distance_rank"]] <- 0
  cells_new[["nn_rank"]] <- 0

  core_parts <- vector("list", length(unique(cells_new[["focus"]])))
  part_index <- 0L
  for (focus_id in unique(cells_new[["focus"]])) {
    focus_rows <- cells_new[["focus"]] == focus_id
    ticket <- cells_new[focus_rows, , drop = FALSE]

    uni_x <- exp(mean(log(ticket[[midpointer_x]])))
    uni_y <- exp(mean(log(ticket[[midpointer_y]])))
    ticket[["distance_to_center"]] <- sqrt(
      (ticket[[midpointer_x]] - uni_x)^2 +
        (ticket[[midpointer_y]] - uni_y)^2
    )
    ticket <- ticket[order(ticket[["distance_to_center"]]), , drop = FALSE]
    ticket[["distance_rank"]] <- seq_len(nrow(ticket))

    ticket <- ticket[order(-ticket[["NNs"]]), , drop = FALSE]
    ticket[["nn_rank"]] <- seq_len(nrow(ticket))
    if (max(ticket[["nn_rank"]]) > 1L && focus_id != 0L) {
      for (row_rank in seq_len(nrow(ticket) - 1L)) {
        if (ticket[["NNs"]][[row_rank]] == ticket[["NNs"]][[row_rank + 1L]]) {
          ticket[["nn_rank"]][[row_rank + 1L]] <- ticket[["nn_rank"]][[row_rank]]
        } else {
          ticket[["nn_rank"]][[row_rank + 1L]] <-
            ticket[["nn_rank"]][[row_rank]] + 1L
        }
      }
    }

    ticket <- ticket[order(-ticket[[variable_intensity]]), , drop = FALSE]
    ticket[["intensity_rank"]] <- seq_len(nrow(ticket))

    nn_var <- stats::var(ticket[["NNs"]])
    distance_var <- stats::var(ticket[["distance_to_center"]])
    intensity_var <- stats::var(ticket[[variable_intensity]])
    total_var <- nn_var + distance_var + intensity_var

    nn_weight <- 1 + nn_var / total_var
    distance_weight <- 1 + distance_var / total_var
    intensity_weight <- 1 + intensity_var / total_var

    ticket[["intensity_rank"]] <- ticket[["intensity_rank"]] * intensity_weight
    ticket[["distance_rank"]] <- ticket[["distance_rank"]] * distance_weight
    ticket[["nn_rank"]] <- ticket[["nn_rank"]] * nn_weight
    ticket[["cap_rank"]] <- (
      ticket[["nn_rank"]] +
        ticket[["distance_rank"]] +
        ticket[["intensity_rank"]]
    ) / nrow(ticket)

    ticket <- ticket[order(ticket[["cap_rank"]]), , drop = FALSE]
    ticket[["core_rank"]] <- seq_len(nrow(ticket))

    part_index <- part_index + 1L
    core_parts[[part_index]] <- ticket
  }

  core_parts <- core_parts[seq_len(part_index)]
  foci <- do.call(rbind, core_parts)
  rownames(foci) <- NULL
  if (!is.null(output_path)) {
    utils::write.csv(foci, output_path, row.names = FALSE)
  }
  cat("Foci compacted.\n")

  foci
}

plaque_data_collect <- function(data_path,
                                plaque_output_path = NULL,
                                cells_output_path = NULL) {
  if (!is.character(data_path) || length(data_path) != 1L || !dir.exists(data_path)) {
    stop("`data_path` must be an existing directory.", call. = FALSE)
  }

  file_list <- list.files(
    path = data_path,
    pattern = "capIFA",
    full.names = TRUE
  )
  plaque_parts <- list()
  cell_parts <- list()

  for (file_path in file_list) {
    print(paste("Opening file:", file_path))
    cells <- utils::read.csv(file_path)
    cells[["top_core"]] <- "Negative"

    for (focus_id in unique(cells[["focus"]])) {
      focus_rows <- cells[["focus"]] == focus_id
      maximum_core <- max(cells[focus_rows, "core_prob"])
      top_rows <- focus_rows & cells[["core_prob"]] == maximum_core
      cells[top_rows, "top_core"] <- "Positive"
    }
    cell_parts[[length(cell_parts) + 1L]] <- cells

    for (focus_id in unique(cells[["focus"]])) {
      focus_rows <- cells[["focus"]] == focus_id
      focus_cells <- cells[focus_rows, , drop = FALSE]
      is_g1 <- focus_cells[["ploidy"]] == "2N" &
        focus_cells[["edu"]] == "Negative"
      is_g2 <- focus_cells[["ploidy"]] == "4N" &
        focus_cells[["edu"]] == "Negative"
      is_s <- focus_cells[["edu"]] != "Negative"

      plaque_parts[[length(plaque_parts) + 1L]] <- data.frame(
        image = unique(focus_cells[["image"]])[[1L]],
        time = unique(focus_cells[["time"]])[[1L]],
        virus = unique(focus_cells[["virus"]])[[1L]],
        infection = unique(focus_cells[["infection"]])[[1L]],
        infected = unique(focus_cells[["focus"]])[[1L]] != 0,
        plaque_number = focus_id,
        plaque_size = max(focus_cells[["focal_size"]]),
        percent_G1 = 100 * sum(is_g1) / nrow(focus_cells),
        percent_G2 = 100 * sum(is_g2) / nrow(focus_cells),
        percent_S = 100 * sum(is_s) / nrow(focus_cells),
        low_prob = min(focus_cells[["core_prob"]]),
        high_prob = max(focus_cells[["core_prob"]]),
        Mean_x_position = mean(focus_cells[["X_ANC_dna"]]),
        Mean_y_position = mean(focus_cells[["Y_ANC_dna"]]),
        stringsAsFactors = FALSE
      )
    }
  }

  plaque_info <- if (length(plaque_parts) > 0L) {
    do.call(rbind, plaque_parts)
  } else {
    data.frame()
  }

  if (length(cell_parts) > 0L) {
    all_names <- unique(unlist(lapply(cell_parts, names), use.names = FALSE))
    cell_parts <- lapply(cell_parts, function(cells) {
      missing_names <- setdiff(all_names, names(cells))
      for (column_name in missing_names) {
        cells[[column_name]] <- 0
      }
      cells[all_names]
    })
    all_cells <- do.call(rbind, cell_parts)
    rownames(all_cells) <- NULL
  } else {
    all_cells <- data.frame()
  }

  if (!is.null(plaque_output_path)) {
    utils::write.csv(plaque_info, plaque_output_path, row.names = FALSE)
  }
  if (!is.null(cells_output_path)) {
    utils::write.csv(all_cells, cells_output_path, row.names = FALSE)
  }

  list(plaque_info = plaque_info, all_cells = all_cells)
}
