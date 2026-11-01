explore <- function(file_name,
                    sort_by = "name_id",
                    categories = TRUE,
                    cell_cycle = TRUE,
                    scheme = "schema.csv",
                    icellate = FALSE,
                    icellate_heat_map = TRUE,
                    v_size = FALSE,
                    random = TRUE,
                    extra_metrics = FALSE,
                    line_analyses = FALSE,
                    save_file = "data/experiment_explored.csv",
                    category_names = NULL,
                    metric_names = NULL,
                    icellate_fun = NULL) {
  if (is.data.frame(file_name)) {
    d_frame <- file_name
  } else if (is.character(file_name) && length(file_name) == 1L && file.exists(file_name)) {
    d_frame <- as.data.frame(readr::read_csv(file_name, show_col_types = FALSE))
  } else {
    stop("`file_name` must be a data frame or an existing data file.", call. = FALSE)
  }
  if (!sort_by %in% names(d_frame)) {
    stop("`sort_by` was not found in the input data.", call. = FALSE)
  }
  if (!file.exists(scheme)) {
    stop("Schema file does not exist: ", scheme, call. = FALSE)
  }

  d_frame[["placeHolder"]] <- as.character(d_frame[[sort_by]])

  if (isTRUE(categories)) {
    if (is.null(category_names)) {
      category_count <- as.integer(readline(
        prompt = "How many categories would you like to explore: "
      ))
      if (is.na(category_count) || category_count < 0L) {
        stop("The category count must be a non-negative integer.", call. = FALSE)
      }
      category_names <- character(category_count)
      for (index in seq_len(category_count)) {
        category_names[[index]] <- readline(
          prompt = paste0("What is the name of category ", index, "? ")
        )
      }
    }
  } else {
    category_names <- character()
  }

  missing_categories <- setdiff(category_names, names(d_frame))
  if (length(missing_categories) > 0L) {
    stop(
      "Missing category columns: ",
      paste(missing_categories, collapse = ", "),
      call. = FALSE
    )
  }
  for (category_name in category_names) {
    d_frame[["placeHolder"]] <- paste0(
      d_frame[["placeHolder"]],
      "_",
      category_name,
      ".",
      d_frame[[category_name]]
    )
  }

  group_ids <- unique(d_frame[["placeHolder"]])
  datum <- data.frame(name_id = group_ids, stringsAsFactors = FALSE)
  schema_data <- utils::read.csv(scheme, check.names = FALSE)
  schema_names <- names(schema_data)
  metadata_names <- if (length(schema_names) >= 8L) {
    schema_names[seq.int(7L, length(schema_names) - 1L)]
  } else {
    character()
  }
  metadata_names <- intersect(metadata_names, names(d_frame))

  for (metadata_name in metadata_names) {
    datum[[metadata_name]] <- vapply(
      group_ids,
      function(group_id) {
        values <- unique(d_frame[[metadata_name]][d_frame[["placeHolder"]] == group_id])
        if (length(values) == 0L) NA_character_ else as.character(values[[1L]])
      },
      character(1)
    )
  }

  if (length(metadata_names) > 0L) {
    datum[["totalName"]] <- do.call(
      paste,
      c(lapply(datum[metadata_names], as.character), sep = "_")
    )
  }

  for (category_name in category_names) {
    datum[[category_name]] <- vapply(
      group_ids,
      function(group_id) {
        rows <- d_frame[["placeHolder"]] == group_id
        as.character(d_frame[[category_name]][which(rows)[[1L]]])
      },
      character(1)
    )
  }

  group_rows <- split(seq_len(nrow(d_frame)), d_frame[["placeHolder"]], drop = TRUE)
  datum[["total"]] <- lengths(group_rows)[match(datum[["name_id"]], names(group_rows))]

  use_icellate <- is.numeric(icellate) && length(icellate) == 1L &&
    !is.na(icellate) && icellate > 0
  if (use_icellate && !is.function(icellate_fun)) {
    stop("Supply `icellate_fun` when `icellate` is a positive number.", call. = FALSE)
  }
  if (use_icellate) {
    for (group_id in names(group_rows)) {
      rows <- group_rows[[group_id]]
      target_cells <- d_frame[rows, , drop = FALSE]
      icellate_fun(
        target_cells = target_cells,
        folder_name = group_id,
        verify_size = v_size,
        heat_intensity = icellate_heat_map,
        sampling_number = min(as.integer(icellate), nrow(target_cells)),
        randomize = random,
        line_analyses = line_analyses
      )
    }
  }

  if (length(category_names) > 0L) {
    datum[["subsetTotal"]] <- datum[["total"]]
    if ("totalName" %in% names(datum)) {
      totals <- stats::ave(datum[["subsetTotal"]], datum[["totalName"]], FUN = sum)
    } else {
      totals <- rep(sum(datum[["subsetTotal"]]), nrow(datum))
    }
    datum[["total"]] <- totals
    datum[["subsetPercent"]] <- round(100 * datum[["subsetTotal"]] / totals, 2)
  }

  if (isTRUE(cell_cycle)) {
    required_cycle_columns <- c("dna_norm", "edu", "ploidy", "edu_norm")
    missing_cycle_columns <- setdiff(required_cycle_columns, names(d_frame))
    if (length(missing_cycle_columns) > 0L) {
      stop(
        "Cell-cycle metrics require columns: ",
        paste(missing_cycle_columns, collapse = ", "),
        call. = FALSE
      )
    }

    cycle_names <- c(
      "MDI", "sPhase", "G1", "G2", "re_sPhase", "G3",
      "S_MFI", "S_MDI", "rS_MFI", "rS_MDI"
    )
    for (cycle_name in cycle_names) {
      datum[[cycle_name]] <- NA_real_
    }

    safe_mean <- function(values) {
      if (length(values) == 0L || all(is.na(values))) NA_real_ else mean(values, na.rm = TRUE)
    }
    for (row_index in seq_len(nrow(datum))) {
      rows <- group_rows[[datum[["name_id"]][[row_index]]]]
      set <- d_frame[rows, , drop = FALSE]
      count <- nrow(set)
      positive <- set[["edu"]] == "Positive"
      negative <- set[["edu"]] == "Negative"
      two_n <- set[["ploidy"]] == "2N"
      four_n <- set[["ploidy"]] == "4N"

      datum[["MDI"]][row_index] <- safe_mean(set[["dna_norm"]])
      datum[["sPhase"]][row_index] <- 100 * sum(positive, na.rm = TRUE) / count
      datum[["G1"]][row_index] <- 100 * sum(negative & two_n, na.rm = TRUE) / count
      datum[["G2"]][row_index] <- 100 * sum(negative & four_n, na.rm = TRUE) / count
      datum[["re_sPhase"]][row_index] <- 100 * sum(positive & !two_n, na.rm = TRUE) / count
      datum[["G3"]][row_index] <- 100 * sum(negative & !two_n, na.rm = TRUE) / count
      datum[["S_MFI"]][row_index] <- safe_mean(set[["edu_norm"]][positive & two_n])
      datum[["S_MDI"]][row_index] <- safe_mean(set[["dna_norm"]][positive & two_n])
      datum[["rS_MFI"]][row_index] <- safe_mean(set[["edu_norm"]][positive & !two_n])
      datum[["rS_MDI"]][row_index] <- safe_mean(set[["dna_norm"]][positive & !two_n])
    }
  }

  if (!identical(extra_metrics, FALSE)) {
    if (is.null(metric_names)) {
      metric_count <- as.integer(readline(
        prompt = "How many extra variables would you like to add: "
      ))
      if (is.na(metric_count) || metric_count < 0L) {
        stop("The metric count must be a non-negative integer.", call. = FALSE)
      }
      metric_names <- character(metric_count)
      for (index in seq_len(metric_count)) {
        metric_names[[index]] <- readline(
          prompt = paste0("Which metric should be added (", index, "): ")
        )
      }
    }
    missing_metrics <- setdiff(metric_names, names(d_frame))
    if (length(missing_metrics) > 0L) {
      stop("Missing metric columns: ", paste(missing_metrics, collapse = ", "), call. = FALSE)
    }

    for (metric_name in unique(metric_names)) {
      values_by_group <- lapply(group_rows[datum[["name_id"]]], function(rows) {
        as.numeric(d_frame[[metric_name]][rows])
      })
      datum[[paste0("Mean_", metric_name)]] <- vapply(
        values_by_group, mean, numeric(1), na.rm = TRUE
      )
      datum[[paste0("SD_", metric_name)]] <- vapply(
        values_by_group, stats::sd, numeric(1), na.rm = TRUE
      )
      datum[[paste0("GMean_", metric_name)]] <- vapply(
        values_by_group,
        function(values) {
          positive_values <- values[is.finite(values) & values > 0]
          if (length(positive_values) == 0L) NA_real_ else exp(mean(log(positive_values)))
        },
        numeric(1)
      )
      datum[[paste0("Median_", metric_name)]] <- vapply(
        values_by_group, stats::median, numeric(1), na.rm = TRUE
      )
    }
  }

  if (!identical(save_file, FALSE) && is.character(save_file) && length(save_file) == 1L) {
    if (file.exists(save_file)) {
      concatenate <- readline(
        prompt = "Save file with same name detected. Concatenate? (Y/n): "
      )
      if (!identical(tolower(concatenate), "n")) {
        datum <- rbind(utils::read.csv(save_file, check.names = FALSE), datum)
      } else {
        save_file <- readline(prompt = "What should the name of the new file be: ")
      }
    }
    dir.create(dirname(save_file), recursive = TRUE, showWarnings = FALSE)
    utils::write.csv(datum, file = save_file, row.names = FALSE)
  }

  datum
}
