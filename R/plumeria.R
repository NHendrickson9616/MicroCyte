plumeria <- function(df,
                     xname = "dna_norm",
                     yname = "edu_norm",
                     color_by = c("Mean_NUC_tag", "Mean_NUC_pmcm2"),
                     norm_by = "zscore",
                     bin_size = 10,
                     groups = "name_id",
                     normalize_within_group = FALSE) {
  if (!is.data.frame(df)) {
    stop("`df` must be a data frame.", call. = FALSE)
  }
  required_columns <- unique(c(xname, yname, color_by))
  missing_columns <- setdiff(required_columns, names(df))
  if (length(missing_columns) > 0L) {
    stop(
      "Missing required columns: ",
      paste(missing_columns, collapse = ", "),
      call. = FALSE
    )
  }
  if (length(bin_size) != 1L || !is.numeric(bin_size) || bin_size < 1) {
    stop("`bin_size` must be a positive number.", call. = FALSE)
  }
  if (!norm_by %in% c("zscore", "mean")) {
    stop("`norm_by` must be either 'zscore' or 'mean'.", call. = FALSE)
  }

  if (identical(groups, FALSE) || length(groups) == 0L) {
    df[["target"]] <- "None"
  } else {
    missing_groups <- setdiff(groups, names(df))
    if (length(missing_groups) > 0L) {
      warning(
        "Ignoring missing grouping columns: ",
        paste(missing_groups, collapse = ", "),
        call. = FALSE
      )
    }
    valid_groups <- intersect(groups, names(df))
    if (length(valid_groups) == 0L) {
      df[["target"]] <- "None"
    } else {
      group_values <- lapply(df[valid_groups], as.character)
      df[["target"]] <- do.call(paste, c(group_values, sep = "_"))
    }
  }

  group_indices <- split(seq_len(nrow(df)), df[["target"]], drop = TRUE)
  neighbour_count <- as.integer(bin_size)

  for (color_name in color_by) {
    feather_name <- paste0("feather_", color_name)
    plume_name <- paste0("plume_", color_name)
    df[[feather_name]] <- NA_real_
    df[[plume_name]] <- NA_real_

    normalization_groups <- if (isTRUE(normalize_within_group)) {
      group_indices
    } else {
      list(all_rows = seq_len(nrow(df)))
    }

    for (row_indices in normalization_groups) {
      values <- as.numeric(df[[color_name]][row_indices])
      if (norm_by == "zscore") {
        spread <- stats::sd(values, na.rm = TRUE)
        normalized <- if (is.finite(spread) && spread > 0) {
          (values - mean(values, na.rm = TRUE)) / spread
        } else {
          rep(NA_real_, length(values))
        }
      } else {
        center <- mean(values, na.rm = TRUE)
        normalized <- if (is.finite(center) && center != 0) {
          values / center
        } else {
          rep(NA_real_, length(values))
        }
      }
      df[[feather_name]][row_indices] <- normalized
    }

    for (row_indices in group_indices) {
      x_values <- as.numeric(df[[xname]][row_indices])
      y_values <- as.numeric(df[[yname]][row_indices])
      feather_values <- df[[feather_name]][row_indices]

      for (local_index in seq_along(row_indices)) {
        distances <- sqrt(
          (x_values[[local_index]] - x_values)^2 +
            (y_values[[local_index]] - y_values)^2
        )
        nearest <- order(distances, na.last = NA)
        nearest <- nearest[seq_len(min(neighbour_count, length(nearest)))]
        df[[plume_name]][row_indices[[local_index]]] <- mean(
          feather_values[nearest],
          na.rm = TRUE
        )
      }
    }

    df[[feather_name]] <- NULL
  }

  df
}
