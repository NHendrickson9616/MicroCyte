plumeria <- function(df = cells,
                     xname = "dna_norm",
                     yname = "edu_norm",
                     color_by = c("Mean_NUC_tag", "Mean_NUC_pmcm2"),
                     norm_by = "zscore",
                     bin_size = 10,
                     groups = c("name_id"),
                     normalize_within_group = FALSE) {
  if (groups != FALSE) {
    df$target <- ""
    for (group in groups) {
      if (group %in% names(df)) {
        df$target <- paste0(df$target, unlist(df[group]))
      } else {
        print(paste0("Variable: ", group, " was not found in dataset. Skipping..."))
      }
    }
  } else {
    df$target <- "None"
  }
  for (color in color_by) {
    print(color)
    feather <- paste0("feather_", color)
    plumage <- paste0("plume_", color)
    df[feather] <- 0
    df[plumage] <- 0
    if (normalize_within_group) {
      for (group in unique(df$target)) {
        interim <- df[df$target == group, ]
        if (norm_by == "zscore") {
          df[df$target == group, feather] <- (
            unlist(interim[color]) - mean(unlist(interim[color]))
          ) / sd(unlist(interim[color_by]))
        } else {
          df[df$target == group, feather] <- unlist(interim[color]) / mean(unlist(interim[color]))
        }
      }
    } else {
      if (norm_by == "zscore") {
        df[, feather] <- (unlist(df[color]) - mean(unlist(df[color]))) / sd(unlist(df[color]))
      } else {
        df[, feather] <- unlist(df[color]) / mean(unlist(df[color]))
      }
    }

    for (i in unique(df$target)) {
      interim <- df[df$target == i, ]

      for (n in seq_len(nrow(interim))) {
        interim$distance <- sqrt(
          (as.numeric(interim[n, xname]) - unlist(interim[, xname]))^2 +
            (as.numeric(interim[n, yname]) - unlist(interim[, yname]))^2
        )
        df[df$target == i, ][n, ][plumage] <- mean(unlist(
          interim[order(interim[, "distance"]), ][1:bin_size, ][feather]
        ))

        if (n %% 100 == 0) {
          print(paste0("Group - ", i, ": ", n, " of ", nrow(interim)))
          print(mean(unlist(df[plumage])))
        }
      }
    }
    df <- df[, -which(names(df) == feather)]
  }
  df
}
