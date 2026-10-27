library(ggpubr)
library(viridis)

smoke_stack <- function(df = cells,
                        norm_by = "zscore",
                        set_x = FALSE,
                        set_y = FALSE,
                        bins = 20,
                        bin_size = FALSE,
                        min_bin_number = 100,
                        lims = FALSE,
                        facet = FALSE,
                        facet_within = TRUE,
                        fixed_coord = TRUE,
                        mid_set = FALSE,
                        save_file = "figures/coloredBar.svg",
                        save_data = FALSE,
                        xdim = 5,
                        ydim = 0.5,
                        viridis_pallet = "inferno",
                        call_icellate = FALSE,
                        fuse_icellate = TRUE,
                        number_icellate = 10) {
  # Asks for a variable to bin on if one isn't given
  if (set_x == FALSE) {
    print(names(df))
    x_name <- readline(prompt = "Which variable should be mapped on the x-axis? ")
  } else {
    x_name <- set_x
  }

  # Makes sure the number of bins is a numeric value
  if (bin_size != FALSE) {
    if (is.numeric(bin_size)) {
      print(paste("Using bin size at", bin_size))
    } else {
      print("Non-numeric bin size detected. Ignoring bin size")
      bin_size <- FALSE
    }
  }

  # Creates the dataframe to store the bin variables
  magic_set <- data.frame("Bin" = c(1:bins))

  # Asks for categorical variables to map on the y-axis if a set isn't supplied
  if (set_y == FALSE) {
    var_numb <- as.numeric(readline(
      prompt = paste0("How many variables should be mapped across ", x_name, "? ")
    ))
    var_mat <- data.frame("var_numb" = c(1:var_numb), "varName" = "HOLD")

    for (i in 1:var_numb) {
      magic_set[, 1 + i] <- readline(prompt = paste0("What is variable number ", i, "? "))
    }
  } else {
    for (i in seq_along(set_y)) {
      magic_set[, 1 + i] <- set_y[i]
    }
  }
  for (i in 2:length(names(magic_set))) {
    names(magic_set)[i] <- unique(magic_set[, i])
    magic_set[, i] <- 0
  }

  # Now we bin the original data set based on the supplied x-axis variable
  # First create a fungible df
  interim <- df
  # Set its bin variable
  interim$bin <- 0
  interim["Name"] <- df[, x_name]
  # Find the bin range and initialize it with the binning variable's global minimum.
  x_name_bin_range <- (max(df[x_name]) - min(df[x_name])) / bins
  bin_set <- c(min(df[x_name]))

  # If a lim double is supplied, the data is trimmed accordingly
  if (is.double(lims)) {
    interim <- interim[interim$Name >= lims[1] & interim$Name <= lims[2], ]
  }

  # Now find the minimum bin value for each bin
  for (i in 1:bins - 1) {
    bin_set <- append(bin_set, rev(bin_set)[1] + x_name_bin_range)
  }
  # Create the bin variable
  bin_id <- 0
  # Iterate through the ranges to bin cells into their respective bins
  for (i in 2:length(bin_set)) {
    bin_id <- bin_id + 1
    if (nrow(interim[interim$Name >= bin_set[i - 1] & interim$Name <= bin_set[i], ]) > 0) {
      interim[interim$Name >= bin_set[i - 1] & interim$Name <= bin_set[i], ]$bin <- bin_id
    }
  }

  if (facet != FALSE) {
    for (silliness in unique(unlist(interim[facet]))) {
      interim_magic <- magic_set

      new_interim <- interim[interim[facet] == silliness, ]

      for (i in names(magic_set)[2:length(magic_set)]) {
        if (norm_by == "zscore") {
          if (facet_within) {
            new_interim$var <-
              (unlist(new_interim[i]) - mean(unlist(interim[i]))) /
              sd(unlist(interim[i]))
          } else {
            new_interim$var <-
              (unlist(new_interim[i]) - mean(unlist(new_interim[i]))) /
              sd(unlist(new_interim[i]))
          }
          legend_name <- "Z-Score"
        } else if (norm_by == "lq") {
          legend_name <- "Lower quartile"
        } else if (norm_by == "uq") {
          legend_name <- "Upper quartile"
        } else if (norm_by == "mean") {
          if (facet_within) {
            new_interim$var <- 100 * unlist(new_interim[i]) / mean(unlist(interim[i]))
          } else {
            new_interim$var <- 100 * unlist(new_interim[i]) / mean(unlist(new_interim[i]))
          }
          legend_name <- "% of Mean"
        } else if (norm_by == "iqr") {
          legend_name <- "Interquartile range"
        } else {
          legend_name <- "Intensity"
        }
        for (j in seq_len(nrow(interim_magic))) {
          if (norm_by == "zscore") {
            interim_magic[i][j, ] <- median(new_interim[new_interim$bin == j, ]$var)
            legend_name <- "Z-Score"
          } else if (norm_by == "lq") {
            interim_magic[i][j, ] <- quantile(unlist(new_interim[new_interim$bin == j, ][i]))[2]
          } else if (norm_by == "uq") {
            legend_name <- "Upper quartile"
            interim_magic[i][j, ] <- quantile(unlist(new_interim[new_interim$bin == j, ][i]))[4]
          } else if (norm_by == "mean") {
            interim_magic[i][j, ] <- median(new_interim[new_interim$bin == j, ]$var)
          } else if (norm_by == "iqr") {
            legend_name <- "Interquartile range"
            interim_magic[i][j, ] <-
              quantile(unlist(new_interim[new_interim$bin == j, ][i]))[4] -
              quantile(unlist(new_interim[new_interim$bin == j, ][i]))[2]
          } else {
            interim_magic[i][j, ] <- median(unlist(new_interim[new_interim$bin == j, ][i]))
            legend_name <- "Intensity"
          }
          if (nrow(new_interim[new_interim$bin == j, ]) < min_bin_number) {
            interim_magic[i][j, ] <- NA
          }
          if (call_icellate && nrow(new_interim[new_interim$bin == j, ]) > 0) {
            if (nrow(new_interim[new_interim$bin == j, ]) < number_icellate) {
              call_it <- nrow(new_interim[new_interim$bin == j, ])
            } else {
              call_it <- number_icellate
            }
            icellate(
              target_cells = new_interim[new_interim$bin == j, ],
              folder_name = paste0(facet, "-", silliness, "_Bin_", j),
              fuse = fuse_icellate,
              sampling_number = call_it,
              line_analyses = FALSE,
              heat_intensity = TRUE
            )
          }
        }
      }
      if (!exists("dumbness")) {
        dumbness <- interim_magic
        dumbness[facet] <- silliness
      } else {
        interim_magic[facet] <- silliness
        dumbness <- rbind(dumbness, interim_magic)
      }
    }
    magic_set <- dumbness
  } else {
    for (i in names(magic_set)[2:length(magic_set)]) {
      if (norm_by == "zscore") {
        interim$var <- (unlist(interim[i]) - mean(unlist(interim[i]))) / sd(unlist(interim[i]))
        legend_name <- "Z-Score"
      } else if (norm_by == "lq") {
        legend_name <- "Lower quartile"
      } else if (norm_by == "uq") {
        legend_name <- "Upper quartile"
      } else if (norm_by == "mean") {
        interim$var <- 100 * unlist(interim[i]) / mean(unlist(interim[i]))
        legend_name <- "% of Mean"
      } else if (norm_by == "iqr") {
        legend_name <- "Interquartile range"
      } else {
        legend_name <- "Intensity"
      }
      for (j in seq_len(nrow(magic_set))) {
        if (norm_by == "zscore") {
          magic_set[i][j, ] <- median(interim[interim$bin == j, ]$var)
          legend_name <- "Z-Score"
        } else if (norm_by == "lq") {
          magic_set[i][j, ] <- quantile(unlist(interim[interim$bin == j, ][i]))[2]
        } else if (norm_by == "uq") {
          legend_name <- "Upper quartile"
          magic_set[i][j, ] <- quantile(unlist(interim[interim$bin == j, ][i]))[4]
        } else if (norm_by == "mean") {
          magic_set[i][j, ] <- median(interim[interim$bin == j, ]$var)
        } else if (norm_by == "iqr") {
          legend_name <- "Interquartile range"
          magic_set[i][j, ] <-
            quantile(unlist(interim[interim$bin == j, ][i]))[4] -
            quantile(unlist(interim[interim$bin == j, ][i]))[2]
        } else {
          magic_set[i][j, ] <- median(unlist(interim[interim$bin == j, ][i]))
          legend_name <- "Intensity"
        }
        if (nrow(interim[interim$bin == j, ]) < min_bin_number) {
          magic_set[i][j, ] <- NA
        }

        if (call_icellate && nrow(interim[interim$bin == j, ]) > 0) {
          if (nrow(interim[interim$bin == j, ]) < number_icellate) {
            call_it <- nrow(interim[interim$bin == j, ])
          } else {
            call_it <- number_icellate
          }
          icellate(
            target_cells = interim[interim$bin == j, ],
            folder_name = paste0("Bin_", j),
            fuse = fuse_icellate,
            sampling_number = call_it,
            line_analyses = FALSE,
            heat_intensity = TRUE,
            save_data = "icellates/bin_data.csv"
          )
        }
      }
    }
  }
  if (facet != FALSE) {
    magic_set <- magic_set |>
      pivot_longer(
        cols = c(2:(ncol(magic_set) - 1)),
        names_to = "Target",
        values_to = "Value"
      )
    magic_set$Target <- paste0(facet, ": ", unlist(magic_set[facet]), " - ", magic_set$Target)
  } else {
    magic_set <- magic_set |>
      pivot_longer(
        cols = c(2:ncol(magic_set)),
        names_to = "Target",
        values_to = "Value"
      )
  }

  c_bin <- c(1:bins)
  c_bin <- c_bin[c(TRUE, FALSE)]

  bin_set <- bin_set[2:length(bin_set)]
  bin_set <- round(bin_set - (x_name_bin_range / 2), 1)
  bin_set <- bin_set[c(TRUE, FALSE)]

  draft <- ggplot(data = magic_set, aes(x = Bin, y = Target, fill = Value)) +
    geom_tile() +
    theme_classic2() +
    scale_fill_viridis(option = viridis_pallet) +
    xlab(x_name) +
    ylab("") +
    guides(fill = guide_legend(title = legend_name)) +
    theme(legend.position = "bottom") +
    scale_x_continuous(breaks = c_bin, labels = bin_set)
  if (fixed_coord == TRUE) {
    draft <- draft + coord_fixed()
  }
  if (mid_set != FALSE && is.numeric(mid_set)) {
    draft <- draft + scale_fill_gradient2(
      low = "blue",
      midpoint = mid_set,
      mid = "black",
      high = "red"
    )
  }

  print(draft)
  if (save_file != FALSE && is.character(save_file)) {
    ggsave(save_file)
  }
  if (save_data != FALSE && is.character(save_data)) {
    write.csv(magic_set, save_data, row.names = FALSE)
  } else if (call_icellate) {
    write.csv(magic_set, "icellates/bin_data.csv", row.names = FALSE)
  }
}
