explore <- function(file_name = cells,
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
                    save_file = "data/experiment_explored.csv") {
  # Use `file_name` directly if it is data, or open it if it names a CSV file.
  if (is.data.frame(file_name)) {
    d_frame <- file_name
  } else if (grepl(".csv", file_name)) {
    cat(paste0("Opening file: ", file_name, "...\n"))
    d_frame <- read_csv(file = file_name)
  } else {
    print("File or dataset not recognized.")
  }

  # A new, temporary copy of the name_id is created for seom reason?
  d_frame["placeHolder"] <- d_frame[sort_by]

  # if categories is true, you are asked home many categories to separate the data by
  if (categories == TRUE) {
    cat_num <- as.numeric(readline(prompt = "How many categories would you like to explore: "))
    if (cat_num > 0) {
      # The names to choose from are listed
      print(names(d_frame)[seq_len(length(d_frame) - 1)])
      # A temporary categories table is generated of an appropriate length
      cat_tab <- data.frame("cat_num" = c(1:cat_num))
      # An image_id name is generated for reasons?
      cat_tab$cat_name <- "imageID"

      # For each category number fiven, you're asked which variable should be used
      for (i in 1:cat_num) {
        cat_name <- readline(prompt = paste0("What is the name of category ", i, "? "))
        # Update each row's placeholder variable with its value for the category.
        d_frame$placeHolder <- paste0(d_frame$placeHolder, "_", cat_name, ".")
        d_frame$placeHolder <- paste0(d_frame$placeHolder, unlist(d_frame[cat_name]))
        cat_tab[cat_tab$cat_num == i, ]$cat_name <- cat_name
      }
      # Use the grouped category values rather than a name ID for sorting.
      sort_by <- "placeHolder"
    }
  } else {
    cat_num <- 0
  }

  # Store frequency data using the unique values of the placeholder variable.
  datum <- data.frame("name_id" = unique(unlist(d_frame[sort_by])))

  # We read the schema to get metadata so we don't have to get it later
  schema <- read.csv(scheme)
  # store the meta names
  meta_names <- names(schema)[7:(length(names(schema)) - 1)]
  # Add the metadata to the final dataset row by row
  for (i in meta_names) {
    datum[i] <- "Holder"
    for (j in seq_len(nrow(datum))) {
      target <- unique(unlist(d_frame[d_frame$placeHolder == datum["name_id"][j, ], ][i]))
      datum[i][j, ] <- target
    }
  }

  # We create another id that is sum of all the metadata so we can pull total data later
  if (categories == TRUE && cat_num > 0) {
    for (noi in names(datum)[2:length(names(datum))]) {
      if (!exists("check")) {
        check <- unlist(datum[noi])
      } else {
        check <- paste(check, unlist(datum[noi]), sep = "_")
      }
    }
    datum$totalName <- check
  }

  # Add each population's unique value for every category.
  if (categories == TRUE && cat_num > 0) {
    for (b in cat_tab$cat_name) {
      datum[b] <- "Holder"
      for (c in seq_len(nrow(datum))) {
        datum[b][c, ] <- strsplit(datum["name_id"][c, ], paste0(b, "."))[[1]][2]
        datum[b][c, ] <- strsplit(datum[b][c, ], "_")[[1]][1]
      }
    }
  }

  # Now we pull the total values for each
  datum$total <- 0
  for (totalSet in unique(datum$name_id)) {
    datum[datum$name_id == totalSet, ]$total <- nrow(d_frame[d_frame["placeHolder"] == totalSet, ])

    if (icellate != FALSE && is.numeric(icellate)) {
      source("bin/scripts/analysis/icellate.R")
      if (icellate >= nrow(d_frame[d_frame["placeHolder"] == totalSet, ])) {
        icellate(
          target_cells = d_frame[d_frame["placeHolder"] == totalSet, ],
          folder_name = totalSet,
          verify_size = v_size,
          heat_intensity = icellate_heat_map,
          sampling_number = nrow(d_frame[d_frame["placeHolder"] == totalSet, ]),
          randomize = random,
          line_analyses = line_analyses
        )
      } else {
        icellate(
          target_cells = d_frame[d_frame["placeHolder"] == totalSet, ],
          folder_name = totalSet,
          verify_size = v_size,
          heat_intensity = icellate_heat_map,
          sampling_number = icellate,
          randomize = random,
          line_analyses = line_analyses
        )
      }
    }
  }
  if (categories == TRUE && cat_num > 0) {
    datum$subsetTotal <- datum$total
    datum$total <- 0
    for (totalSet in unique(datum$totalName)) {
      total_count <- sum(subset(datum, totalName == totalSet)$subsetTotal)
      datum[datum$totalName == totalSet, ]$total <- total_count
    }
    datum$subsetPercent <- round(100 * datum$subsetTotal / datum$total, 2)
  }


  if (cell_cycle == TRUE) {
    datum$MDI <- NA
    datum$sPhase <- NA
    datum$G1 <- NA
    datum$G2 <- NA
    datum$re_sPhase <- NA
    datum$G3 <- NA
    datum$S_MFI <- NA
    datum$S_MDI <- NA
    datum$rS_MFI <- NA
    datum$rS_MDI <- NA
    for (a in unique(datum$name_id)) {
      set <- d_frame[d_frame["placeHolder"] == a, ]
      datum[datum$name_id == a, ]$MDI <- mean(set$dna_norm)
      datum[datum$name_id == a, ]$sPhase <- 100 * nrow(set[set$edu == "Positive", ]) / nrow(set)
      datum[datum$name_id == a, ]$G1 <- 100 * nrow(
        set[set$edu == "Negative" & set$ploidy == "2N", ]
      ) / nrow(set)
      datum[datum$name_id == a, ]$G2 <- 100 * nrow(
        set[set$edu == "Negative" & set$ploidy == "4N", ]
      ) / nrow(set)
      datum[datum$name_id == a, ]$re_sPhase <- 100 * nrow(
        set[set$edu == "Positive" & set$ploidy != "2N", ]
      ) / nrow(set)
      datum[datum$name_id == a, ]$G3 <- 100 * nrow(
        set[set$edu == "Negative" & set$ploidy != "2N", ]
      ) / nrow(set)
      datum[datum$name_id == a, ]$S_MFI <- mean(
        set[set$edu == "Positive" & set$ploidy == "2N", ]$edu_norm
      )
      datum[datum$name_id == a, ]$S_MDI <- mean(
        set[set$edu == "Positive" & set$ploidy == "2N", ]$dna_norm
      )
      datum[datum$name_id == a, ]$rS_MFI <- mean(
        set[set$edu == "Positive" & set$ploidy != "2N", ]$edu_norm
      )
      datum[datum$name_id == a, ]$rS_MDI <- mean(
        set[set$edu == "Positive" & set$ploidy != "2N", ]$dna_norm
      )
    }
  }

  if (extra_metrics != FALSE) {
    metrics_number <- as.numeric(readline(
      prompt = "How many extra variables would you like to add: "
    ))
    metrics_tab <- data.frame("metricNumber" = c(1:metrics_number), "metricName" = "holder")
    for (i in seq_len(nrow(metrics_tab))) {
      metrics_tab[i, ]$metricName <- readline(
        prompt = paste0("Which metric should be added (", i, "): ")
      )
    }
    for (metric in unique(metrics_tab$metricName)) {
      datum[paste0("Mean_", metric)] <- NA
      datum[paste0("SD_", metric)] <- NA
      datum[paste0("GMean_", metric)] <- NA
      datum[paste0("Median_", metric)] <- NA
      for (a in unique(datum$name_id)) {
        set <- d_frame[d_frame["placeHolder"] == a, ]
        datum[datum$name_id == a, ][paste0("Mean_", metric)] <- mean(
          as.numeric(unlist(set[metric]))
        )
        datum[datum$name_id == a, ][paste0("SD_", metric)] <- sd(
          as.numeric(unlist(set[metric]))
        )
        datum[datum$name_id == a, ][paste0("GMean_", metric)] <- exp(mean(log(
          as.numeric(unlist(set[metric]))
        )))
        datum[datum$name_id == a, ][paste0("Median_", metric)] <- median(
          as.numeric(unlist(set[metric]))
        )
      }
    }
  }
  if (file.exists(save_file)) {
    concat_away <- readline(prompt = "Save file with same name detected. Concatenate? (Y/n): ")
    if (concat_away != "n") {
      datum_total <- read.csv(save_file)
      datum <- rbind(datum_total, datum)
    } else {
      save_file <- readline(prompt = "What should the name of the new file be: ")
    }
  }
  write.csv(datum, file = save_file, row.names = FALSE)
}
