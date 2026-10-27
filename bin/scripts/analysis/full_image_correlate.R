suppressPackageStartupMessages(library(imager))
suppressPackageStartupMessages(library(magick))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(stringr))

full_image_correlate <- function(cor_images = c("dna", "edu"),
                                 sub_percent = 50,
                                 background_subtract = TRUE,
                                 threshold_intensity = 0.1,
                                 remove_overexposed = TRUE,
                                 graph_dist = TRUE) {
  scheme <- read_csv("schema.csv")

  setwd("files/")
  folder_set <- list.files()
  for (i in folder_set) {
    if (!grepl(".csv", i)) {
      setwd(i)
      print(paste0("Working in file: ", i))
      image_set <- list.files()
      for (j in image_set) {
        if (!grepl(".csv", j)) {
          print(paste0("Working in image: ", j))
          setwd(paste0(j, "/PNGS"))
          for (k in cor_images) {
            if (paste0(k, ".png") %in% list.files()) {
              old_image <- load.image(paste0(k, ".png"))
              thick <- as.data.frame(old_image)
              if (4 %in% unique(thick$cc)) {
                print("Alpha channel detected. Removing and re-saving.")
                thick <- thick[thick$cc != 4, ]
                new_plot <- suppressWarnings(as.cimg(thick))
                new_plot <- cimg2magick(new_plot, rotate = TRUE)
                new_plot <- image_flop(new_plot)
                image_write(new_plot, path = image_location, format = "png")
                old_image <- load.image(image_location)
                thick <- as.data.frame(old_image)
                print("3-Channel PNG saved.")
              }
              if (background_subtract) {
                thick$value <- thick$value - min(thick$value)
              }

              names(thick)[3] <- k
              if (!exists("image_df")) {
                image_df <- thick
              } else {
                image_df <- cbind(image_df, thick[, 3])
                names(image_df)[length(image_df)] <- k
              }
            } else {
              print(paste0(
                "Unable to find ", k, ".png in file: ", i, "/", j,
                ". Skipping this image."
              ))
            }
          }
          if (exists("image_df")) {
            image_df$image <- j

            # Remove pixels below the threshold
            for (riff in cor_images) {
              image_df <- image_df[image_df[riff] > threshold_intensity, ]
              if (remove_overexposed) {
                image_df <- image_df[image_df[riff] < 1, ]
              }
            }

            # Randomly sample pixels
            row_keep <- round((sub_percent / 100) * nrow(image_df), 0)
            if (row_keep > nrow(image_df)) {
              row_keep <- nrow(image_df)
            }
            image_df <- image_df[sample(nrow(image_df), row_keep), ]

            write_csv(image_df, "image_corr.csv")
            if (!exists("all_image_corr")) {
              all_image_corr <- image_df
            } else {
              all_image_corr <- rbind(all_image_corr, image_df)
            }
          }

          rm(image_df)
          setwd("../../")
        }
      }
      if (exists("all_image_corr")) {
        jimmy <- scheme[scheme$name_id == i, ]
        for (meta in 7:ncol(jimmy)) {
          gin <- as.character(unique(jimmy[, meta]))
          all_image_corr[names(jimmy)[meta]] <- gin
        }
        write_csv(all_image_corr, "all_image_corr.csv")

        if (!exists("experiment_corr")) {
          experiment_corr <- all_image_corr
        } else {
          experiment_corr <- rbind(experiment_corr, all_image_corr)
        }
        rm(all_image_corr)
        cat("\n")
      }
      setwd("../")
    }
  }
  setwd("../")
  if (exists("experiment_corr")) {
    write_csv(experiment_corr, "data/experiment_correlate.csv")
    if (graph_dist) {
      for (riff in cor_images) {
        draft <- ggplot(data = experiment_corr, aes(x = unlist(experiment_corr[, riff]))) +
          geom_density() +
          xlab(paste0("Intensity profile of ", riff)) +
          theme_bw()
        print(draft)
        ggsave(paste0("figures/Corr_plot_", riff, ".png"), dpi = 300)
      }
    }
    rm(experiment_corr)
  }
}
