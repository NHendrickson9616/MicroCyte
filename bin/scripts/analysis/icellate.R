suppressPackageStartupMessages(library(imager))
suppressPackageStartupMessages(library(magick))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(tidyr))

icellate <- function(target_cells,
                     folder_name = "singleCells",
                     d_multiplier = 0.3,
                     verify_size = FALSE,
                     fuse = FALSE,
                     marked = FALSE,
                     heat_intensity = FALSE,
                     heat_color = "inferno",
                     verify_image = "overlay.png",
                     randomize = TRUE,
                     sampling_number = 5,
                     line_analyses = TRUE,
                     number_of_lines = 2,
                     line_graphs = TRUE,
                     initial_angle = 0,
                     random_sampler_denom = 10) {
  # Set a parameter to allow size verification
  size_verified <- FALSE
  fuse_names <- c()

  # First, we'll create the folder in case it doesn't exist, then go to it
  if (!"icellates" %in% list.files()) {
    dir.create("icellates")
  }
  setwd("icellates")
  # Each set of images likely represents a cluster, so create a folder for each cluster.
  if (!folder_name %in% list.files()) {
    dir.create(folder_name)
  }
  setwd(folder_name)

  # Randomly sample the requested number of cells from the population.
  if (randomize == TRUE) {
    target_cells <- target_cells[sample(nrow(target_cells), sampling_number), ]
  }

  for (tCell in seq_len(nrow(target_cells))) {
    target_cell <- target_cells[tCell, ]
    # Now we get the target folder and the target image of the targetCell
    target_name_id <- unique(target_cell$name_id)
    target_image <- unique(target_cell$image)

    # With this, we can make a special folder for JUST THIS CELL
    cell_id_dir <- paste0(target_name_id, "-", unique(target_cell$Number))
    if (!cell_id_dir %in% list.files()) {
      dir.create(cell_id_dir)
      setwd(cell_id_dir)
      # We can also now get a list of the images and move them into the new cellIdDir
      image_list <- list.files(
        path = paste0("../../../files/", target_name_id, "/", target_image, "/PNGS"),
        pattern = ".png", full.names = TRUE
      )
      for (i in image_list) {
        file.copy(i, "./")
      }

      # Now to interate through the images, and them as a dataframe

      if (size_verified == FALSE && verify_size == TRUE) {
        if (!verify_image %in% list.files(pattern = ".png")) {
          print(
            paste(
              "The indicated file for size verification was not found.",
              "Defaulting to general sizing..."
            )
          )
        } else {
          while (TRUE) {
            interim <- load.image(verify_image)
            interim <- as.data.frame(interim)

            # Calculate the required picture area from the image type.
            area_name <- names(cells)[grepl("^Area_ANC", names(target_cell))]
            cell_length <- d_multiplier * round(sqrt(target_cell[area_name] / pi))
            cell_length <- as.integer(round(cell_length))

            # Identify the cell coordinates and trim the data to the calculated length.
            x_name <- names(cells)[grepl("^X_ANC", names(target_cell))]
            x_pos <- as.integer(round(target_cell[x_name]))
            x_min <- round(x_pos - cell_length)
            x_max <- round(x_pos + cell_length)

            y_name <- names(cells)[grepl("^Y_ANC", names(target_cell))]
            y_pos <- as.integer(round(target_cell[y_name]))
            y_min <- round(y_pos - cell_length)
            y_max <- round(y_pos + cell_length)

            # Now we remove all the pixels we don't care about
            interim <- interim[interim$x >= x_min & interim$x <= x_max, ]
            interim <- interim[interim$y >= y_min & interim$y <= y_max, ]

            # Now to reset the corrdinates
            interim$x <- 1 + interim$x - min(interim$x)
            interim$y <- 1 + interim$y - min(interim$y)

            # Now to convert the dataframe to an image object and save the image
            new_plot <- suppressWarnings(as.cimg(interim))
            plot(new_plot)
            print(paste0("dMultiplier is currently set to: ", d_multiplier))
            check_box <- readline(prompt = "Does this look like the appropriate size (Y/n)? ")
            if (check_box == "n") {
              d_multiplier <- as.numeric(readline(prompt = "What should the new dMultiplier be: "))
            } else {
              size_verified <- TRUE
              break
            }
          }
        }
      }

      # For each image...
      for (i in list.files()) {
        # we load the image..
        interim_i <- load.image(i)

        # And convert it to a dataframe.
        interim <- as.data.frame(interim_i)

        # Calculate the required picture area from the image type.
        area_name <- names(cells)[grepl("^Area_ANC", names(target_cell))]
        cell_length <- d_multiplier * round(sqrt(target_cell[area_name] / pi))
        cell_length <- as.integer(round(cell_length))

        # Identify the cell coordinates and trim the data to the calculated length.
        x_name <- names(cells)[grepl("^X_ANC", names(target_cell))]
        x_pos <- as.integer(round(target_cell[x_name]))
        x_min <- round(x_pos - cell_length)
        x_max <- round(x_pos + cell_length)

        y_name <- names(cells)[grepl("^Y_ANC", names(target_cell))]
        y_pos <- as.integer(round(target_cell[y_name]))
        y_min <- round(y_pos - cell_length)
        y_max <- round(y_pos + cell_length)

        # Now we remove all the pixels we don't care about
        interim <- interim[interim$x >= x_min & interim$x <= x_max, ]
        interim <- interim[interim$y >= y_min & interim$y <= y_max, ]

        # Now we save the new centroid positions for possible line analyses
        new_x_pos <- 1 + x_pos - min(interim$x)
        new_y_pos <- 1 + y_pos - min(interim$y)

        # Now to reset the corrdinates
        interim$x <- 1 + interim$x - min(interim$x)
        interim$y <- 1 + interim$y - min(interim$y)


        # Now to convert the dataframe to an image object and save the image
        new_plot <- suppressWarnings(as.cimg(interim))
        new_plot <- cimg2magick(new_plot, rotate = TRUE)
        image_write(new_plot, path = i, format = "png")

        if (marked) {
          crosshair <- cell_length / 20
          if (crosshair < 1) {
            crosshair <- 1
          }
          on_horizontal <- interim$x > new_x_pos - crosshair &
            interim$x < new_x_pos + crosshair &
            interim$y == new_y_pos
          on_vertical <- interim$y > new_y_pos - crosshair &
            interim$y < new_y_pos + crosshair &
            interim$x == new_x_pos
          interim[on_horizontal | on_vertical, ]$value <- 255
          new_plot <- suppressWarnings(as.cimg(interim))
          new_plot <- cimg2magick(new_plot, rotate = TRUE)
          jimmy <- strsplit(i, ".png")[[1]][1]
          jimmy <- paste0(jimmy, "_marked.png")
          image_write(new_plot, path = jimmy, format = "png")
        }

        if (heat_intensity) {
          old_image <- load.image(i)
          thick <- as.data.frame(old_image)
          check_gray <- as.data.frame(grayscale(old_image))

          check_gray$R <- 0
          check_gray$G <- 0
          check_gray$B <- 0

          fake_plot <- ggplot(data = check_gray, aes(fill = value, x = x, y = y)) +
            geom_point() +
            scale_fill_viridis_c(option = heat_color)
          fake_plot <- ggplot_build(fake_plot)
          fake_plot <- as.data.frame(fake_plot$data)
          check_gray$hex <- fake_plot$fill

          for (viro in unique(check_gray$hex)) {
            rgb_code <- col2rgb(viro)
            check_gray[check_gray$hex == viro, ]$R <- rgb_code[1]
            check_gray[check_gray$hex == viro, ]$G <- rgb_code[2]
            check_gray[check_gray$hex == viro, ]$B <- rgb_code[3]
          }

          check_gray$R <- check_gray$R / 255
          check_gray$G <- check_gray$G / 255
          check_gray$B <- check_gray$B / 255

          thick[thick$cc == 1, ]$value <- check_gray$R
          thick[thick$cc == 2, ]$value <- check_gray$G
          thick[thick$cc == 3, ]$value <- check_gray$B

          new_plot <- suppressWarnings(as.cimg(thick))
          new_plot <- cimg2magick(new_plot, rotate = TRUE)
          new_plot <- image_flop(new_plot)
          image_write(new_plot, path = paste0("heatMap_", i), format = "png")
        }

        # If class averaging is called via the 'fuse' parameter AND the image is a square...
        if (fuse == TRUE && max(interim$x) == max(interim$y)) {
          # so we can call all the averaged images later, we'll make a list of each unique image
          if (!strsplit(i, ".png")[[1]][1] %in% fuse_names) {
            fuse_names <- append(fuse_names, strsplit(i, ".png")[[1]][1])
          }
          # Generate a parameter for the data frame unless the variable already exists.
          if (!exists(paste0(strsplit(i, ".png")[[1]][1], "_fused"))) {
            bidet <- interim
            bidet$number <- 1
            assign(paste0(strsplit(i, ".png")[[1]][1], "_fused"), bidet)
            assign(paste0(strsplit(i, ".png")[[1]][1], "_xOri"), new_x_pos)
            assign(paste0(strsplit(i, ".png")[[1]][1], "_yOri"), new_y_pos)
          } else {
            # Otherwise, generate an interim data frame and calculate the average.
            fuse_old <- get(paste0(strsplit(i, ".png")[[1]][1], "_fused"))
            # The number is increased
            fuse_old$number <- fuse_old$number + 1
            fuse_interim <- interim

            # If the images are the same size, then the values are summed
            if (max(fuse_interim$x) == max(get(paste0(strsplit(i, ".png")[[1]][1], "_fused"))$x)) {
              fuse_old$value <- fuse_old$value + fuse_interim$value
              assign(paste0(strsplit(i, ".png")[[1]][1], "_fused"), fuse_old)
            } else {
              # Now we use the midpoints XPos and YPos to align the dataframes

              if (nrow(fuse_old) > nrow(fuse_interim)) {
                x_shift <- get(paste0(strsplit(i, ".png")[[1]][1], "_xOri")) - new_x_pos
                y_shift <- get(paste0(strsplit(i, ".png")[[1]][1], "_yOri")) - new_y_pos
                fuse_interim$x <- fuse_interim$x + x_shift
                fuse_interim$y <- fuse_interim$y + y_shift
              } else {
                x_shift <- new_x_pos - get(paste0(strsplit(i, ".png")[[1]][1], "_xOri"))
                y_shift <- new_y_pos - get(paste0(strsplit(i, ".png")[[1]][1], "_yOri"))
                fuse_old$x <- fuse_old$x + x_shift
                fuse_old$y <- fuse_old$y + y_shift
                assign(paste0(strsplit(i, ".png")[[1]][1], "_xOri"), new_x_pos)
                assign(paste0(strsplit(i, ".png")[[1]][1], "_yOri"), new_y_pos)
              }

              if (nrow(fuse_old[fuse_old$x == 0 | fuse_old$y == 0, ]) >= 1) {
                print(paste0(
                  "A total of",
                  nrow(fuse_old[fuse_old$x == 0 | fuse_old$y == 0, ]),
                  "zero(s) was detected in the old"
                ))
                fuse_old <- fuse_old[fuse_old$x > 0 & fuse_old$y > 0, ]
              }
              if (nrow(fuse_interim[fuse_interim$x == 0 | fuse_interim$y == 0, ]) >= 1) {
                print(paste0(
                  "A total of ",
                  nrow(fuse_interim[fuse_interim$x == 0 | fuse_interim$y == 0, ]),
                  "zero(s) was detected in the interim"
                ))
                fuse_interim <- fuse_interim[fuse_interim$x > 0 & fuse_interim$y > 0, ]
              }
              row.names(fuse_old) <- paste0(fuse_old$x, "_", fuse_old$y, "_", fuse_old$cc)
              row.names(fuse_interim) <- paste0(
                fuse_interim$x, "_", fuse_interim$y, "_", fuse_interim$cc
              )

              fuse_new <- merge(fuse_old[, 1:4], fuse_interim, by = 0, all = TRUE)
              fuse_new[is.na(fuse_new)] <- 0
              the_numbers <- as.numeric(unique(fuse_old$number))

              if (nrow(fuse_old) > nrow(fuse_interim)) {
                fuse_export <- data.frame(
                  "x" = fuse_new$x.x,
                  "y" = fuse_new$y.x,
                  "cc" = fuse_new$cc.x,
                  "value" = fuse_new$value.y + fuse_new$value.x,
                  "number" = the_numbers
                )
              } else {
                fuse_export <- data.frame(
                  "x" = fuse_new$x.y,
                  "y" = fuse_new$y.y,
                  "cc" = fuse_new$cc.y,
                  "value" = fuse_new$value.y + fuse_new$value.x,
                  "number" = the_numbers
                )
              }
              # Now we create the new, averaged image
              assign(paste0(strsplit(i, ".png")[[1]][1], "_fused"), fuse_export)
            }
          }
        }
      }

      # Now to run the line analyses if they're supposed to be
      if (line_analyses == TRUE) {
        line_analysis(
          mid_x = new_x_pos,
          mid_y = new_y_pos,
          image_name = interim,
          lines_number = number_of_lines,
          line_length = cell_length,
          graphing = line_graphs,
          offset = initial_angle,
          total_corr_pix_fraction = random_sampler_denom
        )
      }
      setwd("../")
    }
  }
  if (fuse == TRUE) {
    fuse_data <- data.frame("image" = fuse_names, "number" = 0)
    # this is where you call the assigned things and save them as images
    for (fused in fuse_names) {
      # for each image name that was found, we get the image
      fused_image <- get(paste0(fused, "_fused"))
      # The intensities get rounded and then we drop the image numbers after saving that data
      fused_image$value <- fused_image$value / fused_image$number
      fuse_data[fuse_data$image == fused, ]$number <- max(fused_image$number)
      fused_image <- fused_image[, 1:4]
      # The averaged image is then saved as a png
      fused_plot <- suppressWarnings(as.cimg(fused_image))
      fused_plot <- cimg2magick(fused_plot, rotate = TRUE)
      image_write(fused_plot, path = paste0(fused, ".png"), format = "png")
    }
    # The number of images used is then saved
    write.csv(fuse_data, "ClassAveragingMetrics.csv", row.names = FALSE)
  }
  setwd("../../")
}

line_analysis <- function(mid_x,
                          mid_y,
                          image_name,
                          lines_number = 2,
                          line_length,
                          total_corr_pix_fraction = 10,
                          graphing = TRUE,
                          offset = 0) {
  # First we make a directory to store the graphs and lines
  if (!"line_analyses" %in% list.files()) {
    dir.create("line_analyses")
  }
  setwd("line_analyses")
  if (!"csv" %in% list.files()) {
    dir.create("csv")
  }
  setwd("../")

  # Now we'll pull out the images images to combine them into a single pixel dataframe
  png_list <- list.files(pattern = ".png")
  png_list <- png_list[!grepl("_marked.png", png_list)]
  for (pngImage in png_list) {
    # We'll ignore the overlay image because we may be using more than RBG
    if (pngImage != "overlay.png" && !grepl("heatMap_", pngImage)) {
      # Load the image as an object
      line_image <- load.image(pngImage)
      # Convert it to grayscale dataframe
      line_image <- as.data.frame(grayscale(line_image))
      line_image$value <- line_image$value - min(line_image$value)
      # Rename the pixel value column as the stain ID
      names(line_image)[3] <- strsplit(pngImage, ".png")[[1]][1]

      # Join the iterated extractions into a single dataset
      if (!exists("line_data")) {
        line_data <- line_image
      } else {
        line_data <- cbind(line_data, line_image[strsplit(pngImage, ".png")[[1]][1]])
      }
      # Then we randomly select a fraction of the pixels to save for whole cell analysis
      random_points <- line_data[
        sample(nrow(line_data), nrow(line_data) / total_corr_pix_fraction),
      ]
      write.csv(random_points, "lineAnalyses/WholeCellRandom.csv", row.names = FALSE)
    }
  }

  # now we pull out lines of pixel data that run through the center of the nucleus

  # If the linesNumber variable is 2, it just takes the x and y midpoint lines,
  # otherwise it rotates radially though the center
  if (lines_number < 3) {
    max_x <- mid_x + (line_length)
    min_x <- mid_x - (line_length)
    x_line <- subset(line_data, y == mid_y & x <= max_x & x >= min_x)
    write.csv(x_line, file = "lineAnalyses/csv/xLine.csv", row.names = FALSE)
    max_y <- mid_y + (line_length)
    min_y <- mid_y - (line_length)
    y_line <- subset(line_data, x == mid_x & y <= max_y & y >= min_y)
    write.csv(y_line, file = "lineAnalyses/csv/yLine.csv", row.names = FALSE)
  } else {
    # we'll establish the initial angle as 0
    line_angle <- 0
    # then we figure out the angle between lines
    angle_number <- round(360 / lines_number)
    # now we iterate through each line to extract the line data
    for (line in 1:lines_number) {
      # We add the angle difference from the previous line and convert it radians
      line_angle <- line_angle + angle_number
      angle_radians <- line_angle * pi / 180
      # Using the length parameter, we'll find the starting x and y pixels for the line
      x_start <- round(sin(angle_radians) * line_length + mid_x)
      y_start <- round(cos(angle_radians) * line_length + mid_y)
      x_end <- round(mid_x - sin(angle_radians) * line_length)
      y_end <- round(mid_y - cos(angle_radians) * line_length)

      # The changes from the midpoint to the end mirror those at the start.

      # Now that we've got a mid-point and a starting point, we can get the slope of the line
      if (mid_x != x_start && mid_y != y_start) {
        line_slope <- -(mid_y - y_start) / (mid_x - x_start)
        y_int <- mid_y + (line_slope * mid_x)
        new_line <- subset(
          line_data,
          x <= (max(x_start, x_end)) &
            x >= (min(x_start, x_end)) &
            y == -round(line_slope * x - y_int)
        )
      } else {
        new_line <- subset(
          line_data,
          x <= (max(x_start, x_end)) &
            x >= (min(x_start, x_end)) &
            y <= (max(y_start, y_end)) &
            y >= (min(y_start, y_end))
        )
      }
      # Use the slope and endpoints to select pixels along the line.
      write.csv(
        new_line,
        file = paste0("lineAnalyses/csv/lineAngle_", line_angle, ".csv"),
        row.names = FALSE
      )
    }
  }

  # Now that we've got the line data, we can graph them because - why not?
  if (graphing == TRUE) {
    setwd("lineAnalyses/csv")
    line_files <- list.files(pattern = ".csv")
    for (lineFile in line_files) {
      line_data <- read.csv(lineFile)
      if (!"pix" %in% names(line_data)) {
        line_data$pix <- seq_len(nrow(line_data))
        col_length <- ncol(line_data) - 1
        line_data <- line_data |> pivot_longer(
          cols = 3:col_length,
          names_to = "color",
          values_to = "intensity"
        )
        line_data <- line_data[, c(3, 1, 2, 4, 5)]
        write.csv(line_data, file = paste0("longer_", lineFile), row.names = FALSE)
      }
      linera <- ggplot(
        data = line_data,
        aes(x = pix, y = intensity, color = color)
      ) +
        geom_step(linewidth = 1.5) +
        xlab("Relative pixel position") +
        ylab("Relative pixel intentsity") +
        theme_classic() +
        theme(
          axis.line = element_line(colour = "black", linewidth = 1.5),
          axis.ticks = element_line(colour = "black", linewidth = 1.5),
          axis.text = element_text(face = "bold", color = "black", size = 12),
          axis.title = element_text(face = "bold", color = "black", size = 20),
          legend.position = "top"
        )
      print(linera)
      pdf_name <- strsplit(lineFile, ".csv")[[1]][1]
      ggsave(filename = paste0("../", pdf_name, ".pdf"), height = 5, width = 5, units = "in")
    }
    setwd("../../")
  }
}
