# SirMixaPlot version 8.3 - Less is more
#------------------------------------------------------------------------------------------------

# This program takes a CSV file and generates scatterplots based on its columns.

# Changes from version 8.2:
#   -Removed redundant functions, mostly gate() and negPopCorrect, which no longer serve a purpose
#   -Streamlined some functions to remove unnecessary global variables
#   -Changed normalizer to use LMean_ANC_edu instead of the now defunct ALIMean_ANC_edu

# Changes from version 8.1:
#   -Changed the default EdU neg setting on normalizer()
#   -Added a work-around for concater() should rbinding fail due to missing columns
#   -imaGen() now creates a concatenated CSV of all WholeCell ROIs with the _ROI_ tag
#   -Added runPCA() function to take a dataframe and generate a PCA biplot; still in beta testing

# Changes from SirMixaPlot version 8.0:
#   -Added some comments and typos
#   -Added the reMap() function for re-creating a facsimile of the original image

# Changes from SirMixaPlot version 7.0:
#   -Takes data in from the YggData.imj macro
#   -imaGen() is now much bigger and generates combined CSVs for nuclear vs whole cell comparisons

# Things still that need to be coded:

# Program opens required packages
# If these packages aren't install, use the '0_setup.R' function to install them.

suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(ggpubr))
suppressPackageStartupMessages(library(MASS))
suppressPackageStartupMessages(library(viridis))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(RecordLinkage))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(plotly))


#------------------------------------------------------------------------------------------------

# Functions for your pleasure

#----------------------------------------------------------------
# sirmixaplot() loads filo and calls joiner() to graph when run is TRUE.
# This may lead to negPopCorrect if the appropriate names are not present.
# To re-run the negPopCorrect, simply call adjust as TRUE
sirmixaplot <- function(filo,
                        run = TRUE,
                        adjust = TRUE) {
  # The script, m'lord

  # Begin the prompts for desired output conditions and add safety columns.
  thingee <<- filo
  cat(paste("Opening file: ", thingee, sep = ""))
  cat("\n")

  # A '_cells' filename signals prior processing; otherwise, create a cells file.
  if (grepl("_cells.csv", thingee, fixed = TRUE)) {
    cells <<- read_csv(file = thingee)
  } else {
    cells <<- read_csv(file = thingee)
    cat("Creating a new file, sir or madam.")
    cat("\n")
    thingee <<- paste0(substr(thingee, 1, nchar(thingee) - 4), "_cells.csv")
    cat(paste(thingee, "created."))
    if ("IntDen_ANC_dna" %in% names(cells) && "IntDen_ANC_anchor" %in% names(cells)) {
      log2_name <- readline(
        prompt = paste0(
          "Both 'DNA' and 'anchor' found. ",
          "Which should be used for log2_dna? (ANCHOR/dna)? "
        )
      )
      if (log2_name == "dna") {
        cells$log2_dna <<- log(cells$IntDen_ANC_dna, 2)
      } else {
        cells$log2_dna <<- log(cells$IntDen_ANC_anchor, 2)
      }
    } else if ("IntDen_ANC_dna" %in% names(cells)) {
      cells$log2_dna <<- log(cells$IntDen_ANC_dna, 2)
    }

    for (i in seq_len(ncol(cells))) {
      if (grepl("Mean_", names(cells)[i])) {
        cells[paste0("L", names(cells)[i])] <<- log(cells[, i], 10)
      }
    }
    # writes the file
    write_csv(cells, thingee)
  }
  cat("\n")

  # Checks that some size correctiong has occured and begins making the graph
  if (run == TRUE) {
    joiner()
    # writes the file
    write_csv(cells, thingee)
  }
  if (!"dna_norm" %in% names(cells) && TRUE %in% grepl("ANC_edu", names(cells))) {
    normalizer()
    # writes the file
    write_csv(cells, thingee)
  }
}

#--------------------------------------------------------------------------------------
# joiner() connects sirmixaplot() with grapho() and can graph without sirmixaplot().
# The default dataset is cells, but any dataset can be used.
joiner <- function(df = cells) {
  # newAxis is reset for the changeAxis() function
  new_axis <<- FALSE

  # cats the dataset columns
  cat("\n")
  cat("Data categories detected:")
  cat("\n")
  print(sort(names(df)))
  cat("\n")

  # This sets the x and y variables
  px <<- readline(prompt = "px - What parameter is the x-axis: ")
  xname <<- readline(prompt = "What is the name of the x axis: ")
  cat("\n")
  py <<- readline(prompt = "py - What parameter is the y-axis: ")
  yname <<- readline(prompt = "What is the name of the y axis: ")
  cat("Working...")
  cat("\n")

  # calls the graphing function
  grapho(
    df = cells,
    x = px,
    y = py,
    x_name = xname,
    y_name = yname
  )
}

# grapho() actually makes the scatter plots
grapho <- function(df = cells,
                   x = px,
                   y = py,
                   x_name = xname,
                   y_name = yname) {
  # backup dataframe is generated
  interim <<- df

  cat("Data frame generated")
  cat("\n")

  # Initial scatter plot is generated
  qq <<- ggplot(data = df, aes_string(x, y)) +
    theme_linedraw() +
    theme(plot.title = element_text(size = 40, face = "bold")) +
    theme(axis.title = element_text(size = 30, face = "bold")) +
    geom_point(size = 2) +
    scale_x_continuous(name = x_name) +
    scale_y_continuous(name = y_name) +
    theme(axis.text = element_text(face = "bold", size = 16), axis.ticks = element_blank())

  if (!exists("new_axis")) {
    new_axis <<- FALSE
  }
  if (new_axis == TRUE) {
    changeAxis(store[1], store[2], store[3], store[4])
  }

  # Program finishes up
  cat("Plot generated")
  cat("\n")
  print(qq)
}

#----------------------------------------------------------------------------------
# This function creates quadrants and calculates the percentage in each quadrant
#   Geom_line function can be altered to change the aesthetics of the quadrant lines

modthequad <- function(xi = 1.5,
                       yi = 1.5,
                       df = cells,
                       saveit = FALSE,
                       pause = TRUE,
                       graphit = TRUE) {
  if (pause == TRUE) {
    xi <- readline(prompt = "What is the x-intercept: ")
    yi <- readline(prompt = "What is the y-intercept: ")
  }
  quad <- c("Q1", "Q2", "Q3", "Q4")
  todos <- nrow(df)
  qt1 <<- nrow(subset(df, df[px] < xi & df[py] > yi))
  qt2 <<- nrow(subset(df, df[px] > xi & df[py] > yi))
  qt3 <<- nrow(subset(df, df[px] < xi & df[py] < yi))
  qt4 <<- nrow(subset(df, df[px] > xi & df[py] < yi))
  percent_total <- c(
    format(round(qt1 / todos * 100, 2), nsmall = 2),
    format(round(qt2 / todos * 100, 2), nsmall = 2),
    format(round(qt3 / todos * 100, 2), nsmall = 2),
    format(round(qt4 / todos * 100, 2), nsmall = 2)
  )
  quads <- data.frame(quad, percent_total)
  xmin <- min(df[px])
  xmax <- max(df[px])
  ymin <- min(df[py])
  ymax <- max(df[py])
  qqmod <<- qq
  qqmod <<- qqmod +
    geom_hline(yintercept = as.numeric(yi), color = "red", linetype = "dashed", linewidth = 2) +
    geom_vline(xintercept = as.numeric(xi), color = "red", linetype = "dashed", linewidth = 2) +
    annotate(
      "text", x = xmin * 1.05, y = ymax * 0.95,
      label = paste0("Q1: ", format(round(qt1 / todos * 100, 2), nsmall = 2), "%")
    ) +
    annotate(
      "text", x = xmax * 0.95, y = ymax * 0.95,
      label = paste0("Q2: ", format(round(qt2 / todos * 100, 2), nsmall = 2), "%")
    ) +
    annotate(
      "text", x = xmin * 1.05, y = ymin * 1.05,
      label = paste0("Q3: ", format(round(qt3 / todos * 100, 2), nsmall = 2), "%")
    ) +
    annotate(
      "text", x = xmax * 0.95, y = ymin * 1.05,
      label = paste0("Q4: ", format(round(qt4 / todos * 100, 2), nsmall = 2), "%")
    )
  print(quads)
  if (graphit == TRUE) {
    print(qqmod + geom_density_2d())
    if (saveit == TRUE) {
      ggsave(paste0(px, "_by_", py, "_mod.pdf"))
    }
    cat("\n")
    cat("call 'qqmod' to see plot with quadrants.")
    cat("\n")
  }
}

#---------------------------------------------------------------------------
# This function will remove certain points from the cells dataframe.
cull <- function() {
  print(sort(names(cells)))
  para <- readline(prompt = "Which parameter should be culled: ")
  b_or_s <- readline(prompt = "Data with values (b)igger, (s)maller, or (=) to the target: ")
  vale <- readline(prompt = "What is the target value: ")

  if (b_or_s == "b") {
    cells <<- cells[unlist(cells[, para]) < as.numeric(vale), ]
  } else if (b_or_s == "s") {
    cells <<- cells[unlist(cells[, para]) > as.numeric(vale), ]
  } else if (b_or_s == "=") {
    cells <<- cells <<- cells[unlist(cells[, para]) == as.numeric(vale), ]
  } else {
    grapho(cells)
  }
  grapho(cells)
  cat("\n")
}

#----------------------------------------------------------------
# Automatically pseudo-color the density of a cell-cycle profile.
get_density <- function(px,
                        py,
                        ...) {
  theme_set(theme_bw(base_size = 16))
  dens <<- MASS::kde2d(px, py, ...)
  ix <<- findInterval(px, dens$x)
  iy <<- findInterval(py, dens$y)
  ii <<- cbind(ix, iy)
  dens$z[ii]
}

color_density <- function() {
  cells$density <<- get_density(cells[px][, 1], cells[py][, 1], n = 100)
  qq <<- ggplot(data = cells, aes_string(px, py)) +
    geom_point(aes(color = density)) +
    scale_color_viridis() +
    xlab(xname) +
    ylab(yname)

  if (new_axis == TRUE) {
    changeAxis(store[1], store[2], store[3], store[4])
  }
  print(qq)
}


#----------------------------------------------------------------
# This function asks you to define a G1 peak and re-assigns that value to  '1'
normalizer <- function() {
  # create a column holder (just in case...)
  cells$edu <<- "Negative"
  cells$ploidy <<- "2N"
  if (TRUE %in% is.na(cells$LMean_ANC_edu)) {
    print(paste0(
      "Detected ", nrow(cells[is.na(cells$LMean_ANC_edu), ]),
      " NA values in LMean_ANC_edu (",
      round(100 * nrow(cells[is.na(cells$LMean_ANC_edu), ]) / nrow(cells), 2),
      "% of tota cells)."
    ))
    remove_nas <- readline(prompt = "Should these be removed? (Y/n) ")
    if (remove_nas != "n") {
      print("Removing rows with NA values")
      cells <<- cells[!is.na(cells$LMean_ANC_edu), ]
    }
  }
  if (TRUE %in% is.na(cells$log2_dna)) {
    print(paste0(
      "Detected ", nrow(cells[is.na(cells$log2_dna), ]),
      " NA values in log2_dna (",
      round(100 * nrow(cells[is.na(cells$log2_dna), ]) / nrow(cells), 2),
      "% of tota cells)."
    ))
    remove_nas <- readline(prompt = "Should these be removed? (Y/n) ")
    if (remove_nas != "n") {
      print("Removing rows with NA values")
      cells <<- cells[!is.na(cells$log2_dna), ]
    }
  }
  cat("Determining the G1 population now.")
  cat("\n")

  # This part calls an EdU negative popuation and creates the normalized values
  edu_values <- cells[
    !is.na(cells$LMean_ANC_edu) & !is.infinite(cells$LMean_ANC_edu),
  ]$LMean_ANC_edu
  big_bin <- which.max(density(edu_values)$y)
  edu_neg <- density(edu_values)$x[big_bin] + 0.1
  edu_hist <- ggplot(cells, aes(LMean_ANC_edu)) +
    geom_density() +
    geom_vline(xintercept = edu_neg) +
    xlab("Log EdU")
  print(edu_hist)
  g1_good <- readline(
    prompt = paste0(
      "Is ", format(round(edu_neg, 2), nsmall = 2),
      " representative of the EdU cutoff? (Y/n) "
    )
  )
  if (g1_good == "n") {
    edu_neg <- as.numeric(readline(prompt = "What value should EdU be cutoff as negative? "))
    edu_hist <- ggplot(cells, aes(LMean_ANC_edu)) +
      geom_density() +
      geom_vline(xintercept = edu_neg) +
      xlab("Log EdU")
    print(edu_hist)
  }
  cat(paste0("Using ", format(round(edu_neg, 2), nsmall = 2), " as the EdU cutoff."))
  cat("\n")

  cells[cells$LMean_ANC_edu > edu_neg & !is.na(cells$LMean_ANC_edu), ]$edu <<- "Positive"

  edu_neg_cells <<- subset(cells, LMean_ANC_edu < edu_neg)
  cells$edu_norm <<- (cells$LMean_ANC_edu + 1) - mean(edu_neg_cells$LMean_ANC_edu)

  # This part calls the 2N peak that is EdU-negative
  big_bin <- which.max(density(edu_neg_cells$log2_dna)$y)
  diploid <- density(edu_neg_cells$log2_dna)$x[big_bin]
  dna_hist <- ggplot(edu_neg_cells, aes(log2_dna)) +
    geom_density() +
    geom_vline(xintercept = diploid, color = "green", size = 1.5) +
    xlab("DNA content")
  print(dna_hist)
  g1_good <- readline(
    prompt = paste0(
      "Is ", format(round(diploid, 2), nsmall = 2),
      " representative of the 2N peak? (Y/n/manual) "
    )
  )
  if (g1_good == "n") {
    diploid <- as.numeric(readline(prompt = "What is the value of the 2N peak? "))
    dna_hist <- ggplot(edu_neg_cells, aes(log2_dna)) +
      geom_density() +
      geom_vline(xintercept = diploid, color = "green", size = 1.5) +
      xlab("DNA content")
    print(dna_hist)
  }


  if (g1_good == "manual") {
    diploid_old <- diploid
    diploid <- readline(prompt = "What is the value of the 2N peak? ")
    if (diploid == "") {
      diploid <- diploid_old
    } else {
      diploid <- as.numeric(diploid)
    }
    dna_hist <- ggplot(edu_neg_cells, aes(log2_dna)) +
      geom_density() +
      geom_vline(xintercept = diploid, color = "green", size = 1.5) +
      xlab("DNA content")
    print(dna_hist)
    peak_2n <- as.numeric(readline(prompt = "What should be the 2N - 4N trough value: "))

    dna_hist <- dna_hist + geom_vline(xintercept = peak_2n, color = "blue", size = 1.5)
    print(dna_hist)
    peak_4n <- as.numeric(readline(prompt = "What should be the 4N - >4N trough value: "))
    dna_hist <- dna_hist + geom_vline(xintercept = peak_4n, color = "red", size = 1.5)
    print(dna_hist)

    cells$dna_norm <<- (cells$log2_dna + 1) - diploid
    cells[cells$log2_dna > peak_2n, ]$ploidy <<- "4N"
    cells[cells$log2_dna > peak_4n, ]$ploidy <<- ">4N"
  } else {
    cat(paste0("Using ", format(round(diploid, 2), nsmall = 2), " as the 2N peak value."))
    cat("\n")
    # This part bins each point's ploidy
    cells$dna_norm <<- (cells$log2_dna + 1) - diploid
    cells[cells$dna_norm > 1.6, ]$ploidy <<- "4N"
    cells[cells$dna_norm > 2.6, ]$ploidy <<- ">4N"
  }


  ck_f <- readline(prompt = paste0("Is the filename ", thingee, "? (Y/n) "))
  if (ck_f != "n") {
    write_csv(cells, thingee)
  } else {
    thingee <<- readline(prompt = "What is the name of the file: ")
    write_csv(cells, thingee)
  }
  px <<- "dna_norm"
  py <<- "edu_norm"
  xname <<- "DNA content (Log 2)"
  yname <<- "EdU content (Log 10)"
  grapho(cells)
}


#-----------------------------------------------------------------
# Compensate for spectral overlap after graphing the offending colors against each other.
# IMPORTANT - The skewed color must be graphed on the y axis
# IMPORTANT - This should be obsolete if purgo.R has been run correctly
compensate <- function(df = cells) {
  # Generate a data frame to hold points for a possible future shape-fitting gate.
  cat(paste0(
    "Thank you for choosing coompensate(). To begin, lets assign a first point. ",
    "Be sure the intended line follows the overlap line:"
  ))
  # Then you pick a point on the line of the overlapping population
  point_one_x <- as.numeric(readline(prompt = "What is the x-value of the first point? "))
  point_one_y <- as.numeric(readline(prompt = "What is the y-value of the first point? "))
  print(point_one_x)
  print(point_one_y)
  point_one <- data.frame(x = point_one_x, y = point_one_y)
  print(point_one)
  cat("\n")
  print(qq + geom_point(data = point_one, aes(x = x, y = y, color = "red", size = 16)))
  cat(paste0("Great, I have added that point. Now lets move on to the second point"))

  # Then you pick a second point
  point_one_x <- as.numeric(readline(prompt = "What is the x-value of the second point? "))
  point_one_y <- as.numeric(readline(prompt = "What is the y-value of the second point? "))
  point_two <- data.frame(x = point_one_x, y = point_one_y)
  point_one <- rbind(point_one, point_two)

  print(qq + geom_point(data = point_one, aes(x = x, y = y, color = "red", size = 16)))

  y_or_n <- readline(prompt = "Does this look right (Y/n)? ")
  if (y_or_n == "n") {
    cat("Sorry, lets try that again.")
    cat("\n")
    compensate()
  } else {
    cat("Excellent, generating the linear equation now:")
    cat("\n")
    rise <- point_one[2, 2] - point_one[1, 2]
    run <- point_one[2, 1] - point_one[1, 1]
    slop <- rise / run
    yint <- point_one[1, 2] - (point_one[1, 1] * slop)
    print(
      qq + geom_point(data = point_one, aes(x = x, y = y, color = "red", size = 16)) +
        geom_abline(slope = slop, intercept = yint, size = 2, color = "blue")
    )
    cat(paste0("Slope has been determined to be: ", slop))
    cat("\n")
    cat(paste0("Y-intercept has been determined to be: ", yint))
    y_or_n <- readline(prompt = "Does this look right (Y/n)? ")
    if (y_or_n == "n") {
      cat("Sorry, lets try that again.")
      cat("\n")
      compensate()
    } else {
      cat("Applying compensation.")
      check <- df
      check[py] <- check[py] - (slop * check[px] + yint)
      cat("Regraphing:")
      grapho(check)
      y_or_n <- readline(prompt = "Does this look right (Y/n)? ")
      if (y_or_n == "n") {
        cat("Sorry, lets try that again.")
        cat("\n")
        compensate()
      } else {
        cells <<- check
        cat(paste0("Cells saved with compensated ", py))
        y_or_n <- readline(prompt = "Would you like to recalculate the means (y/n)? ")
        if (y_or_n == "y") {
          for (i in seq_len(ncol(cells))) {
            if (grepl(py, names(cells)[i])) {
              cells[paste0("I", names(cells)[i])] <<- cells[, i] * cells$Area
            }
          }
          for (i in seq_len(ncol(cells))) {
            if (grepl(py, names(cells)[i])) {
              cells[paste0("L", names(cells)[i])] <<- log(cells[, i] + 1, 10)
            }
          }
          cat("Means recalculated. Thank you for your patience and patronage.")
          write_csv(cells, thingee)
        }
      }
    }
  }
}

#------------------------------------------------------

# This function take the nuclear locations and sizes and re-maps the original image
#   df is the dataframe, defaulted to 'cells'
re_map <- function(df = cells,
                   # highlight contains the Number(s) to highlight in the remapped image
                   highlight = FALSE,
                   # X and Y begin the position names; defaults use nuclear geometry
                   y = "X_ANC_",
                   x = "Y_ANC_",
                   #   S is the size of the ROI, defaulted to nucelar size
                   size = "Area_ANC_",
                   #   Xn and Yn are the axis names and don't really need to be changed
                   x_name = "X position",
                   y_name = "Y position",
                   # colz supplies one color for each highlighted population
                   colz = cb_black) {
  # These set the variable names using regex to anchor tot he front of the variable name
  x <- names(df)[str_detect(names(df), paste0("^", x))]
  y <- names(df)[str_detect(names(df), paste0("^", y))]
  size <- names(df)[str_detect(names(df), paste0("^", size))]
  # This creates a ggplot object for graphing
  test <- ggplot(data = df, aes(y = -unlist(df[x]), x = unlist(df[y]), size = unlist(df[size])))
  # If the highlight variable is a vector and is not F, its graphed here
  if (is.vector(highlight) && highlight != FALSE) {
    check <- df
    check$highlight <- check$Number %in% highlight
    qq <<- test + geom_point(aes(color = check$highlight)) +
      ylab(y_name) +
      xlab(x_name) +
      scale_color_manual(values = colz) +
      theme_classic() +
      theme(axis.text = element_text(face = "bold", size = 16), axis.ticks = element_blank()) +
      theme(plot.title = element_text(size = 40, face = "bold")) +
      theme(axis.title = element_text(size = 30, face = "bold"))
  } else if (length(highlight) != 1 && highlight != FALSE) {
    # If highlight is not a vector (such as a variable), graph it with color.
    qq <<- test + geom_point(aes(color = highlight)) +
      ylab(y_name) +
      xlab(x_name) +
      theme_classic() +
      theme(axis.text = element_text(face = "bold", size = 16), axis.ticks = element_blank()) +
      theme(plot.title = element_text(size = 40, face = "bold")) +
      theme(axis.title = element_text(size = 30, face = "bold")) +
      scale_color_manual(values = colz)
  } else {
    # If highlight is FALSE, graph it without a color option.
    qq <<- test + geom_point() +
      ylab(y_name) +
      xlab(x_name) +
      theme_classic() +
      theme(axis.text = element_text(face = "bold", size = 16), axis.ticks = element_blank()) +
      theme(plot.title = element_text(size = 40, face = "bold")) +
      theme(axis.title = element_text(size = 30, face = "bold"))
  }
  # The plot is graphed
  print(qq)
}

#-----------------------------------------------------------------
# Runs either a PCA of t-SNE
run_ca <- function(df = cells,
                   plot_type = "tsne",
                   save_file = "figures/pcaPlot.png",
                   subz = FALSE,
                   groupz = TRUE) {
  cellist <- df

  # Get rid of NAs and Infs
  cellist[is.na(cellist)] <- 0
  cellist[cellist == Inf] <- 0
  cellist[cellist == -Inf] <- 0

  if (groupz == TRUE) {
    elly <- TRUE
    print(names(cellist))
    grope <- readline(
      prompt = paste0(
        "What should the ", plot_type,
        " be grouped by from the 'cellist' dataframe: "
      )
    )
    group_list <- unlist(df[grope])
  }
  # get rid of non-numerics
  num <- unlist(lapply(df, is.numeric))
  cellist_nums <- df[, num]

  # Get only the nuclear columns & remove NA
  if (subz != FALSE) {
    cellist_nums <- cellist_nums[
      str_detect(names(cellist_nums), subz) | str_detect(names(cellist_nums), "Number")
    ]
    cellist <- cellist[str_detect(names(cellist), subz) | str_detect(names(cellist), "Number")]
  }

  # Get rid of non-unique rows
  cellist <- unique(cellist)
  cellist_nums <- unique(cellist_nums)

  if (plot_type != "pca" && plot_type != "tsne") {
    plot_type <- readline(prompt = "Are you trying to make a tsne or pca plot: ")
  }

  # Make the pca object
  if (plot_type == "pca") {
    # Get rid of columns with no variance
    x <- foo(cellist_nums)
    cellist_nums <- cellist_nums[, -x]
    cellist_nums[is.na(cellist_nums)] <- 0

    cell_pca <- prcomp(cellist_nums, center = TRUE, scale = TRUE)
    percentage <- round(cell_pca$sdev / sum(cell_pca$sdev) * 100, 2)

    if (groupz == TRUE) {
      pca_out <- data.frame(x = cell_pca$x[, 1], y = cell_pca$x[, 2], groupBy = cellist[, grope])
      pp <<- ggplot(data = pca_out, aes(x = x, y = y)) +
        geom_point(aes(color = groupBy), size = 3, alpha = 0.5) +
        geom_density_2d(aes(color = groupBy)) +
        xlab(paste0("PC1 (", percentage[1], " %)")) +
        ylab(paste0("PC2 (", percentage[2], " %)")) +
        theme_classic2()
    } else {
      pca_out <- data.frame(x = cell_pca$x[, 1], y = cell_pca$x[, 2])
      pp <<- ggplot(data = pca_out, aes(x = x, y = y)) +
        geom_point(size = 3, alpha = 0.5) +
        geom_density_2d() +
        xlab(paste0("PC1 (", percentage[1], " %)")) +
        ylab(paste0("PC2 (", percentage[2], " %)")) +
        theme_classic2()
    }
  } else {
    cellist_m <- as.matrix(cellist_nums)
    cellist_m[is.na(cellist_m)] <- 0
    tsne_out <- Rtsne::Rtsne(cellist_m, perplexity = 100)
    if (groupz == TRUE) {
      tsne_plot <- data.frame(x = tsne_out$y[, 1], y = tsne_out$y[, 2], groupBy = cellist[, grope])
      pp <<- ggplot(data = tsne_plot, aes(x = x, y = y)) +
        geom_point(aes(color = groupBy), size = 3, alpha = 0.5) +
        geom_density_2d(aes(color = groupBy)) +
        theme_classic2()
    } else {
      tsne_plot <- data.frame(x = tsne_out$y[, 1], y = tsne_out$y[, 2])
      pp <<- ggplot(data = tsne_plot, aes(x = x, y = y)) +
        geom_point(size = 3, alpha = 0.5) +
        geom_density_2d() +
        theme_classic2()
    }
  }
  print(pp)

  # Report how much of the data is shown
  how_many_left <- round(
    (nrow(cellist_nums) * ncol(cellist_nums)) / (nrow(df) * ncol(df)) * 100,
    digits = 1
  )
  cat(paste0("Approximately ", how_many_left, "% of the data is represented after trimming."))
  cat("\n")


  ggsave(save_file, dpi = 300, height = 4, width = 4)
}

foo <- function(dat) {
  out <- lapply(dat, function(x) length(unique(x)))
  want <- which(!out > 1)
  unlist(want)
}

#-----------------------------------------------------------------
# Find the nearest cell in another population.
nearest_n <- function(home = cells, target = cells, label = "CellToCell") {
  cells[paste0("nearest_", label)] <<- "NotHome"
  cells[paste0("nearest_", label, "_distance")] <<- 0
  for (i in seq_len(nrow(cells))) {
    if (cells$Number[i] %in% unique(home$Number)) {
      pool <- subset(target, Number != cells$Number[i])
      pool$distance <- sqrt(
        (pool$X_ANC_dna - cells$X_ANC_dna[i])**2 +
          (pool$Y_ANC_dna - cells$Y_ANC_dna[i])**2
      )
      nearest <- subset(pool, distance == min(pool$distance))
      cells[paste0("nearest_", label)][i, ] <<- nearest$Number[1]
      cells[paste0("nearest_", label, "_distance")][i, ] <<- nearest$distance[1]
    }
  }
}

#-----------------------------------------------------------------
# Cluster cells using the nearestN principle and a gated variable.
# Sort other cells by distance, then record selected-category cells until reaching the edge.
# Record the size and cell-number IDs of the cluster centered around each cell.
# After characterizing the clusters of each cell, it then compares the various clusters
cluster_count <- function(df = cells,
                          clust_pop = "tag",
                          center = "Positive",
                          edge = "Negative",
                          cluster_cutoff = 0.5) {
  cells["cluster_size"] <<- 0
  cells["cluster_mates"] <<- "None"
  cells["cluster_id"] <<- "None"
  cluster_id <- 1
  home <- cells[cells[clust_pop] == center, ]
  for (i in seq_len(nrow(cells))) {
    if (cells$Number[i] %in% unique(home$Number)) {
      pool <- subset(cells, Number != cells$Number[i])
      pool$distance <- sqrt(
        (pool$X_ANC_dna - cells$X_ANC_dna[i])**2 +
          (pool$Y_ANC_dna - cells$Y_ANC_dna[i])**2
      )
      pool <- pool[order(pool$distance), ]
      cluster_size <- 1
      cluster_mates <- c(cells$Number[i])
      for (j in seq_len(nrow(pool))) {
        if (pool[clust_pop][j, ] != edge) {
          cluster_size <- cluster_size + 1
          cluster_mates <- append(cluster_mates, as.character(pool[j, ]$Number))
        } else {
          break
        }
      }
      cluster_mates <- sort(cluster_mates)
      cells["cluster_size"][i, ] <<- cluster_size
      cells["cluster_mates"][i, ] <<- paste(unlist(cluster_mates), collapse = "_")
    }
  }
  for (i in seq_len(nrow(cells))) {
    if (cells$Number[i] %in% unique(home$Number & cells$cluster_id[i] == "None")) {
      pool <- subset(home, Number != cells$Number[i])
      pool$lScore <-
        p_cluster <- subset(pool, lScore > cluster_cutoff)
      cells$cluster_id[cells$Number %in% p_cluster$Number, ] <- cluster_id
      cluster_id <- cluster_id + 1
    }
  }
}

#-----------------------------------------------------------------
#-----------------------------------------------------------------
# Color palettes:
color_blind <- c(
  "#999999", "#E69F00", "#56B4E9", "#009E73",
  "#F0E442", "#0072B2", "#D55E00", "#CC79A7"
)
cb_black <- c(
  "#000000", "#E69F00", "#56B4E9", "#009E73",
  "#F0E442", "#0072B2", "#D55E00", "#CC79A7"
)

# The cellCycle_colors palette are the hexdecimal codes for:
# G1, S, G2/M, mitosis, and premature mitosis, respectively, from Justice et al., 2019
cell_cycle_colors <- c("#D4D4D4", "#98C84C", "#23B8CC", "#F16B1A", "#E5001C")


#------------------------------------------------------------------------------------------------
