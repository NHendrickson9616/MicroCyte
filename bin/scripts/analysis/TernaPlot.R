# Welcome to TernaPlot_v2

# This script takes a csv file and generates a ternary plot.

# For help customizing this script, contact Jason Needhan at jneedham@uab.edu.

#--------------------------------------------------------------------------------------

# Requires the ggplot2 and ggtern packages.

#--------------------------------------------------------------------------------------

# Opens the ternary package
library("ggplot2")
library("ggtern")

#--------------------------------------------------------------------------------------
# Asks for the target CSV, opens, reports detected column names, and imports data
ternaplot <- function(file_name = datum, customize_it = TRUE, save_it = TRUE) {
  if (is.data.frame(file_name)) {
    things <<- file_name
  } else if (!is.data.frame(file_name) && grepl(".csv", file_name)) {
    things <<- read.csv(file = thingee, header = TRUE)
  } else {
    print(paste(
      "Neither a csv nor a dataframe was supplied.",
      "This script will now fail spectacularly. Please try again but better..."
    ))
  }
  print("Data categories detected:")
  print(names(things))
  samples <<- readline(prompt = "Which category is the sample ID: ")
  cat("\n")
  lname <<- readline(prompt = "Which category is left point: ")
  llabel <<- readline(prompt = "What should the left label be: ")
  cat("\n")
  tname <<- readline(prompt = "Which category is top point: ")
  tlabel <<- readline(prompt = "What should the top label be: ")
  cat("\n")
  rname <<- readline(prompt = "Which category is right point: ")
  rlabel <<- readline(prompt = "What should the right label be: ")
  cat("\n")

  # Generates a dataframe of the csv
  intera <<- data.frame(
    sampleName = things[samples],
    lpoint = things[, lname],
    tpoint = things[, tname],
    rpoint = things[, rname]
  )
  print("Dataframe generated.")
  cat("\n")

  # Makes a ternary plot called 'tex'
  tex <<- ggtern(data = intera, aes(x = lpoint, y = tpoint, z = rpoint)) +
    geom_point(size = 3) +
    theme_bw() +
    Tlab(tlabel) +
    Rlab(rlabel) +
    Llab(llabel)

  print(tex)
  if (customize_it == TRUE) {
    customize()
    tex <<- texq
  } else {
    cat("\n")
    print("Call customize() to make things fancy.")
  }
  if (save_it == TRUE) {
    save_png()
  } else {
    print("Call savePNG() to save as a high-resolution PNG")
    cat("\n")
  }
}

#-------------------------------------------------------------------------------------------------
# Saves plot as high-resolution PNG

save_png <- function() {
  fm <- readline(prompt = "What name should the image file be saved as: ")
  ftif <- paste0(fm, ".png")
  print(tex)
  ggsave(ftif, dpi = 300, units = "in", width = 11, height = 8.5)
  print("PNG of tex has been saved.")
}

#-------------------------------------------------------------------------------------------------
# Makes a customized ternary plot with colors, shapes, and adds a 95% confidencen variable
# The 95% confidence variable groups by both the color and shape discernments.

customize <- function() {
  print("Data categories detected:")
  print(names(things))

  colz <- readline(prompt = "Would you like to color the dots based on a variable (y/n): ")
  if (colz == "y") {
    colv <- readline(prompt = "Which category should color be based on: ")
    coln <- readline(prompt = "What is the name of this variable: ")
    intera$coln <- things[, colv]
  }
  cat("\n")
  shapez <- readline(
    prompt = "Would you like to alter shape based on a categorical variable (y/n): "
  )
  if (shapez == "y") {
    shapev <- readline(prompt = "Which category should shape be based on: ")
    shapen <- readline(prompt = "What is the name of this variable: ")
    intera$shapen <- things[, shapev]
  }
  cat("\n")
  conz <- readline(prompt = "Would you like to add a confidence interval (y/n): ")
  if (conz == "y") {
    conv <- as.integer(readline(prompt = "What percentage of confidence should be used: ")) / 100
  }

  if (colz == "y" && shapez == "y") {
    # Makes a ternary plot called 'tex'
    texq <<- ggtern(
      data = intera,
      aes(x = lpoint, y = tpoint, z = rpoint, color = coln, shape = shapen)
    ) +
      geom_point(size = 3) +
      Tlab(tlabel) +
      Rlab(rlabel) +
      Llab(llabel) +
      guides(shape = guide_legend(shapen), color = guide_legend(coln))
  } else if (colz == "y" && shapez != "y") {
    # Makes a ternary plot called 'tex'
    texq <<- ggtern(data = intera, aes(x = lpoint, y = tpoint, z = rpoint, color = coln)) +
      geom_point(size = 3) +
      Tlab(tlabel) +
      Rlab(rlabel) +
      Llab(llabel) +
      guides(color = guide_legend(coln))
  } else if (colz != "y" && shapez == "y") {
    # Makes a ternary plot called 'tex'
    texq <<- ggtern(data = intera, aes(x = lpoint, y = tpoint, z = rpoint, shape = shapen)) +
      geom_point(size = 3) +
      Tlab(tlabel) +
      Rlab(rlabel) +
      Llab(llabel) +
      guides(shape = guide_legend(shapen))
  } else {
    # Makes a ternary plot called 'tex'
    texq <<- ggtern(data = intera, aes(x = lpoint, y = tpoint, z = rpoint)) +
      geom_point(size = 3) +
      Tlab(tlabel) +
      Rlab(rlabel) +
      Llab(llabel)
  }

  if (conz == "y") {
    texq <<- texq + geom_confidence_tern(breaks = c(conv))
  }

  texq <<- texq + theme_cycle

  print(texq)
  print("Call savetif() to save as a high-resolution tiff.")
  cat("\n")
  print("Call customize() to make things fancy in a different way.")
  cat("\n")
}

# Cell-cycle colors for G1, S, G2/M, mitosis, and premature mitosis,
# respectively, as used in Justice et al., 2019.
cell_cycle_colors <- c("#D4D4D4", "#98C84C", "#23B8CC", "#F16B1A", "#E5001C")

theme_cycle <- theme_custom(
  base_size = 20,
  col.T = cell_cycle_colors[2],
  col.L = cell_cycle_colors[1],
  col.R = cell_cycle_colors[3],
  tern.panel.background = "white",
  tern.plot.background = NULL
)
