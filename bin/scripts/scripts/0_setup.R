#!/usr/bin/Rscript

MICROCYTE_REQUIRED_PACKAGES <- c( # nolint: object_name_linter.
  "ggplot2",
  "ggpubr",
  "ggtern",
  "tidyverse",
  "MASS",
  "viridis",
  "Rtsne",
  "stringr",
  "imager",
  "magick",
  "RecordLinkage",
  "hash",
  "plotly"
)


missing_microcyte_packages <- function(packages = MICROCYTE_REQUIRED_PACKAGES) {
  setdiff(packages, rownames(installed.packages()))
}

install_microcyte_packages <- function(packages = MICROCYTE_REQUIRED_PACKAGES,
                                       repos = getOption("repos")) {
  print("Checking installed packages...")
  missing_packages <- missing_microcyte_packages(packages)

  if (length(missing_packages) == 0L) {
    print("All required packages are installed.")
    return(invisible(character()))
  }

  print(paste("Installing:", paste(missing_packages, collapse = ", ")))

  if (is.null(repos) || identical(unname(repos["CRAN"]), "@CRAN@")) {
    repos <- c(CRAN = "https://cloud.r-project.org")
  }

  is_mac <- identical(Sys.info()[["sysname"]], "Darwin")
  is_cran_mac_r <- is_mac && grepl("^/Library/Frameworks/R.framework", R.home())

  # Avoid linking Homebrew R packages against incompatible Anaconda libraries.
  if (is_mac && grepl("anaconda", Sys.which("xml2-config"), ignore.case = TRUE)) {
    path_entries <- strsplit(Sys.getenv("PATH"), .Platform$path.sep, fixed = TRUE)[[1]]
    path_entries <- path_entries[!grepl("anaconda", path_entries, ignore.case = TRUE)]
    Sys.setenv(PATH = paste(path_entries, collapse = .Platform$path.sep))
  }

  install_args <- list(pkgs = missing_packages, repos = repos, dependencies = NA)
  if (is_cran_mac_r) {
    install_args$type <- "binary"
  }
  do.call(install.packages, install_args)

  still_missing <- missing_microcyte_packages(packages)
  if (length(still_missing) > 0L) {
    stop("Package installation failed for: ", paste(still_missing, collapse = ", "))
  }

  print("Installation complete.")
  invisible(missing_packages)
}

if (sys.nframe() == 0L) {
  install_microcyte_packages()
}
