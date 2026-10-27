pull_n_zip <- function(image_folder = "./files/",
                       schema = "./schema.csv",
                       outputname = "raw_images") {
  if (!outputname %in% list.dirs()) {
    dir.create(outputname)
    run_it <- TRUE
  } else {
    cat(paste0("File name '", outputname, "' already present. Choose another name for the output."))
    run_it <- FALSE
  }
  for (i in list.files(path = image_folder)) {
    dir.create(paste0(outputname, "/", i))
    for (j in list.files(pattern = "image_", path = paste0(image_folder, "/", i))) {
      dir.create(paste0(outputname, "/", i, "/", j))
      for (k in list.files(path = paste0(image_folder, "/", i, "/", j), pattern = "tif")) {
        file.copy(
          from = paste0(image_folder, "/", i, "/", j, "/", k),
          to = paste0(outputname, "/", i, "/", j, "/", k)
        )
      }
    }
  }
  file.copy(schema, outputname)
  zip(paste0(outputname, ".zip"), outputname)
}
