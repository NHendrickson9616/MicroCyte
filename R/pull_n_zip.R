pull_n_zip <- function(image_folder = "./files/",
                       schema = "./schema.csv",
                       outputname = "raw_images") {
  if (!dir.exists(image_folder)) {
    stop("`image_folder` does not exist: ", image_folder, call. = FALSE)
  }
  if (!file.exists(schema)) {
    stop("`schema` does not exist: ", schema, call. = FALSE)
  }
  if (dir.exists(outputname) || file.exists(outputname)) {
    stop(
      "Output path already exists; choose another `outputname`: ",
      outputname,
      call. = FALSE
    )
  }

  dir.create(outputname, recursive = TRUE, showWarnings = FALSE)
  copied_files <- character()
  sample_names <- list.files(image_folder)

  for (sample_name in sample_names) {
    sample_source <- file.path(image_folder, sample_name)
    if (!dir.exists(sample_source)) {
      next
    }

    image_names <- list.files(sample_source, pattern = "^image_")
    for (image_name in image_names) {
      image_source <- file.path(sample_source, image_name)
      if (!dir.exists(image_source)) {
        next
      }

      image_destination <- file.path(outputname, sample_name, image_name)
      dir.create(image_destination, recursive = TRUE, showWarnings = FALSE)
      tif_names <- list.files(image_source, pattern = "\\.tiff?$", ignore.case = TRUE)
      for (tif_name in tif_names) {
        source_file <- file.path(image_source, tif_name)
        destination_file <- file.path(image_destination, tif_name)
        copied <- file.copy(source_file, destination_file, overwrite = FALSE)
        if (!isTRUE(copied)) {
          stop("Failed to copy image: ", source_file, call. = FALSE)
        }
        copied_files <- c(copied_files, destination_file)
      }
    }
  }

  schema_destination <- file.path(outputname, basename(schema))
  if (!isTRUE(file.copy(schema, schema_destination, overwrite = FALSE))) {
    stop("Failed to copy schema: ", schema, call. = FALSE)
  }

  archive <- paste0(outputname, ".zip")
  zip_status <- utils::zip(zipfile = archive, files = outputname, flags = "-r9X")
  if (!identical(zip_status, 0L) || !file.exists(archive)) {
    stop("Failed to create archive: ", archive, call. = FALSE)
  }

  invisible(list(
    output_directory = outputname,
    archive = archive,
    schema = schema_destination,
    images = copied_files
  ))
}
