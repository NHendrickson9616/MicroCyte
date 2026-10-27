#!/usr/bin/Rscript

# Generate the project directories described by schema.csv.
dir_gen <- function(id_type = "name",
                    single_image = "auto",
                    image_number = 5,
                    number_of_targets = 5,
                    project_dir = ".") {
  if (!is.character(project_dir) || length(project_dir) != 1L || !dir.exists(project_dir)) {
    stop("project_dir must identify an existing directory.")
  }

  id_type <- tolower(id_type)
  single_image <- tolower(single_image)
  if (length(id_type) != 1L || !id_type %in% c("name", "number")) {
    stop("id_type must be either 'name' or 'number'.")
  }
  if (length(single_image) != 1L || !single_image %in% c("auto", "single")) {
    stop("single_image must be either 'auto' or 'single'.")
  }
  if (length(image_number) != 1L || !is.numeric(image_number) || is.na(image_number) ||
        image_number < 1 || image_number != as.integer(image_number)) {
    stop("image_number must be a positive whole number.")
  }
  if (length(number_of_targets) != 1L || !is.numeric(number_of_targets) ||
        is.na(number_of_targets) || number_of_targets < 1 ||
        number_of_targets != as.integer(number_of_targets)) {
    stop("number_of_targets must be a positive whole number.")
  }

  project_dir <- normalizePath(project_dir, mustWork = TRUE)
  schema_path <- file.path(project_dir, "schema.csv")
  if (!file.exists(schema_path)) {
    stop("No schema.csv found in project_dir: ", project_dir)
  }

  output_dirs <- file.path(project_dir, c("files", "figures", "data"))
  for (output_dir in output_dirs) {
    if (!dir.exists(output_dir) && !dir.create(output_dir, recursive = TRUE)) {
      stop("Could not create directory: ", output_dir)
    }
  }

  schema <- read.csv(schema_path)
  if (!"notebook_id" %in% names(schema)) {
    stop("schema.csv must contain a notebook_id column.")
  }
  if (nrow(schema) == 0L) {
    stop("schema.csv must contain at least one experiment row.")
  }

  starting_point <- as.integer(number_of_targets) + 2L
  if (id_type == "number" || ncol(schema) < starting_point) {
    x_list <- unique(schema$notebook_id)
  } else {
    x_list <- NULL
    for (columnIndex in seq.int(starting_point, ncol(schema))) {
      if (names(schema)[columnIndex] != "name_id") {
        if (is.null(x_list)) {
          x_list <- schema[, columnIndex]
        } else {
          x_list <- paste(x_list, schema[, columnIndex], sep = "_")
        }
        schema$name_id <- x_list
        write.csv(schema, schema_path, row.names = FALSE)
      }
    }
    if (is.null(x_list)) {
      stop("schema.csv has no columns available to construct name_id values.")
    }
  }

  if (anyNA(x_list) || any(!nzchar(as.character(x_list)))) {
    stop("Directory identifiers in schema.csv must not be missing or empty.")
  }

  condition_dirs <- file.path(project_dir, "files", as.character(x_list))
  for (condition_dir in condition_dirs) {
    if (!dir.exists(condition_dir) && !dir.create(condition_dir, recursive = TRUE)) {
      stop("Could not create directory: ", condition_dir)
    }
    if (single_image == "single") {
      image_dirs <- file.path(condition_dir, paste0("image_", seq_len(image_number)))
      for (image_dir in image_dirs) {
        if (!dir.exists(image_dir) && !dir.create(image_dir, recursive = TRUE)) {
          stop("Could not create directory: ", image_dir)
        }
      }
    }
  }

  invisible(condition_dirs)
}

run_dir_gen_cli <- function(project_dir = ".", prompt = readline) {
  if (!is.function(prompt)) {
    stop("prompt must be a function.")
  }

  ask <- function(message) prompt(message)

  run_type <- tolower(trimws(ask(
    "Generating directories. Should this be done by name or number (NAME/number): "
  )))
  if (!nzchar(run_type)) {
    run_type <- "name"
  }

  single_images <- tolower(trimws(ask(
    "Are you taking single images, or automated capture (AUTO/single): "
  )))
  if (!nzchar(single_images)) {
    single_images <- "auto"
  }

  images_number <- 5
  if (single_images == "single") {
    image_input <- suppressWarnings(as.numeric(trimws(ask(
      "How many images per condition (Default = 5): "
    ))))
    if (!is.na(image_input)) {
      images_number <- image_input
    }
    print(paste0("Generating ", images_number, " image folders in each condition..."))
  }

  target_input <- suppressWarnings(as.numeric(trimws(ask(
    "How many targets are listed in the schema (Default 5): "
  ))))
  number_of_targets <- if (is.na(target_input)) 5 else target_input

  invisible(dir_gen(
    id_type = run_type,
    single_image = single_images,
    image_number = images_number,
    number_of_targets = number_of_targets,
    project_dir = project_dir
  ))
}

if (sys.nframe() == 0L) {
  run_dir_gen_cli()
}
