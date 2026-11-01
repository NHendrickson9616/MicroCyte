#!/usr/bin/Rscript

.microcyte_stage_files <- c(
  setup = "0_setup.R",
  dir_gen = "1_dirGen.R",
  purgo = "2_purgo.R",
  re_name = "3_reName.R",
  ima_gen = "5_imaGen.R",
  concato = "6_concato.R",
  unite = "7_unite.R"
)

prompt_user <- function(prompt = "", input_connection = NULL) {
  if (interactive() && is.null(input_connection)) {
    return(base::readline(prompt = prompt))
  }

  cat(prompt)
  close_connection <- FALSE
  if (is.null(input_connection)) {
    input_connection <- file("stdin", open = "r")
    close_connection <- TRUE
  }
  if (close_connection) {
    on.exit(close(input_connection), add = TRUE)
  }

  response <- readLines(input_connection, n = 1L, warn = FALSE)
  if (length(response) == 0L) {
    stop("Standard input ended before the workflow completed.")
  }
  response
}

load_microcyte_stages <- function(project_dir = ".",
                                  envir = new.env(parent = parent.frame())) {
  project_dir <- normalizePath(project_dir, mustWork = TRUE)
  script_dir <- file.path(project_dir, "bin", "scripts", "scripts")
  script_paths <- file.path(script_dir, .microcyte_stage_files)
  missing_scripts <- script_paths[!file.exists(script_paths)]
  if (length(missing_scripts) > 0L) {
    stop("Missing workflow scripts: ", paste(missing_scripts, collapse = ", "))
  }

  for (scriptPath in script_paths) {
    sys.source(scriptPath, envir = envir)
  }
  envir
}

collect_unite_options <- function(project_dir = ".", prompt = readline) {
  full_response <- tolower(trimws(prompt(
    "Should the entire dataset be used (FULL), or equivalent cells per sample (partial): "
  )))
  use_full_dataset <- full_response != "partial"
  sample_size <- 1000
  if (!use_full_dataset) {
    sample_size <- suppressWarnings(as.numeric(trimws(prompt(
      "How many cells should be pulled: "
    ))))
    if (length(sample_size) != 1L || is.na(sample_size) || sample_size < 1 ||
          sample_size != as.integer(sample_size)) {
      stop("The number of cells must be a positive whole number.")
    }
  }

  grouping <- tolower(trimws(prompt(
    "Should samples be united on a variable (VAR), or altogether (all): "
  )))
  if (grouping == "all") {
    return(list(full = use_full_dataset, sample_size = sample_size, patterns = NULL))
  }

  grouping_source <- tolower(trimws(prompt(
    "Is the variable a part of the (name), or a (SCHEMA) variable: "
  )))
  if (grouping_source == "name") {
    pattern <- trimws(prompt("What is the exact name you wish to unite: "))
    if (!nzchar(pattern)) {
      stop("The name pattern must not be empty.")
    }
    return(list(full = use_full_dataset, sample_size = sample_size, patterns = pattern))
  }

  schema_path <- file.path(project_dir, "schema.csv")
  schema <- read.csv(schema_path)
  print(names(schema))
  schema_name <- trimws(prompt("Which schema variable would you like to unite by: "))
  if (!schema_name %in% names(schema)) {
    stop("Unknown schema variable: ", schema_name)
  }
  patterns <- unique(as.character(schema[[schema_name]]))
  patterns <- patterns[!is.na(patterns) & nzchar(patterns)]
  if (length(patterns) == 0L) {
    stop("The selected schema variable contains no usable values: ", schema_name)
  }

  list(full = use_full_dataset, sample_size = sample_size, patterns = patterns)
}

run_unite_options <- function(stages, options, project_dir = ".") {
  if (is.null(options$patterns)) {
    return(invisible(stages$unite(
      full = options$full,
      sample_size = options$sample_size,
      project_dir = project_dir
    )))
  }

  outputs <- lapply(options$patterns, function(pattern) {
    stages$unite(
      full = options$full,
      pattern = pattern,
      sample_size = options$sample_size,
      project_dir = project_dir,
      fixed = TRUE
    )
  })
  invisible(outputs)
}

run_unite_cli <- function(stages, project_dir = ".", prompt = readline) {
  print("Uniting sample datasets...")
  options <- collect_unite_options(project_dir = project_dir, prompt = prompt)
  outputs <- run_unite_options(stages, options, project_dir = project_dir)
  print("Sample datasets combined")
  invisible(outputs)
}

ask_for_normalization <- function(prompt = readline, default = FALSE) {
  default_label <- if (default) "y/n" else "y/N"
  response <- tolower(trimws(prompt(paste0(
    "Should the numeric values be normalized among image sets (", default_label, ")? "
  ))))
  if (default) response != "n" else response == "y"
}

run_concato_cli <- function(stages,
                            project_dir = ".",
                            prompt = readline,
                            offer_unite = TRUE) {
  print("Concatenating image sets...")
  normalization <- ask_for_normalization(prompt = prompt, default = FALSE)
  outputs <- stages$concato(normalization = normalization, project_dir = project_dir)
  print("Image sets combined")

  if (offer_unite) {
    continue <- tolower(trimws(prompt("Would you like to now run Unite (y/n)? ")))
    if (continue != "n") {
      run_unite_cli(stages, project_dir = project_dir, prompt = prompt)
    }
  }
  invisible(outputs)
}

run_ima_gen_cli <- function(stages,
                            project_dir = ".",
                            prompt = readline,
                            offer_concato = TRUE) {
  print("Running ImaGen now...")
  outputs <- stages$run_ima_gen(project_dir = project_dir, prompt = prompt)
  print("Datasets combined.")

  if (offer_concato) {
    continue <- tolower(trimws(prompt("Would you like to now run concato (y/n)? ")))
    if (continue != "n") {
      run_concato_cli(stages, project_dir = project_dir, prompt = prompt)
    }
  }
  invisible(outputs)
}

run_full_workflow <- function(stages, project_dir = ".", prompt = readline) {
  print("Running the entire process")
  normalization <- ask_for_normalization(prompt = prompt, default = TRUE)
  unite_options <- collect_unite_options(project_dir = project_dir, prompt = prompt)

  generate_directories <- tolower(trimws(prompt("Should directories be generated (y/n)? ")))
  if (generate_directories != "n") {
    print("Generating directories...")
    stages$run_dir_gen_cli(project_dir = project_dir, prompt = prompt)
    print("Directories complete")
    print("Please take your images and continue when ready...")
    prompt("Continue: ")
  }

  images_taken <- tolower(trimws(prompt("Have images been taken (y/n)? ")))
  if (images_taken != "n") {
    print("Beginning the purge...")
    stages$run_purgo_cli(project_dir = project_dir, prompt = prompt)
    print("Purge complete")

    print("Renaming images per the schema file...")
    stages$re_name(project_dir = project_dir)
    print("Images renamed. Please run YggData and continue when finished...")
    prompt("Continue: ")
  }

  run_ima_gen_cli(stages, project_dir = project_dir, prompt = prompt, offer_concato = FALSE)

  print("Concatenating image sets...")
  stages$concato(normalization = normalization, project_dir = project_dir)
  print("Image sets combined")

  print("Uniting sample datasets...")
  outputs <- run_unite_options(stages, unite_options, project_dir = project_dir)
  print("Sample datasets combined")
  invisible(outputs)
}

load_analysis_scripts <- function(project_dir = ".", envir = parent.frame()) {
  if (!requireNamespace("MicroCyte", quietly = TRUE)) {
    stop(
      "Install the MicroCyte package before loading analysis tools. ",
      "Run `R CMD INSTALL .` from: ", normalizePath(project_dir, mustWork = TRUE),
      call. = FALSE
    )
  }

  exports <- getNamespaceExports("MicroCyte")
  for (name in exports) {
    assign(name, getExportedValue("MicroCyte", name), envir = envir)
  }
  message("Loaded ", length(exports), " MicroCyte analysis functions.")
  invisible(exports)
}

run_microcyte <- function(project_dir = ".", prompt = NULL, analysis_env = parent.frame()) {
  project_dir <- normalizePath(project_dir, mustWork = TRUE)

  input_connection <- NULL
  if (is.null(prompt)) {
    if (!interactive()) {
      input_connection <- file("stdin", open = "r")
      on.exit(close(input_connection), add = TRUE)
    }
    prompt <- function(message = "") {
      prompt_user(message, input_connection = input_connection)
    }
  }
  if (!is.function(prompt)) {
    stop("prompt must be a function.")
  }

  stages <- load_microcyte_stages(project_dir = project_dir)
  missing_packages <- stages$missing_microcyte_packages()
  if (length(missing_packages) > 0L) {
    cat("Missing packages:", paste(missing_packages, collapse = ", "), "\n")
    cat("Choose option 0 to install them before running the workflow.\n")
  }

  scripts <- data.frame(
    Number = 0:9,
    Script_Name = c(
      "Setup",
      "dir_gen",
      "purgo",
      "re_name",
      "ima_gen",
      "concato",
      "unite",
      "full run",
      "reset names",
      "analyze"
    )
  )
  print(paste(
    "Welcome to MicroCyte. Be sure to have edited your schema file before beginning.",
    "Here are your options:"
  ))
  print(scripts)

  repeat {
    option <- tolower(trimws(prompt(paste0(
      "What would you like to run (0-", max(scripts$Number), ", or 'end'): "
    ))))

    if (option == "end") {
      print("Thank you for using MicroCyte.")
      break
    }

    if (option == "0") {
      stages$install_microcyte_packages()
    } else if (option == "1") {
      print("Generating directories...")
      stages$run_dir_gen_cli(project_dir = project_dir, prompt = prompt)
      print("Directories complete")
    } else if (option == "2") {
      print("Beginning the purge...")
      stages$run_purgo_cli(project_dir = project_dir, prompt = prompt)
      print("Purge complete")
      continue <- tolower(trimws(prompt("Would you like to now run re_name (y/n)? ")))
      if (continue != "n") {
        print("Renaming images per the schema file...")
        stages$re_name(project_dir = project_dir)
        print("Images renamed. Please run YggData to continue.")
      }
    } else if (option == "3") {
      print("Renaming images per the schema file...")
      stages$re_name(project_dir = project_dir)
      print("Images renamed. Please run YggData to continue.")
    } else if (option == "4") {
      run_ima_gen_cli(stages, project_dir = project_dir, prompt = prompt)
    } else if (option == "5") {
      run_concato_cli(stages, project_dir = project_dir, prompt = prompt)
    } else if (option == "6") {
      run_unite_cli(stages, project_dir = project_dir, prompt = prompt)
    } else if (option == "7") {
      run_full_workflow(stages, project_dir = project_dir, prompt = prompt)
    } else if (option == "8") {
      stages$reset_names(project_dir = project_dir)
    } else if (option == "9") {
      load_analysis_scripts(project_dir = project_dir, envir = analysis_env)
      cat("\n")
      print("Good luck with your analysis!")
      break
    } else {
      print("Sorry, I didn't catch that. Please try again, or type 'end' to finish.")
      print(scripts)
    }
  }

  invisible(stages)
}

if (sys.nframe() == 0L) {
  run_microcyte()
}
