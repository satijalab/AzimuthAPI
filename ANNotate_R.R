#!/usr/bin/env Rscript 
# nolint start



########################################################################
##### hardware configuration ###########################################
########################################################################
########################################################################


#### function to check availability of an nvidia gpu through nvidia-smi
if_gpu <- function() {
  if (.Platform$OS.type == "windows") {
    # Windows-specific check
    gpu_check <- tryCatch(
      system("where nvidia-smi", intern = TRUE, ignore.stderr = TRUE),
      error = function(e) character(0)
    )
  } else {
    # Unix-like systems (Linux, macOS)
    gpu_check <- tryCatch(
      system("which nvidia-smi", intern = TRUE, ignore.stderr = TRUE),
      error = function(e) character(0)
    )
  }
  
  return(length(gpu_check) > 0)
}




########################################################################
####### environment and dependencies ###################################
########################################################################
########################################################################


#### creating and setting up conda environment
if (!requireNamespace("reticulate", quietly = TRUE) || 
    packageVersion("reticulate") < "1.40.0") {
  message(
    "Installing or updating reticulate to version 1.40.0 "
    "alongwith dependencies Rcpp and RcppTOML..."
    )
  install.packages("Rcpp")
  install.packages("RcppTOML")
  install.packages("reticulate")
  
  if (!requireNamespace("reticulate", quietly = TRUE) || 
      packageVersion("reticulate") < "1.40.0") {
    stop(
      "Failed to install or update reticulate to "
      "version 1.40.0. Please resolve this issue externally."
      )
  } else {
    message("reticulate successfully updated to version 1.40.0.")
  }
} 

library(reticulate)

#### create conda env and install dependencies
setup_conda_env <- function(
  yml_file, 
  requirements_file = NULL, 
  conda_path = NULL, 
  force_create = FALSE,
  tensorflow_gpu = "tensorflow[and-cuda]==2.17",
  tensorflow_cpu = "tensorflow-cpu==2.17.0"
) {
  # Validate conda installation
  if (is.null(conda_path)) {
    conda_path <- reticulate::conda_binary()
  }
  if (is.null(conda_path)) {
    stop("Conda not found. Please install conda or miniconda and try again.")
  }
  
  # Parse and validate YAML environment file
  library(yaml)
  if (!file.exists(yml_file)) {
    stop(sprintf("YML file does not exist at the path: %s", yml_file))
  }
  yml_data <- tryCatch(
    yaml::read_yaml(yml_file),
    error = function(e) {
      stop(sprintf("Error reading YAML file: %s", e$message))
    }
  )
  env_name <- yml_data$name
  
  # Check if environment exists
  env_exists <- env_name %in% reticulate::conda_list()$name
  
  # Handle force_create flag
  if (force_create && env_exists) {
    message(sprintf(
      "Conda environment '%s' exists, but force_create is TRUE. Deleting it...",
      env_name
    ))
    system2(conda_path, c("env", "remove", "--name", env_name, "--yes"))
    message(sprintf("Environment '%s' deleted successfully.", env_name))
    env_exists <- FALSE
  }
  
  # Create environment if it doesn't exist
  if (!env_exists) {
    message(sprintf(
      "Creating conda environment '%s' from '%s'...", 
      env_name, 
      yml_file
    ))
    system2(conda_path, c("env", "create", "-f", yml_file))
    message(sprintf("Environment '%s' created.", env_name))
  }
  
  # Activate the environment
  message(sprintf("Activating environment '%s'...", env_name))
  reticulate::use_condaenv(
    condaenv = env_name,
    conda = conda_path,
    required = TRUE
  )
  message(sprintf("Environment '%s' is now active.", env_name))
  
  # Verify environment connection
  python_binary <- ifelse(
    .Platform$OS.type == "windows",
    file.path("Scripts", "python.exe"),
    "bin/python"
  )
  
  active_env <- reticulate::py_config()$python
  expected_env <- file.path(
    dirname(dirname(conda_path)), 
    "envs", 
    env_name, 
    python_binary
  )
  
  verify_environment <- function() {
    tryCatch({
      if (normalizePath(active_env) != normalizePath(expected_env)) {
        stop(sprintf(
          paste0(
            "Error: Reticulate is not connected to the expected conda ",
            "environment properly.\nExpected: %s\nActive: %s"
          ),
          expected_env, 
          active_env
        ))
      }
      message(sprintf("Reticulate is connected to: %s", active_env))
    }, error = function(e) {
      warning(sprintf(
        paste0(
          "Could not verify environment paths due to: %s\n",
          "Raw paths - Expected: %s, Active: %s"
        ), 
        e$message, expected_env, active_env
      ))
    })
  }
  
  verify_environment()
  
  # Install additional requirements if environment was just created
  if (!env_exists && 
      !is.null(requirements_file) && 
      file.exists(requirements_file)) {
    
    message(sprintf("Installing dependencies from %s...", requirements_file))
    
    # Install requirements file
    system2(conda_path, c(
      "run", "-n", env_name, "pip", "install", "-r", requirements_file
    ))
    
    # Upgrade pip
    system2(conda_path, c(
      "run", "-n", env_name, "pip", "install", "--upgrade", "pip"
    ))
    
    # Upgrade keras
    system2(conda_path, c(
      "run", "-n", env_name, "pip", "install", "--upgrade", "keras"
    ))
    
    # Install appropriate TensorFlow version
    if (if_gpu()) {
      message("GPU detected, installing GPU version of TensorFlow...")
      system2(conda_path, c(
        "run", "-n", env_name, "pip", "install", tensorflow_gpu
      ))
    } else {
      message("No GPU detected, installing CPU version of TensorFlow...")
      system2(conda_path, c(
        "run", "-n", env_name, "pip", "install", tensorflow_cpu
      ))
    }
    
    message("All dependencies installed through pip.")
  } else if (!env_exists && 
             (is.null(requirements_file) || !file.exists(requirements_file))) {
    message("Requirements file not provided or doesn't exist.")
  }
  
  # Print final Python configuration
  print(reticulate::py_config())
}

yml_path <- "utils/env_min.yml" 
requirements_path <- "utils/requirements_min.txt" 

#conda_binary <- "/path/to/conda/bin/conda" 
setup_conda_env(yml_file = yml_path, requirements_file = requirements_path)


#### R dependencies
if (!requireNamespace("Seurat", quietly = TRUE)) {
  stop("Seurat is not installed. Please install it to proceed.")
} else {
  seurat_version <- as.numeric(substr(
    as.character(packageVersion("Seurat")),
    1,
    1
    ))  
}

library(Seurat)
library(Matrix)


#### python dependencies
home_dir <- ifelse(.Platform$OS.type == "windows", 
                   Sys.getenv("USERPROFILE"), 
                   file.path("/home", Sys.getenv("USER")))
python_module_path <- file.path(home_dir, "panhumanpy/")
python_module_src <- file.path(home_dir, "panhumanpy/src")

py_run_string(paste(
  "import sys; sys.path.append('", python_module_src, "')", 
  sep = ""
  ))
annotate <- import("panhumanpy.ANNotate")
sp <- import("scipy.sparse")




########################################################################
######### util functions ###############################################
########################################################################
########################################################################




get_data <- function(object, assay, layer= "data") {
  if (packageVersion("Seurat") >= "5.0.0") {
    return(LayerData(object[[assay]], layer = layer))
  } else {
    assay_obj <- object@assays[[assay]]
    return(slot(assay_obj, layer))
  }
}


########################################################################


read_obj_min <- function(query_obj, feature_names_col, assay_default='RNA') {
  
  if (!(assay_default %in% names(query_obj@assays))){
    stop(paste("Assay", assay_default, "not present in the seurat object."))
  }
  
  DefaultAssay(query_obj) <- assay_default
  
  if (packageVersion("Seurat") >= "5.0.0"){
    layers_obj <- Layers(query_obj, assay = assay_default)
  } else{
    layers_obj <- slotNames(query_obj[[assay_default]])
  }
  
  if ("data" %in% layers_obj) {
    normalized_data <- get_data(query_obj, assay = assay_default)
  } else {
    query_obj <- NormalizeData(query_obj[[assay_default]])
    normalized_data <- get_data(query_obj, assay = assay_default)
  }
  
  X_query <- t(normalized_data)
  X_query <- Matrix(X_query, sparse = TRUE)
  
  cell_metadata <- query_obj@meta.data
  query_cells_df <- as.data.frame(cell_metadata)
  
  if (!is.null(feature_names_col)) {
    feature_metacols <- colnames(query_obj[[assay_default]][[]])
    if (feature_names_col %in% feature_metacols){
      query_features <- query_obj[[assay_default]][[feature_names_col]]
      query_features <- as.list(query_features[[feature_names_col]])
    } else {
      stop(paste(
        feature_names_col, 
        "not found as a column in the df returned by "
        "object[[",assay_default,"]][[]]"
        ))
    }
    
  } else {
    query_features <- as.list(rownames(normalized_data))
  }
  return(list(
    X_query = X_query, 
    query_features = query_features, 
    query_cells_df = query_cells_df
    ))
}

########################################################################


package_obj <- function(
  extract_embeddings, 
  embeddings_dict, 
  umap_embeddings, 
  umap_embeddings_dict, 
  query_cells_df, 
  query_obj
  ) {
  
  if (extract_embeddings) {
    for (em_name in names(embeddings_dict)) {
      
      em_matrix <- as.matrix(embeddings_dict[[em_name]])
      
      if (nrow(em_matrix) != length(Cells(query_obj))) {
        stop(paste(
          "Dimension mismatch:", em_name, " does not have as many "
          "cells as the query obj."
          ))
      }
      
      rownames(em_matrix) <- Cells(query_obj)
      
      dimreduc_obj <- CreateDimReducObject(
        embeddings = em_matrix, 
        key = paste0(em_name, "_"), 
        assay = DefaultAssay(query_obj)
        )
      
      query_obj[[paste0(em_name)]] <- dimreduc_obj
    }
  }
  
  if (umap_embeddings) {
    for (em_name in names(umap_embeddings_dict)) {
      
      em_matrix <- as.matrix(umap_embeddings_dict[[em_name]])
      
      if (nrow(em_matrix) != length(Cells(query_obj))) {
        stop(paste(
          "Dimension mismatch:", em_name, " does not "
          "have as many cells as the query obj."
          ))
      }
      
      rownames(em_matrix) <- Cells(query_obj)
      
      dimreduc_obj <- CreateDimReducObject(
        embeddings = em_matrix, 
        key = paste0(em_name, "_"), 
        assay = DefaultAssay(query_obj)
        )
      
      query_obj[[paste0("umapANN", em_name)]] <- dimreduc_obj
    }
  }
  
  query_obj@meta.data <- as.data.frame(query_cells_df)  
  
  return(query_obj)
}


########################################################################


PrepLabel <- function(
  object, 
  label_id = 'final_level_label', 
  newid = 'PrepLabel', 
  cutid = 'Other', 
  cutoff=10
  ) {
  rejected_names <- names(which(table(object@meta.data[,label_id])<cutoff))
  object@meta.data[,newid]=as.character(object@meta.data[,label_id])
  rejected_cells <- which(object@meta.data[,label_id]%in%rejected_names)
  object@meta.data[rejected_cells,newid]=cutid
  return(object)
}


########################################################################


ANNotate <- function(
  query_obj,
  feature_names_col = NULL,
  annotation_pipeline = 'supervised',
  eval_batch_size = 40960,
  normalization_override = FALSE,
  norm_check_batch_size = 1000,
  output_mode = 'minimal',
  refine_labels = TRUE,
  extract_embeddings = TRUE,
  umap_embeddings = TRUE,
  n_neighbors = 30,
  n_components = 2,
  metric = "cosine",
  min_dist = 0.3,
  umap_lr = 1.0,
  umap_seed = 42,
  spread = 1.0,
  verbose = TRUE,
  init = "spectral",
  process_obj = TRUE,
  cutoff_abs = 5,
  cutoff_frac = 0.001
) {
  options(warn = -1)
  cat("Running Pan-Human Azimuth:\n")
  cat("\n")
  
  # Convert integers
  eval_batch_size <- as.integer(eval_batch_size)
  n_neighbors <- as.integer(n_neighbors)
  n_components <- as.integer(n_components)
  umap_seed <- as.integer(umap_seed)
  
  # Read and process the Seurat object
  query <- read_obj_min(query_obj, feature_names_col)
  X_query <- sp$csr_matrix(r_to_py(query$X_query))
  query_features <- query$query_features
  cells_meta <- query$query_cells_df
  
  # Run annotation core
  core_outputs <- annotate$annotate_core(
    X_query,
    query_features,
    cells_meta,
    annotation_pipeline,
    eval_batch_size,
    normalization_override,
    norm_check_batch_size,
    output_mode,
    refine_labels,
    extract_embeddings,
    umap_embeddings,
    n_neighbors, 
    n_components, 
    metric, 
    min_dist, 
    umap_lr, 
    umap_seed, 
    spread,
    verbose,
    init
  )
  
  
  embeddings_dict <- core_outputs$embeddings_dict
  umap_embeddings_dict <- core_outputs$umap_dict
  cells_meta <- core_outputs$cells_meta
  
  annotated_obj = package_obj(
    extract_embeddings, 
    embeddings_dict, 
    umap_embeddings, 
    umap_embeddings_dict, 
    cells_meta, 
    query_obj
    )
  
  
  
  if (process_obj){
    annotated_obj <- PrepLabel(
      annotated_obj,
      label_id = 'final_level_label',
      cutoff = min(cutoff_abs, cutoff_frac*ncol(annotated_obj)),
      cutid = 'Other',
      newid = 'azimuth_label'
      )
    Idents(annotated_obj) <- 'azimuth_label'
  }
  
  
  
  return(annotated_obj)
  
}



########################################################################
########## functions for command-line usage ############################
########################################################################
########################################################################



read_seurat_object <- function(filepath, assay_name = "RNA") {
  # Check if file exists
  if (!file.exists(filepath)) {
    stop(paste("File does not exist at the specified path:", filepath))
  }
  
  # Check file extension
  if (!grepl("\\.rds$", filepath, ignore.case = TRUE)) {
    warning("File does not have .rds extension. Attempting to read anyway.")
  }
  
  # Load the RDS file
  tryCatch({
    seurat_obj <- readRDS(filepath)
  }, error = function(e) {
    stop(paste("Error reading RDS file:", e$message))
  })
  
  # Check if it's a Seurat object
  if (!inherits(seurat_obj, "Seurat")) {
    stop("The file does not contain a Seurat object.")
  }
  
  # Check Seurat version
  seurat_version <- packageVersion("Seurat")
  if (seurat_version < "4.4.0") {
    warning(paste(
      "Current Seurat version:", seurat_version, 
      "is below 4.4.0. Some functionality may "
      "not work as expected."
      ))
  }
  
  # Check if specified assay exists
  if (!assay_name %in% names(seurat_obj@assays)) {
    available_assays <- paste(names(seurat_obj@assays), collapse = ", ")
    stop(paste("Assay", assay_name, "not found in the Seurat object.",
               "Available assays:", available_assays))
  }
  
  # Set default assay
  DefaultAssay(seurat_obj) <- assay_name
  
  return(seurat_obj)
}



########################################################################

save_seurat_object <- function(seurat_obj, filepath) {
  # Check if object is a Seurat object
  if (!inherits(seurat_obj, "Seurat")) {
    stop("The object is not a Seurat object.")
  }
  
  # Construct output filepath
  dir_path <- dirname(filepath)
  file_name <- basename(filepath)
  file_name_ann <- sub("\\.rds$", "_ANN.rds", file_name)
  if (file_name == file_name_ann) {
    # If no .rds extension was found, add it
    file_name_ann <- paste0(file_name, "_ANN.rds")
  }
  output_path <- file.path(dir_path, file_name_ann)
  
  # This is a joke about massively dieting the object
  # so that we don't redownload expression or embedding data
  # in this API case we know the object has an RNA assay with counts and
  # data layers
  keto_object=TRUE
  if (keto_object) {
    if ("RNA" %in% names(seurat_obj@assays)) {
      # Handle Seurat v5+ and v4 differently
      if (packageVersion("Seurat") >= "5.0.0") {
        # For Seurat v5+
        if ("data" %in% names(seurat_obj[["RNA"]]@layers)) {
          seurat_obj[
            ["RNA"]
            ]$data <- Matrix::Matrix(
              0, 
              nrow = nrow(seurat_obj[["RNA"]]$data),
              ncol = ncol(seurat_obj[["RNA"]]$data),
              sparse = TRUE
              )
        }
        if ("counts" %in% names(seurat_obj[["RNA"]]@layers)) {
          seurat_obj[
            ["RNA"]
            ]$counts <- Matrix::Matrix(
              0, 
              nrow = nrow(seurat_obj[["RNA"]]$counts),
              ncol = ncol(seurat_obj[["RNA"]]$counts),
              sparse = TRUE
              )
        }
      } else {
        # For Seurat v4
        if ("data" %in% slotNames(seurat_obj[["RNA"]])) {
          slot(
            seurat_obj[["RNA"]], 
            "data"
            ) <- Matrix::Matrix(
              0, 
              nrow = nrow(slot(seurat_obj[["RNA"]], "data")),
              ncol = ncol(slot(seurat_obj[["RNA"]], "data")),
              sparse = TRUE
            )
        }
        if ("counts" %in% slotNames(seurat_obj[["RNA"]])) {
          slot(
            seurat_obj[["RNA"]], 
            "counts"
            ) <- Matrix::Matrix(
              0, 
              nrow = nrow(slot(seurat_obj[["RNA"]], "counts")),
              ncol = ncol(slot(seurat_obj[["RNA"]], "counts")),
              sparse = TRUE
            )
        }
      }
    }
    
    # Remove ANNshallow_embeddings
    if ("azimuth_embed" %in% names(seurat_obj@reductions)) {
      seurat_obj[["azimuth_embed"]] <- NULL
    }
  }
  
  # Save the object
  message("Saving annotated object to ", output_path)
  saveRDS(seurat_obj, file = output_path)
  
  # Return the path invisibly
  invisible(output_path)
}


########################################################################
############################ main executable function ##################

ANNotate_R <- function() {
  # Parse command line arguments
  args <- parse_annotate_args()
  formatted_args <- format_annotate_args(args)
  
  # Check if running interactively or from command line
  is_interactive <- interactive()
  
  if (!is_interactive) {
    message("ANNotate_R is being executed from the command line...")
  } else {
    message("ANNotate_R is being run interactively.")
  }
  
  # Extract filepath and read the Seurat object
  filepath <- formatted_args$filepath
  assay_name <- formatted_args$assay
  
  message(paste("Reading Seurat object from", filepath))
  query_obj <- read_seurat_object(filepath, assay_name)
  
  # Filter arguments to only include those expected by ANNotate
  # Get ANNotate function arguments
  annotate_args <- formals(ANNotate)
  annotate_arg_names <- names(annotate_args)
  
  # Filter formatted_args to only include arguments that ANNotate expects
  args_for_annotate <- formatted_args[names(formatted_args) %in% annotate_arg_names]
  
  # Add query_obj to the arguments
  args_for_annotate$query_obj <- query_obj
  
  # Print the function call for debugging
  if (formatted_args$verbose) {
    message("Calling ANNotate function with the following arguments:")
    annotate_call <- paste0("ANNotate(", 
                           paste(sprintf("%s = %s", 
                                        names(args_for_annotate), 
                                        sapply(args_for_annotate, function(x) {
                                          if (is.character(x)) {
                                            return(paste0("'", x, "'"))
                                          } else if (is.null(x)) {
                                            return("NULL")
                                          } else {
                                            return(as.character(x))
                                          }
                                        })), 
                                 collapse = ", "), 
                           ")")
    message(annotate_call)
  }
  
  # Extract filepath from arguments
  filepath <- formatted_args$filepath
  
  # Call ANNotate function with the filtered arguments
  message("Running ANNotate function...")
  annotated_obj <- do.call(ANNotate, args_for_annotate)
  
  # Save the annotated object
  save_seurat_object(seurat_obj = annotated_obj, filepath = filepath)
  
  # Return the annotated object
  message("ANNotate_R completed successfully.")
  return(annotated_obj)
}

# Execute the main function if run from command line
if (!interactive()) {
  ANNotate_R()
}
