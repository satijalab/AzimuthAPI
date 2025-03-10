#!/usr/bin/env Rscript 
# nolint start


#### function to check availability of an nvidia gpu through nvidia-smi
if_gpu <- function() {
  gpu_check <- system("which nvidia-smi", intern = TRUE, ignore.stderr = TRUE)
  
  if (length(gpu_check) == 0) {
    return(FALSE)
  } else{
    return(TRUE)
  }
}


#### creating and setting up conda environment
if (!requireNamespace("reticulate", quietly = TRUE) || 
    packageVersion("reticulate") < "1.40.0") {
  message("Installing or updating reticulate to version 1.40.0 alongwith dependencies Rcpp and RcppTOML...")
  install.packages("Rcpp")
  install.packages("RcppTOML")
  install.packages("reticulate")
  
  if (!requireNamespace("reticulate", quietly = TRUE) || 
      packageVersion("reticulate") < "1.40.0") {
    stop("Failed to install or update reticulate to version 1.40.0. Please resolve this issue externally.")
  } else {
    message("reticulate successfully updated to version 1.40.0.")
  }
} else {
  message("reticulate is already installed and up to date.")
}
library(reticulate)

#### create conda env and install dependencies
setup_conda_env <- function(yml_file, requirements_file, conda_path=NULL, force_create=FALSE) {
  
  if (is.null(conda_path)) {
    conda_path <- reticulate::conda_binary()
  }
  if (is.null(conda_path)) {
    stop("Conda not found. Please install conda or miniconda and try again.")
  }
  
  env_name <- "AzimuthNN_test" 
  
  env_check_command <- sprintf('%s env list | grep "%s"', conda_path, env_name)
  env_check <- system(env_check_command, intern = TRUE)
  
  if (force_create && length(env_check) > 0) {
    message(sprintf("Conda environment '%s' exists, but force_create is TRUE. Deleting it first...", env_name))
    system2(conda_path, c("env", "remove", "--name", env_name, "--yes"))
    message(sprintf("Environment '%s' deleted successfully.", env_name))
    env_check <- character(0) 
  }
  
  
  if (length(env_check) == 0) {
    message(sprintf("Conda environment '%s' does not exist. Creating it first...", env_name))
    
    if (!file.exists(yml_file)) {
      stop(sprintf("YML file does not exist at the path: %s", yml_file))
    }
    
    library(yaml)
    yml_data <- yaml::read_yaml(yml_file)
    env_name <- yml_data$name
    
    message(sprintf("Creating conda environment '%s' from '%s'...", env_name, yml_file))
    
    cmd <- c("env", "create", "-f", yml_file) 
    message(sprintf("Running command: %s", paste(cmd, collapse = " ")))
    
    
    system2(conda_path, cmd)
    
    message(sprintf("Environment '%s' created successfully.", env_name))
    
    reticulate::use_condaenv(condaenv = env_name, conda = conda_path, required = TRUE)
    message(sprintf("Environment '%s' is now active.", env_name))
    
    active_env <- reticulate::py_config()$python
    expected_env <- file.path(dirname(dirname(conda_path)), "envs", env_name, "bin", "python")
    
    if (normalizePath(active_env) != normalizePath(expected_env)) {
      stop(sprintf(
        "Error: Reticulate is not connected to the expected conda environment properly.\nExpected: %s\nActive: %s",
        expected_env, active_env
      ))
    } else {
      message(sprintf("Reticulate is correctly connected to: %s", active_env))
    }
    
    if (!is.null(requirements_file) && file.exists(requirements_file)) {
      message(sprintf("Installing pip dependencies from %s...", requirements_file))
      system2(conda_path, c("run", "-n", env_name, "pip", "install", "-r", requirements_file))
      system2(conda_path, c("run", "-n", env_name, "pip", "install", "--upgrade", "pip"))
      system2(conda_path, c("run", "-n", env_name, "pip", "install", "--upgrade", "keras"))
      if (if_gpu()){
        system2(conda_path, c("run", "-n", env_name, "pip", "install", "tensorflow[and-cuda]==2.17"))
      }else{
        system2(conda_path, c("run", "-n", env_name, "pip", "install", "tensorflow-cpu==2.17.0"))
      }
      
      message("All pip dependencies installed successfully.")
    } else {
      message("No requirements file provided or found.")
    }
    
    print(reticulate::py_config())
  } else {
    message(sprintf("Conda environment '%s' found.", env_name))
    
    reticulate::use_condaenv(condaenv = env_name, conda = conda_path, required = TRUE)
    message(sprintf("Environment '%s' is now active.", env_name))
    
    active_env <- reticulate::py_config()$python
    expected_env <- file.path(dirname(dirname(conda_path)), "envs", env_name, "bin", "python")
    
    if (normalizePath(active_env) != normalizePath(expected_env)) {
      stop(sprintf(
        "Error: Reticulate is not connected to the expected conda environment properly.\nExpected: %s\nActive: %s",
        expected_env, active_env
      ))
    } else {
      message(sprintf("Reticulate is correctly connected to: %s", active_env))
    }
    
    print(reticulate::py_config())
  }
  
}

yml_path <- "utils/env_min.yml" 
requirements_path <- "utils/requirements_min.txt" 

#conda_binary <- "/path/to/conda/bin/conda" 
setup_conda_env(yml_file = yml_path, requirements_file = requirements_path)


#### R dependencies
library(Seurat)
library(Matrix)

ensure_argparse <- function() {
  if (!requireNamespace("argparse", quietly = TRUE)) {
    message("The 'argparse' package is not installed. Installing now...")
    install.packages("argparse")
    message("'argparse' package installed successfully.")
  } else {
    message("'argparse' package is already installed.")
  }
}

ensure_argparse()
library(argparse)

#### other python dependencies
python_module_path <- file.path("/home", Sys.getenv("USER"), "panhumanpy/src")
py_run_string(paste("import sys; sys.path.append('", python_module_path, "')", sep = ""))
annotate <- import("panhumanpy.core.ANNotate")
sp <- import("scipy.sparse")

#### functions

get_data <- function(object, assay, layer= "data") {
  if (packageVersion("Seurat") >= "5.0.0") {
    return(LayerData(object[[assay]], layer = layer))
  } else {
    assay_obj <- object@assays[[assay]]
    return(slot(assay_obj, layer))
  }
}

read_obj_R <- function(query_filepath, feature_names_col) {
  query_obj <- readRDS(query_filepath)
  if ("data" %in% names(query_obj@assays$RNA)) {
    normalized_data <- LayerData(query_obj,layer = 'data')
  } else {
    query_obj <- NormalizeData(query_obj)
    normalized_data <- LayerData(query_obj,layer = 'data')
  }
  X_query <- t(normalized_data)
  X_query <- Matrix(X_query, sparse = TRUE)
  cell_metadata <- query_obj@meta.data
  query_cells_df <- as.data.frame(cell_metadata)
  
  if (!is.null(feature_names_col)) {
    query_features <- as.list(query_obj@assays$RNA@meta.features[[feature_names_col]])
  } else {
    query_features <- as.list(rownames(normalized_data))
  }
  
  return(list(X_query = X_query, query_features = query_features, query_cells_df = query_cells_df, query_obj = query_obj))
}


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
      stop(paste(feature_names_col, "not found as a column in the df returned by object[[",assay_default,"]][[]]"))
    }
    
  } else {
    query_features <- as.list(rownames(normalized_data))
  }
  return(list(X_query = X_query, query_features = query_features, query_cells_df = query_cells_df))
}

package_obj <- function(embeddings_mode, embeddings_dict, if_umap_embeddings, umap_embeddings_dict, query_cells_df, query_obj) {
  
  if (!is.null(embeddings_mode)) {
    for (em_name in names(embeddings_dict)) {
      
      em_matrix <- as.matrix(embeddings_dict[[em_name]])
      
      if (nrow(em_matrix) != length(Cells(query_obj))) {
        stop(paste("Dimension mismatch:", em_name, " does not have as many cells as the query obj."))
      }
      
      rownames(em_matrix) <- Cells(query_obj)
      
      dimreduc_obj <- CreateDimReducObject(embeddings = em_matrix, key = paste0("ANN", em_name, "_"), assay = DefaultAssay(query_obj))
      
      query_obj[[paste0("ANN", em_name)]] <- dimreduc_obj
    }
  }
  
  if (if_umap_embeddings) {
    for (em_name in names(umap_embeddings_dict)) {
      
      em_matrix <- as.matrix(umap_embeddings_dict[[em_name]])
      
      if (nrow(em_matrix) != length(Cells(query_obj))) {
        stop(paste("Dimension mismatch:", em_name, " does not have as many cells as the query obj."))
      }
      
      rownames(em_matrix) <- Cells(query_obj)
      
      dimreduc_obj <- CreateDimReducObject(embeddings = em_matrix, key = paste0("umapANN", em_name, "_"), assay = DefaultAssay(query_obj))
      
      query_obj[[paste0("umapANN", em_name)]] <- dimreduc_obj
    }
  }
  
  query_obj@meta.data <- as.data.frame(query_cells_df)  
  
  return(query_obj)
}

PrepLabel <- function(object, label_id = 'final_level_label', newid = 'PrepLabel', cutid = 'Other', cutoff=10) {
  rejected_names <- names(which(table(object@meta.data[,label_id])<cutoff))
  object@meta.data[,newid]=as.character(object@meta.data[,label_id])
  rejected_cells <- which(object@meta.data[,label_id]%in%rejected_names)
  object@meta.data[rejected_cells,newid]=cutid
  return(object)
}




ANNotate <- function(
                    query_obj,
                    feature_names_col=NULL,
                    source_data_dir="/data/kfold_data",
                    features_txt="features.txt",
                    split_mode="cumulative",
                    model="M0.2",
                    loss_name="level_wt_focal_loss",
                    epochs=55,
                    train_seed=100,
                    data_seed=414,
                    data_source="data/kfold_data/datasets/fold10_02_26_2025_17_53_139",
                    data_split=c(7,1,2),
                    mask_seed=NULL,
                    tm_frac=NULL,
                    lm_frac=NULL,
                    save=TRUE,
                    batch_size=256,
                    eval_batch_size=40960,
                    optimizer_name="adam",
                    lr=NULL,
                    l1=NULL,
                    l2=0.01,
                    dropout=0.1,
                    normalization_override=FALSE,
                    embeddings_mode="shallow",
                    if_knn_scores=FALSE,
                    if_umap_embeddings=TRUE,
                    if_refine_labels=TRUE,
                    n_neighbors=30,
                    n_components=2,
                    metric="cosine",
                    min_dist=0.3,
                    umap_lr=1.0,
                    umap_seed=42,
                    spread=1.0,
                    verbose=TRUE,
                    init="spectral",
                    object_disk=FALSE,
                    out_file_disk=FALSE,
                    process_obj=TRUE,
                    cutoff_abs=5,
                    cutoff_frac=0.001){

  options(warn = -1)
  
  cat("Running Pan-Human Azimuth:\n")
  cat("\n")
  
  #### make sure that integers are passed to python as integers
  epochs <- as.integer(epochs)
  train_seed <- as.integer(train_seed)
  data_seed <- as.integer(data_seed)
  data_split <- lapply(data_split, as.integer)
  if (!is.null(mask_seed)){
    mask_seed <- as.integer(mask_seed)
  }
  batch_size <- as.integer(batch_size)
  eval_batch_size <- as.integer(eval_batch_size)
  n_neighbors <- as.integer(n_neighbors)
  n_components <- as.integer(n_components)
  umap_seed <- as.integer(umap_seed)

   
  
  #### reading the seurat object 
  query <- read_obj_min(query_obj, feature_names_col)
  X_query <- sp$csr_matrix(r_to_py(query$X_query))
  query_features <- query$query_features
  query_cells_df <- query$query_cells_df

  
  
  
  #### run annotate_core
  core_outputs <- annotate$annotate_core(
    X_query,
    query_features,
    source_data_dir,
    features_txt,
    split_mode,
    model,
    epochs,
    train_seed,
    loss_name,
    data_seed,
    data_source,
    data_split,
    mask_seed,
    tm_frac,
    lm_frac,
    batch_size,
    optimizer_name,
    lr,
    l1,
    l2, 
    dropout,
    save,
    eval_batch_size,
    normalization_override,
    embeddings_mode,
    query_cells_df,
    if_knn_scores,
    if_umap_embeddings,
    if_refine_labels,
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
  
  
  embeddings_mode <- core_outputs[[3]]
  embeddings_dict <- core_outputs[[4]]
  query_cells_df <- core_outputs[[10]]
  if_umap_embeddings <- core_outputs[[11]]
  umap_embeddings_dict <- core_outputs[[12]]
  
  annotated_obj = package_obj(embeddings_mode, embeddings_dict, if_umap_embeddings, umap_embeddings_dict, query_cells_df, query_obj)
  
  
  
  if (process_obj){
    annotated_obj <- PrepLabel(annotated_obj,label_id = 'final_level_label',cutoff = min(cutoff_abs, cutoff_frac*ncol(annotated_obj)),cutid = 'Other',newid = 'azimuth_label')
    Idents(annotated_obj) <- 'azimuth_label'
  }
  
  
  
  return(annotated_obj)
  
}

arg_parse_in_R <- function() {
  cat("Capturing arguments passed to R scipt...")
  cat("\n")
  
  parser <- ArgumentParser(description = "Argument parser for the script")
  
  parser$add_argument(
    "filepath",
    help = "Enter abs file path to the query. Query should be in h5ad format.",
    type = "character"
  )
  
  parser$add_argument(
    "-fn", "--feature_names_col",
    default = NULL,
    help = "Enter the column name where the feature names are stored in query.var where query is the anndata object read from the h5ad.",
    type = "character"
  )
  
  parser$add_argument(
    "-sdd", "--source_data_dir",
    default = "/brahms/sarkars/AzimuthNN_clone/AzimuthNN/sarkars/data/dataset_main",
    help = "Source data directory",
    type = "character"
  )
  
  parser$add_argument(
    "-ft", "--features_txt",
    default = "features_02_26_25_17_50.txt",
    help = "Features text file",
    type = "character"
  )
  
  parser$add_argument(
    "-m", "--mode",
    default = "cumulative",
    help = "Enter the label mode: cumulative or independent",
    type = "character"
  )
  
  parser$add_argument(
    "-md", "--model",
    default = "M0.2",
    help = "Enter the model name",
    type = "character"
  )
  
  parser$add_argument(
    "-l", "--loss",
    default = "level_wt_focal_loss",
    help = "Enter the loss function used for optimization",
    type = "character"
  )
  
  parser$add_argument(
    "--epochs",
    default = 55,
    help = "Enter the number of epochs the model has been trained for",
    type = "integer"
  )
  
  parser$add_argument(
    "-ts", "--train_seed",
    default = 100,
    help = "Enter the training seed used",
    type = "integer"
  )
  
  parser$add_argument(
    "-ds", "--data_seed",
    default = 414,
    help = "Enter the data prep seed that was used",
    type = "integer"
  )
  
  parser$add_argument(
    "-dso", "--data_source",
    default = "data/kfold_data/datasets/fold10_02_26_2025_17_53_139",
    help = "Enter the source dataset with no / at either end",
    type = "character"
  )
  
  parser$add_argument(
    "-dsp", "--data_split",
    nargs = 3,
    default = c(7, 1, 2),
    help = "Enter the train:valid:test split as ints separated by a space",
    type = "integer"
  )
  
  parser$add_argument(
    "-ms", "--mask_seed",
    default = NULL,
    help = "Enter the seed used in the masking processes",
    type = "integer"
  )
  
  parser$add_argument(
    "-tm", "--tail_mask",
    nargs = 2,
    default = NULL,
    help = "Enter the tail masking parameters as floats separated by a space: [frac of cells masked, max frac of tail depth masked]",
    type = "double"
  )
  
  parser$add_argument(
    "-slm", "--single_level_mask",
    default = NULL,
    help = "Enter the fraction of cells in which a single random level was masked",
    type = "double"
  )
  
  parser$add_argument(
    "-bs", "--batch_size",
    default = 256,
    help = "Enter the batch size used in training",
    type = "integer"
  )
  
  parser$add_argument(
    "-ebs", "--eval_batch_size",
    default = 40960,
    help = "Enter the evaluation batch size suitable to your hardware, defaults to 40960",
    type = "integer"
  )
  
  parser$add_argument(
    "-opt", "--optimizer",
    default = "adam",
    help = "Enter the name of the optimizer used",
    type = "character"
  )
  
  parser$add_argument(
    "-lr", "--lr",
    default = NULL,
    help = "Enter the learning rate used",
    type = "double"
  )
  
  parser$add_argument(
    "-l1", "--l1",
    default = NULL,
    help = "Enter L1 reg strength used if any",
    type = "double"
  )
  
  parser$add_argument(
    "-l2", "--l2",
    default = 0.01,
    help = "Enter L2 reg strength used if any",
    type = "double"
  )
  
  parser$add_argument(
    "-dp", "--dropout",
    default = 0.1,
    help = "Enter dropout rate used if any",
    type = "double"
  )
  
  parser$add_argument(
    "-norm", "--normalization_override", 
    default = FALSE,
    help = "is the counts data lop1p normalized after scaling to 10k? defaults to False", 
    type = "logical"
  )
  
  parser$add_argument(
    "-em", "--embeddings",
    default = "shallow",
    help = "extract embeddings? defaults to 'shallow', other options: ['deep', 'both']",
    type = "character"
  )
  
  parser$add_argument(
    "-knn", "--knn_scores",
    default = FALSE,
    help = "specify if you want scores based on k nearest neighbours, defaults to False",
    type = "logical"
  )
  
  parser$add_argument(
    "-umap", "--umap_embeddings",
    default = TRUE,
    help = "specify if you want umap embeddings, defaults to True",
    type = "logical"
  )
  
  parser$add_argument(
    "-irl", "--if_refine_labels",
    default = TRUE,
    help = "whether to refine labels, defaults to True",
    type = "logical"
  )
  
  parser$add_argument(
    "-nnbrs", "--n_neighbors", 
    default = 30, 
    help = "n_neighbors param for umaps, defaults to Seurat default 30", 
    type = "integer"
  )
  
  parser$add_argument(
    "-nc", "--n_components", 
    default = 2, 
    help = "n_components param for umaps, defaults to Seurat default 2", 
    type = "integer"
  )
  
  parser$add_argument(
    "-me", "--metric", 
    default = "cosine", 
    help = "metric param for umaps, defaults to Seurat default 'cosine'", 
    type = "character"
  )
  
  parser$add_argument(
    "-mdt", "--min_dist", 
    default = 0.3, 
    help = "min_dist param for umaps, defaults to Seurat default 0.3", 
    type = "numeric"
  )
  
  parser$add_argument(
    "-ulr", "--umap_lr", 
    default = 1.0, 
    help = "learning_rate param for umaps, defaults to Seurat default 1.0", 
    type = "numeric"
  )
  
  parser$add_argument(
    "-useed", "--umap_seed", 
    default = 42, 
    help = "random_state param for reproducibility of umaps, defaults to Seurat default 42", 
    type = "integer"
  )
  
  parser$add_argument(
    "-sp", "--spread", 
    default = 1.0, 
    help = "spread param for umaps, defaults to Seurat default 1.0", 
    type = "numeric"
  )
  
  parser$add_argument(
    "-uv", "--umap_verbose", 
    default = TRUE, 
    help = "verbose param for umaps, defaults to TRUE", 
    type = "logical"
  )
  
  parser$add_argument(
    "-uin", "--umap_init", 
    default = "spectral", 
    help = "init param for umaps, defaults to 'spectral', the other option is 'random'", 
    type = "character"
  )
  
  parser$add_argument(
    "-objd", "--object_disk", 
    default = TRUE, 
    help = "do you want to write object to disk? default is TRUE", 
    type = "logical"
  )
  
  parser$add_argument(
    "-ofd", "--out_file_disk", 
    default = TRUE, 
    help = "do you want to write separate files to disk? default is TRUE", 
    type = "logical"
  )
  
  parser$add_argument(
    "-po", "--process_obj",
    default = TRUE,
    help = "whether to process the object with PrepLabel, defaults to TRUE",
    type = "logical"
  )
  
  parser$add_argument(
    "-ca", "--cutoff_abs",
    default = 5,
    help = "absolute cutoff for PrepLabel, defaults to 5",
    type = "integer"
  )
  
  parser$add_argument(
    "-cf", "--cutoff_frac",
    default = 0.001,
    help = "fractional cutoff for PrepLabel, defaults to 0.001",
    type = "numeric"
  )
  
  args <- parser$parse_args()
  return(args)
}


arg_parse_out_R <- function(args) {
  cat("Reading arguments... \n\n")
  
  query_filepath <- args$filepath
  feature_names_col <- args$feature_names_col
  source_data_dir <- args$source_data_dir
  features_txt <- args$features_txt
  
  if (args$mode %in% c("independent", "cumulative")) {
    split_mode <- args$mode
  } else {
    stop("Mode should either be 'independent' or 'cumulative'")
  }
  
  model <- args$model
  loss_name <- args$loss
  epochs <- args$epochs
  train_seed <- args$train_seed
  data_seed <- args$data_seed
  data_source <- args$data_source
  data_split <- args$data_split
  mask_seed <- args$mask_seed
  tm_frac <- args$tail_mask
  lm_frac <- args$single_level_mask
  save <- TRUE  
  batch_size <- args$batch_size
  eval_batch_size <- args$eval_batch_size
  optimizer_name <- args$optimizer
  lr <- args$lr
  l1 <- args$l1
  l2 <- args$l2
  dropout <- args$dropout
  normalization_override <- args$normalization_override
  
  embeddings_mode <- args$embeddings
  if_knn_scores <- args$knn_scores
  if_umap_embeddings <- args$umap_embeddings
  if_refine_labels <- args$if_refine_labels
  n_neighbors <- args$n_neighbors
  n_components <- args$n_components
  metric <- args$metric
  min_dist <- args$min_dist
  umap_lr <- args$umap_lr
  umap_seed <- args$umap_seed
  spread <- args$spread
  verbose <- args$umap_verbose
  init <- args$umap_init
  
  object_disk <- args$object_disk
  out_file_disk <- args$out_file_disk
  process_obj <- args$process_obj
  cutoff_abs <- args$cutoff_abs
  cutoff_frac <- args$cutoff_frac
  
  arguments <- list(
    query_filepath, 
    feature_names_col,
    source_data_dir,
    features_txt,
    split_mode,
    model,
    loss_name,
    epochs,
    train_seed,
    data_seed,
    data_source,
    data_split,
    mask_seed,
    tm_frac,
    lm_frac,
    save,
    batch_size,
    eval_batch_size,
    optimizer_name,
    lr,
    l1,
    l2,
    dropout,
    normalization_override,
    embeddings_mode,
    if_knn_scores,
    if_umap_embeddings,
    if_refine_labels,
    n_neighbors,
    n_components,
    metric,
    min_dist,
    umap_lr,
    umap_seed,
    spread,
    verbose,
    init,
    object_disk,
    out_file_disk,
    process_obj,
    cutoff_abs,
    cutoff_frac
  )
  
  return(arguments)
}


annotate_R <- function(){
  
  #### parsing arguments
  cat("\n")
  
  
  args <- arg_parse_in_R()
  
  #### reading parsed arguments
  arguments <- arg_parse_out_R(args)
  
  query_filepath <- arguments[[1]]
  feature_names_col <- arguments[[2]]
  source_data_dir <- arguments[[3]]
  features_txt <- arguments[[4]]
  split_mode <- arguments[[5]]
  model <- arguments[[6]]
  loss_name <- arguments[[7]]
  epochs <- arguments[[8]]
  train_seed <- arguments[[9]]
  data_seed <- arguments[[10]]
  data_source <- arguments[[11]]
  data_split <- arguments[[12]]
  mask_seed <- arguments[[13]]
  tm_frac <- arguments[[14]]
  lm_frac <- arguments[[15]]
  save <- arguments[[16]]
  batch_size <- arguments[[17]]
  eval_batch_size <- arguments[[18]]
  optimizer_name <- arguments[[19]]
  lr <- arguments[[20]]
  l1 <- arguments[[21]]
  l2 <- arguments[[22]]
  dropout <- arguments[[23]]
  normalization_override <- arguments[[24]]
  embeddings_mode <- arguments[[25]]
  if_knn_scores <- arguments[[26]]
  if_umap_embeddings <- arguments[[27]]
  if_refine_labels <- arguments[[28]]
  n_neighbors <- arguments[[29]]
  n_components <- arguments[[30]]
  metric <- arguments[[31]]
  min_dist <- arguments[[32]]
  umap_lr <- arguments[[33]]
  umap_seed <- arguments[[34]]
  spread <- arguments[[35]]
  verbose <- arguments[[36]]
  init <- arguments[[37]]
  object_disk <- arguments[[38]]
  out_file_disk <- FALSE
  process_obj <- arguments[[39]]
  cutoff_abs <- arguments[[40]]
  cutoff_frac <- arguments[[41]]
  
  #### reading the seurat object 
  query <- read_obj_R(query_filepath, feature_names_col)
  X_query <- sp$csr_matrix(r_to_py(query$X_query))
  query_features <- query$query_features
  query_cells_df <- query$query_cells_df
  query_obj <- query$query_obj
  
  
  #### reading the seurat object 
  query <- read_obj_min(query_obj, feature_names_col)
  X_query <- sp$csr_matrix(r_to_py(query$X_query))
  query_features <- query$query_features
  query_cells_df <- query$query_cells_df

  
  
  
  #### run annotate_core
  core_outputs <- annotate$annotate_core(
    X_query,
    query_features,
    source_data_dir,
    features_txt,
    split_mode,
    model,
    epochs,
    train_seed,
    loss_name,
    data_seed,
    data_source,
    data_split,
    mask_seed,
    tm_frac,
    lm_frac,
    batch_size,
    optimizer_name,
    lr,
    l1,
    l2, 
    dropout,
    save,
    eval_batch_size,
    normalization_override,
    embeddings_mode,
    query_cells_df,
    if_knn_scores,
    if_umap_embeddings,
    if_refine_labels,
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
  
  
  embeddings_mode <- core_outputs[[3]]
  embeddings_dict <- core_outputs[[4]]
  query_cells_df <- core_outputs[[10]]
  if_umap_embeddings <- core_outputs[[11]]
  umap_embeddings_dict <- core_outputs[[12]]
  
  annotated_obj = package_obj(embeddings_mode, embeddings_dict, if_umap_embeddings, umap_embeddings_dict, query_cells_df, query_obj)
  
  # this is a joke about massively dieting the object
  # so that we don't redownload expression data
  # in this API case we know the object has an RNA assay with counts and data layers
  keto_object=TRUE
  if (keto_object) {
    annotated_obj[["RNA"]]$data <- Matrix(0, nrow=nrow(annotated_obj[["RNA"]]$data),ncol=ncol(annotated_obj[["RNA"]]$data))
  }
  if (object_disk){
    dir_path <- dirname(query_filepath)
    file_name <- basename(query_filepath)
    
    file_name_ann <- sub("\\.rds$", "_ANN.rds", file_name)
    
    new_path <- file.path(dir_path, file_name_ann)
    
    message("Saving annotated object to ", new_path)
    
    saveRDS(annotated_obj, file = new_path)
  }
  
  return(annotated_obj)
  
}




args <- commandArgs(trailingOnly = TRUE)

if (length(args) > 0) {
  cat("\n")
  cat("ANNotate_R is executed from the terminal...")
  annotate_R() 
} else {
  cat("\n")
  cat("ANNotate_R is being run interactively.")
}
