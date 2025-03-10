#!/usr/bin/env Rscript 


#### creating and setting up conda environment
if (!requireNamespace("reticulate", quietly = TRUE)) {
  install.packages("reticulate")
}
library(reticulate)
reticulate::use_condaenv("AzimuthNN_min", required=TRUE)
#### R dependencies
library(Seurat)
library(Matrix)
library(argparse)

#### other python dependencies
python_module_path <- "/brahms/sarkars/AzimuthNN_clone/AzimuthNN/sarkars"
py_run_string(paste("import sys; sys.path.append('", python_module_path, "')", sep = ""))
annotate <- import("ANNotate")
sp <- import("scipy.sparse")

#### functions

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

read_obj_min <- function(query_obj, feature_names_col) {
  
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
    "-m", "--mode",
    default = "cumulative",
    help = "Enter the label mode: cumulative or independent",
    type = "character"
  )
  
  parser$add_argument(
    "-md", "--model",
    default = "M0.1_cumul",
    help = "Enter the model name",
    type = "character"
  )
  
  parser$add_argument(
    "-fn", "--feature_names_col",
    default = NULL,
    help = "Enter the column name where the feature names are stored in query.var where query is the anndata object read from the h5ad.",
    type = "character"
  )
  
  parser$add_argument(
    "--epochs",
    default = 25,
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
    "-l", "--loss",
    default = "level_wt_focal_loss",
    help = "Enter the loss function used for optimization",
    type = "character"
  )
  
  parser$add_argument(
    "-ds", "--data_seed",
    default = 101,
    help = "Enter the data prep seed that was used",
    type = "integer"
  )
  
  parser$add_argument(
    "-dso", "--data_source",
    default = "data/kfold_data/datasets/fold10_11_28_2024_23_10_101",
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
    "-l2", "--l2",
    default = NULL,
    help = "Enter L2 reg strength used if any",
    type = "double"
  )
  
  parser$add_argument(
    "-dp", "--dropout",
    default = NULL,
    help = "Enter dropout rate used if any",
    type = "double"
  )
  
  parser$add_argument(
    "-norm", "--normalization_override", 
    default=FALSE,
    help="is the counts data lop1p normalized after scaling to 10k? defaults to False", 
    type="logical"
  )
  
  parser$add_argument(
    "-em", "--embeddings",
    default=NULL,
    help="extract embeddings? defaults to None, other options: ['shallow', 'deep', 'both']",
    type="character"
  )
  
  parser$add_argument(
    "-knn", "--knn_scores",
    default=FALSE,
    help="specify if you want scores based on k nearest neighbours, defaults to False",
    type="logical"
  )
  
  parser$add_argument(
    "-umap", "--umap_embeddings",
    default=FALSE,
    help="specify if you want umap embeddings, defaults to False",
    type="logical"
  )
  
  parser$add_argument(
    "-nnbrs", "--n_neighbors", 
    default = 30, 
    help = "n_neighbors param for umaps, defaults to Seurat default 30", 
    type = "integer")
  
  parser$add_argument(
    "-nc", "--n_components", 
    default = 2, 
    help = "n_components param for umaps, defaults to Seurat default 2", 
    type = "integer")
  
  parser$add_argument(
    "-me", "--metric", 
    default = "cosine", 
    help = "metric param for umaps, defaults to Seurat default 'cosine'", 
    type = "character")
  
  parser$add_argument(
    "-mdt", "--min_dist", 
    default = 0.3, 
    help = "min_dist param for umaps, defaults to Seurat default 0.3", 
    type = "numeric")
  
  parser$add_argument(
    "-ulr", "--umap_lr", 
    default = 1.0, 
    help = "learning_rate param for umaps, defaults to Seurat default 1.0", 
    type = "numeric")
  
  parser$add_argument(
    "-useed", "--umap_seed", 
    default = 42, 
    help = "random_state param for reproducibility of umaps, defaults to Seurat default 42", 
    type = "integer")
  
  parser$add_argument(
    "-sp", "--spread", 
    default = 1.0, 
    help = "spread param for umaps, defaults to Seurat default 1.0", 
    type = "numeric")
  
  parser$add_argument(
    "-uv", "--umap_verbose", 
    default = TRUE, 
    help = "verbose param for umaps, defaults to TRUE", 
    type = "logical")
  
  parser$add_argument(
    "-uin", "--umap_init", 
    default = "spectral", 
    help = "init param for umaps, defaults to 'spectral', the other option is 'random'", 
    type = "character")
  
  parser$add_argument(
    "-objd", "--object_disk", 
    default = TRUE, 
    help = "do you want to write object to disk? default is TRUE", 
    type = "logical")
  
  parser$add_argument(
    "-ofd", "--out_file_disk", 
    default = TRUE, 
    help = "do you want to write separate files to disk? default is TRUE", 
    type = "logical")
  
  args <- parser$parse_args()
  return(args)
}

arg_parse_out_R <- function(args) {
  cat("Reading arguments... \n\n")
  
  query_filepath <- args$filepath
  feature_names_col <- args$feature_names_col
  
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
  l2 <- args$l2
  dropout <- args$dropout
  normalization_override <- args$normalization_override
  
  embeddings_mode <- args$embeddings
  if_knn_scores <- args$knn_scores
  if_umap_embeddings <- args$umap_embeddings
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
  
  arguments <- list(
    query_filepath, 
    feature_names_col, 
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
    l2,
    dropout,
    normalization_override,
    embeddings_mode,
    if_knn_scores,
    if_umap_embeddings,
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
    out_file_disk
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
  split_mode <- arguments[[3]]
  model <- arguments[[4]]
  loss_name <- arguments[[5]]
  epochs <- arguments[[6]]
  train_seed <- arguments[[7]]
  data_seed <- arguments[[8]]
  data_source <- arguments[[9]]
  data_split <- arguments[[10]]
  mask_seed <- arguments[[11]]
  tm_frac <- arguments[[12]]
  lm_frac <- arguments[[13]]
  save <- arguments[[14]]
  batch_size <- arguments[[15]]
  eval_batch_size <- arguments[[16]]
  optimizer_name <- arguments[[17]]
  lr <- arguments[[18]]
  l2 <- arguments[[19]]
  dropout <- arguments[[20]]
  normalization_override <- arguments[[21]]
  embeddings_mode <- arguments[[22]]
  if_knn_scores <- arguments[[23]]
  if_umap_embeddings <- arguments[[24]]
  n_neighbors <- arguments[[25]]
  n_components <- arguments[[26]]
  metric <- arguments[[27]]
  min_dist <- arguments[[28]]
  umap_lr <- arguments[[29]]
  umap_seed <- arguments[[30]]
  spread <- arguments[[31]]
  verbose <- arguments[[32]]
  init <- arguments[[33]]
  object_disk <- arguments[[34]]
  out_file_disk <- FALSE
  
  #### reading the seurat object 
  query <- read_obj_R(query_filepath, feature_names_col)
  X_query <- sp$csr_matrix(r_to_py(query$X_query))
  query_features <- query$query_features
  query_cells_df <- query$query_cells_df
  query_obj <- query$query_obj
  
  
  #### run annotate_core
  core_outputs <- annotate$annotate_core(
    X_query,
    query_features,
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
    l2, 
    dropout,
    save,
    eval_batch_size,
    normalization_override,
    embeddings_mode,
    query_cells_df,
    if_knn_scores,
    if_umap_embeddings,
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
  query_cells_df <- core_outputs[[9]]
  if_umap_embeddings <- core_outputs[[10]]
  umap_embeddings_dict <- core_outputs[[11]]
  
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



ANNotate <- function(query_obj,
                     feature_names_col=NULL,
                     split_mode="cumulative",
                     model="M0.1_cumul",
                     loss_name="level_wt_focal_loss",
                     epochs=25,
                     train_seed=100,
                     data_seed=101,
                     data_source="data/kfold_data/datasets/fold10_11_28_2024_23_10_101",
                     data_split=c(7,1,2),
                     mask_seed=NULL,
                     tm_frac=NULL,
                     lm_frac=NULL,
                     save=TRUE,
                     batch_size=256,
                     eval_batch_size=40960,
                     optimizer_name="adam",
                     lr=NULL,
                     l2=NULL,
                     dropout=0.2,
                     normalization_override=FALSE,
                     embeddings_mode=NULL,
                     if_knn_scores=FALSE,
                     if_umap_embeddings=FALSE,
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
                     out_file_disk=FALSE){
  
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
    l2, 
    dropout,
    save,
    eval_batch_size,
    normalization_override,
    embeddings_mode,
    query_cells_df,
    if_knn_scores,
    if_umap_embeddings,
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
  query_cells_df <- core_outputs[[9]]
  if_umap_embeddings <- core_outputs[[10]]
  umap_embeddings_dict <- core_outputs[[11]]
  
  annotated_obj = package_obj(embeddings_mode, embeddings_dict, if_umap_embeddings, umap_embeddings_dict, query_cells_df, query_obj)
  
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
