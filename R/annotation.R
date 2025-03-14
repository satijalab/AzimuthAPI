#' Run Pan-Human Azimuth annotation
#'
#' @param query_obj Seurat object to annotate
#' @param feature_names_col Column name for feature names
#' @param source_data_dir Source data directory
#' @param features_txt Features text file
#' @param split_mode Split mode for annotation
#' @param model Model name
#' @param loss_name Loss function name
#' @param epochs Number of epochs
#' @param train_seed Training seed
#' @param data_seed Data seed
#' @param data_source Data source
#' @param data_split Train/validation/test split
#' @param mask_seed Mask seed
#' @param tm_frac Training mask fraction
#' @param lm_frac Label mask fraction
#' @param save Whether to save results
#' @param batch_size Batch size
#' @param eval_batch_size Evaluation batch size
#' @param optimizer_name Optimizer name
#' @param lr Learning rate
#' @param l1 L1 regularization
#' @param l2 L2 regularization
#' @param dropout Dropout rate
#' @param normalization_override Whether to override normalization
#' @param embeddings_mode Embeddings mode
#' @param if_knn_scores Whether to include KNN scores
#' @param if_umap_embeddings Whether to include UMAP embeddings
#' @param if_refine_labels Whether to refine labels
#' @param n_neighbors Number of neighbors for UMAP
#' @param n_components Number of components for UMAP
#' @param metric Distance metric for UMAP
#' @param min_dist Minimum distance for UMAP
#' @param umap_lr Learning rate for UMAP
#' @param umap_seed Seed for UMAP
#' @param spread Spread parameter for UMAP
#' @param verbose Whether to show progress
#' @param init Initialization method for UMAP
#' @param object_disk Whether to save object to disk
#' @param out_file_disk Whether to save output to disk
#' @param process_obj Whether to process the object
#' @param cutoff_abs Absolute cutoff for label filtering
#' @param cutoff_frac Fractional cutoff for label filtering
#' @return Annotated Seurat object
#' @export
ANNotate <- function(
  query_obj,
  feature_names_col = NULL,
  source_data_dir = "data/kfold_data",
  features_txt = "features.txt",
  split_mode = "cumulative",
  model = "M0.2",
  loss_name = "level_wt_focal_loss",
  epochs = 55,
  train_seed = 100,
  data_seed = 414,
  data_source = "data/kfold_data",
  data_split = c(7, 1, 2),
  mask_seed = NULL,
  tm_frac = NULL,
  lm_frac = NULL,
  save = TRUE,
  batch_size = 256,
  eval_batch_size = 40960,
  optimizer_name = "adam",
  lr = NULL,
  l1 = NULL,
  l2 = 0.01,
  dropout = 0.1,
  normalization_override = FALSE,
  embeddings_mode = "shallow",
  if_knn_scores = FALSE,
  if_umap_embeddings = TRUE,
  if_refine_labels = TRUE,
  n_neighbors = 30,
  n_components = 2,
  metric = "cosine",
  min_dist = 0.3,
  umap_lr = 1.0,
  umap_seed = 42,
  spread = 1.0,
  verbose = TRUE,
  init = "spectral",
  object_disk = FALSE,
  out_file_disk = FALSE,
  process_obj = TRUE,
  cutoff_abs = 5,
  cutoff_frac = 0.001
) {
  options(warn = -1)
  source_data_dir <- paste0(python_module_path, source_data_dir)
  cat("Running Pan-Human Azimuth:\n")
  cat("\n")
  
  # Convert integers
  epochs <- as.integer(epochs)
  train_seed <- as.integer(train_seed)
  data_seed <- as.integer(data_seed)
  data_split <- lapply(data_split, as.integer)
  if (!is.null(mask_seed)) {
    mask_seed <- as.integer(mask_seed)
  }
  batch_size <- as.integer(batch_size)
  eval_batch_size <- as.integer(eval_batch_size)
  n_neighbors <- as.integer(n_neighbors)
  n_components <- as.integer(n_components)
  umap_seed <- as.integer(umap_seed)
  
  # Read and process the Seurat object
  query <- read_obj_min(query_obj, feature_names_col)
  X_query <- sp$csr_matrix(r_to_py(query$X_query))
  query_features <- query$query_features
  query_cells_df <- query$query_cells_df
  
  # Run annotation core
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
  
  annotated_obj <- package_obj(embeddings_mode, embeddings_dict, if_umap_embeddings, 
                             umap_embeddings_dict, query_cells_df, query_obj)
  
  if (process_obj) {
    annotated_obj <- PrepLabel(annotated_obj, 
                             label_id = 'final_level_label',
                             cutoff = min(cutoff_abs, cutoff_frac * ncol(annotated_obj)),
                             cutid = 'Other',
                             newid = 'azimuth_label')
    Idents(annotated_obj) <- 'azimuth_label'
  }
  
  return(annotated_obj)
} 