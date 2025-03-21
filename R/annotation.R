#' Run Pan-Human Azimuth annotation
#'
#' @param query_obj Seurat object to annotate
#' @param feature_names_col Column name for feature names
#' @param annotation_pipeline Set to 'supervised' as a default
#' @param eval_batch_size Evaluation batch size
#' @param normalization_override Whether to override normalization
#' @param norm_check_batch_size Batch size to inspect normalization of data
#' @param output_mode Output mode for annotated cell metadata
#' @param refine_labels Whether to refine labels
#' @param extract_embeddings Whether to azimuth embeddings
#' @param umap_embeddings Whether to include UMAP embeddings
#' @param n_neighbors Number of neighbors for UMAP
#' @param n_components Number of components for UMAP
#' @param metric Distance metric for UMAP
#' @param min_dist Minimum distance for UMAP
#' @param umap_lr Learning rate for UMAP
#' @param umap_seed Seed for UMAP
#' @param spread Spread parameter for UMAP
#' @param verbose Whether to show progress
#' @param init Initialization method for UMAP
#' @param process_obj Whether to process the object
#' @param cutoff_abs Absolute cutoff for label filtering
#' @param cutoff_frac Fractional cutoff for label filtering
#' @return Annotated Seurat object
#' @export
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
  #source_data_dir <- paste0(python_module_path, source_data_dir)
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