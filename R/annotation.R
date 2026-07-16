#' Run Pan-human Azimuth annotation locally
#' 
#' This function runs the Pan-human Azimuth model on a Seurat object to annotate cell types, via reticulate and the `panhumanpy` Python package. **We recommend using the `CloudAzimuth` function, which runs cloud-based annotation, can handle large datasets, and performs robust error handling.**
#'
#' This function requires the `panhumanpy` Python package to be installed and accessible via reticulate.
#'
#' @param query_obj Seurat object to annotate
#' @param feature_names_col Column name for feature names
#' @param annotation_pipeline Set to 'supervised' as a default
#' @param eval_batch_size Evaluation batch size
#' @param normalization_override Whether to override normalization
#' @param norm_check_batch_size Batch size to inspect normalization of data
#' @param output_mode Output mode for annotated cell metadata
#' @param refine_labels Whether to refine labels
#' @param map_to_cl One or more annotation columns to map to Cell Ontology labels
#' @param include_cl_id Whether to add Cell Ontology IDs to the output metadata
#' @param extract_embeddings Whether to extract Azimuth embeddings
#' @param umap_embeddings Whether to include UMAP of Azimuth embeddings
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
#' @param model_version Version of the model to use
#' @param assay Assay to use for annotation
#' 
#' @importFrom SeuratObject Idents<-
#' @importFrom reticulate r_to_py
#' @importFrom rlang %||%
#'
#' @concept annotation 
#' @return Annotated Seurat object
#' 
#' @export
ANNotate <- function(
  query_obj,
  feature_names_col = NULL,
  annotation_pipeline = 'supervised',
  eval_batch_size = 40960L,
  normalization_override = FALSE,
  norm_check_batch_size = 1000L,
  output_mode = 'minimal',
  refine_labels = TRUE,
  map_to_cl = NULL,
  include_cl_id = FALSE,
  extract_embeddings = TRUE,
  umap_embeddings = TRUE,
  n_neighbors = 30L,
  n_components = 2L,
  metric = "cosine",
  min_dist = 0.3,
  umap_lr = 1.0,
  umap_seed = 42L,
  spread = 1.0,
  verbose = TRUE,
  init = "spectral",
  model_version = "v0",
  process_obj = TRUE,
  cutoff_abs = 5,
  cutoff_frac = 0.001,
  assay = NULL
) {
  options(warn = -1)
  # python dependencies
  annotate <- reticulate::import("panhumanpy.ANNotate")
  sp <- reticulate::import("scipy.sparse")

  cat("Running Pan-Human Azimuth:\n")
  cat("\n")
  
  # Convert integers
  eval_batch_size <- as.integer(eval_batch_size)
  n_neighbors <- as.integer(n_neighbors)
  n_components <- as.integer(n_components)
  umap_seed <- as.integer(umap_seed)

  assay <- assay %||% DefaultAssay(query_obj)
  
  # Read and process the Seurat object
  query <- read_obj_min(query_obj, feature_names_col, assay_default = assay)
  X_query <- sp$csr_matrix(reticulate::r_to_py(query$X_query))
  query_features <- query$query_features
  cells_meta <- query$query_cells_df
  assay_cells <- query$assay_cells

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
    map_to_cl,
    include_cl_id,
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
    init,
    model_version
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
    assay_cells,
    query_obj
  )
  
  if (process_obj){
    annotated_obj <- PrepLabel(annotated_obj,
                              label_id = 'final_level_labels',
                              cutoff = min(cutoff_abs, cutoff_frac * ncol(annotated_obj)),
                              cutid = 'Other',
                              newid = 'azimuth_label')
    Idents(annotated_obj) <- 'azimuth_label'
  }
  
  return(annotated_obj)
} 

#' Prepare labels for annotation
#'
#' @param object Seurat object
#' @param label_id Column name for labels
#' @param newid New column name for processed labels
#' @param cutid Label for rejected cells
#' @param cutoff Minimum count threshold
#' @return Updated Seurat object
#' @concept annotation 
#' @export
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