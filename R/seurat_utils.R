#' Get data from a Seurat object layer
#'
#' @param object Seurat object
#' @param assay Assay name
#' @param layer Layer name (default: "data")
#' @return Layer data
#' @export
get_data <- function(object, assay, layer = "data") {
  if (packageVersion("Seurat") >= "5.0.0") {
    return(LayerData(object[[assay]], layer = layer))
  } else {
    assay_obj <- object@assays[[assay]]
    return(slot(assay_obj, layer))
  }
}

#' Read and process a Seurat object from file
#'
#' @param query_filepath Path to the RDS file
#' @param feature_names_col Column name for feature names
#' @return List containing processed data
#' @export
read_obj_R <- function(query_filepath, feature_names_col) {
  query_obj <- readRDS(query_filepath)
  if ("data" %in% names(query_obj@assays$RNA)) {
    normalized_data <- LayerData(query_obj, layer = 'data')
  } else {
    query_obj <- NormalizeData(query_obj)
    normalized_data <- LayerData(query_obj, layer = 'data')
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
  
  return(list(X_query = X_query, query_features = query_features, 
              query_cells_df = query_cells_df, query_obj = query_obj))
}

#' Read and process a Seurat object
#'
#' @param query_obj Seurat object
#' @param feature_names_col Column name for feature names
#' @param assay_default Default assay name
#' @return List containing processed data
#' @export
read_obj_min <- function(query_obj, feature_names_col, assay_default = 'RNA') {
  if (!(assay_default %in% names(query_obj@assays))) {
    stop(paste("Assay", assay_default, "not present in the seurat object."))
  }
  
  DefaultAssay(query_obj) <- assay_default
  
  if (packageVersion("Seurat") >= "5.0.0") {
    layers_obj <- Layers(query_obj, assay = assay_default)
  } else {
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
    if (feature_names_col %in% feature_metacols) {
      query_features <- query_obj[[assay_default]][[feature_names_col]]
      query_features <- as.list(query_features[[feature_names_col]])
    } else {
      stop(paste(feature_names_col, "not found as a column in the df returned by object[[", assay_default, "]][[]]"))
    }
  } else {
    query_features <- as.list(rownames(normalized_data))
  }
  
  return(list(X_query = X_query, query_features = query_features, query_cells_df = query_cells_df))
}

#' Package embeddings into a Seurat object
#'
#' @param embeddings_mode Mode of embeddings
#' @param embeddings_dict Dictionary of embeddings
#' @param if_umap_embeddings Whether to include UMAP embeddings
#' @param umap_embeddings_dict Dictionary of UMAP embeddings
#' @param query_cells_df Cell metadata
#' @param query_obj Seurat object
#' @return Updated Seurat object
#' @export
package_obj <- function(embeddings_mode, embeddings_dict, if_umap_embeddings, 
                       umap_embeddings_dict, query_cells_df, query_obj) {
  if (!is.null(embeddings_mode)) {
    for (em_name in names(embeddings_dict)) {
      em_matrix <- as.matrix(embeddings_dict[[em_name]])
      
      if (nrow(em_matrix) != length(Cells(query_obj))) {
        stop(paste("Dimension mismatch:", em_name, " does not have as many cells as the query obj."))
      }
      rownames(em_matrix) <- Cells(query_obj)
      dimreduc_obj <- CreateDimReducObject(embeddings = em_matrix, 
                                         key = paste0("ANN", em_name, "_"), 
                                         assay = DefaultAssay(query_obj))
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
      dimreduc_obj <- CreateDimReducObject(embeddings = em_matrix, 
                                         key = paste0("umapANN", em_name, "_"), 
                                         assay = DefaultAssay(query_obj))
      query_obj[[paste0("umapANN", em_name)]] <- dimreduc_obj
    }
  }
  
  query_obj@meta.data <- as.data.frame(query_cells_df)  
  return(query_obj)
}

#' Prepare labels for annotation
#'
#' @param object Seurat object
#' @param label_id Column name for labels
#' @param newid New column name for processed labels
#' @param cutid Label for rejected cells
#' @param cutoff Minimum count threshold
#' @return Updated Seurat object
#' @export
PrepLabel <- function(object, label_id = 'final_level_label', newid = 'PrepLabel', 
                     cutid = 'Other', cutoff = 10) {
  rejected_names <- names(which(table(object@meta.data[, label_id]) < cutoff))
  object@meta.data[, newid] <- as.character(object@meta.data[, label_id])
  rejected_cells <- which(object@meta.data[, label_id] %in% rejected_names)
  object@meta.data[rejected_cells, newid] <- cutid
  return(object)
} 