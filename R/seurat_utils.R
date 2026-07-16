#' @importFrom rlang %||% .data
#' @importFrom dplyr %>%
NULL


#' Read a Seurat object from an RDS file
#'
#' Reads a Seurat object from an RDS file, ensuring
#' compatibility with Seurat versions >= 4.4.0
#'
#' @param filepath Character string specifying the path to the RDS file
#' @param assay_name Character string specifying the assay to use (default: "RNA")
#' @importFrom utils packageVersion
#' @importFrom Seurat DefaultAssay<-
#'
#' @return A valid Seurat object
#'
#' @keywords internal
#' @noRd
#'
#' @examples
#' \dontrun{
#' seu_obj <- read_seurat_object("path/to/seurat_object.rds")
#' }
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
    warning(paste("Current Seurat version:", seurat_version, 
                  "is below 4.4.0. Some functionality may not work as expected."))
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

#' Save a Seurat object to an RDS file
#'
#' Saves a Seurat object to an RDS file after applying "keto diet" to reduce size
#'
#' @param seurat_obj A Seurat object to save
#' @param filepath Character string specifying the output path
#' @importFrom utils packageVersion
#' @importFrom methods slotNames slot slot<-
#'
#' @return Invisible filepath where the object was saved
#' @keywords internal
#' @noRd
#'
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
  # in this API case we know the object has an RNA assay with counts and data layers
  keto_object=TRUE
  if (keto_object) {
    if ("RNA" %in% names(seurat_obj@assays)) {
      # Handle Seurat v5+ and v4 differently
      if (packageVersion("Seurat") >= "5.0.0") {
        # For Seurat v5+
        if ("data" %in% names(seurat_obj[["RNA"]]@layers)) {
          seurat_obj[["RNA"]]$data <- Matrix::Matrix(0, 
                                                    nrow = nrow(seurat_obj[["RNA"]]$data),
                                                    ncol = ncol(seurat_obj[["RNA"]]$data),
                                                    sparse = TRUE)
        }
        if ("counts" %in% names(seurat_obj[["RNA"]]@layers)) {
          seurat_obj[["RNA"]]$counts <- Matrix::Matrix(0, 
                                                      nrow = nrow(seurat_obj[["RNA"]]$counts),
                                                      ncol = ncol(seurat_obj[["RNA"]]$counts),
                                                      sparse = TRUE)
        }
      } else {
        # For Seurat v4
        if ("data" %in% slotNames(seurat_obj[["RNA"]])) {
          slot(seurat_obj[["RNA"]], "data") <- Matrix::Matrix(0, 
                                                             nrow = nrow(slot(seurat_obj[["RNA"]], "data")),
                                                             ncol = ncol(slot(seurat_obj[["RNA"]], "data")),
                                                             sparse = TRUE)
        }
        if ("counts" %in% slotNames(seurat_obj[["RNA"]])) {
          slot(seurat_obj[["RNA"]], "counts") <- Matrix::Matrix(0, 
                                                               nrow = nrow(slot(seurat_obj[["RNA"]], "counts")),
                                                               ncol = ncol(slot(seurat_obj[["RNA"]], "counts")),
                                                               sparse = TRUE)
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

#' Get data from a Seurat object layer
#'
#' @param object Seurat object
#' @param assay Assay name
#' @param layer Layer name (default: "data")
#' @importFrom methods slot
#' @importFrom utils packageVersion
#' @importFrom SeuratObject LayerData
#' @return Layer data
#' @keywords internal
#' @noRd
get_data <- function(object, assay, layer = "data") {
  if (packageVersion("Seurat") >= "5.0.0") {
    return(LayerData(object[[assay]], layer = layer))
  } else {
    assay_obj <- object@assays[[assay]]
    return(slot(assay_obj, layer))
  }
}

#' Read and process a Seurat object
#'
#' @param query_obj Seurat object
#' @param feature_names_col Column name for feature names
#' @param assay_default Default assay name
#' @importFrom Seurat DefaultAssay<- NormalizeData
#' @importFrom SeuratObject Layers
#' @importFrom Matrix Matrix
#' @importFrom utils packageVersion
#' @importFrom methods slotNames
#' @return List containing processed data
#' @keywords internal
#' @noRd
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
  X_query <- Matrix::t(normalized_data)
  X_query <- Matrix(X_query, sparse = TRUE)
  assay_cells <- colnames(normalized_data)
  
  cell_metadata <- query_obj@meta.data
  query_cells_df <- as.data.frame(cell_metadata[assay_cells, , drop = FALSE])
  
  if (!is.null(feature_names_col)) {
    feature_metacols <- colnames(query_obj[[assay_default]][[]])
    if (feature_names_col %in% feature_metacols){
      query_features <- query_obj[[assay_default]][[feature_names_col]]
      query_features <- as.list(query_features[[feature_names_col]])
    } else {
      stop(paste(
        feature_names_col, 
        "not found as a column in the df returned by ",
        "object[[",assay_default,"]][[]]"
        ))
    }
    
  } else {
    query_features <- as.list(rownames(normalized_data))
  }
  return(list(
    X_query = X_query, 
    query_features = query_features, 
    query_cells_df = query_cells_df,
    assay_cells = assay_cells
    ))
}



#' Package embeddings into a Seurat object
#'
#' @param extract_embeddings Whether embeddings are extracted 
#' @param embeddings_dict Dictionary of embeddings
#' @param umap_embeddings Whether UMAP embeddings are computed
#' @param umap_embeddings_dict Dictionary of UMAP embeddings
#' @param query_cells_df Cell metadata
#' @param query_cells Cells present in the processed assay
#' @param query_obj Seurat object
#' @importFrom SeuratObject CreateDimReducObject
#' @importFrom Seurat DefaultAssay
#' @return Updated Seurat object
#' @keywords internal
#' @noRd
package_obj <- function(
  extract_embeddings, 
  embeddings_dict, 
  umap_embeddings, 
  umap_embeddings_dict, 
  query_cells_df, 
  query_cells,
  query_obj
  ) {
  if (nrow(query_cells_df) != length(query_cells)) {
    stop("Dimension mismatch: query cells does not match the processed assay cells.")
  }

  if (!identical(rownames(query_cells_df), query_cells)) {
    if (setequal(rownames(query_cells_df), query_cells)) {
      query_cells_df <- query_cells_df[query_cells, , drop = FALSE]
    } else {
      rownames(query_cells_df) <- query_cells
    }
  }
  
  if (extract_embeddings) {
    for (em_name in names(embeddings_dict)) {
      
      em_matrix <- as.matrix(embeddings_dict[[em_name]])
      
      if (nrow(em_matrix) != length(query_cells)) {
        stop(paste(
          "Dimension mismatch:", em_name, " does not have as many ",
          "cells as the processed assay."
          ))
      }
      
      rownames(em_matrix) <- query_cells
      
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
      
      if (nrow(em_matrix) != length(query_cells)) {
        stop(paste(
          "Dimension mismatch:", em_name, " does not ",
          "have as many cells as the processed assay."
          ))
      }
      
      rownames(em_matrix) <- query_cells
      
      dimreduc_obj <- CreateDimReducObject(
        embeddings = em_matrix, 
        key = paste0(em_name, "_"), 
        assay = DefaultAssay(query_obj)
        )
      
      query_obj[[paste0("umapANN", em_name)]] <- dimreduc_obj
    }
  }
  
  for (md_col in colnames(query_cells_df)) {
    if (!md_col %in% colnames(query_obj@meta.data)) {
      query_obj@meta.data[, md_col] <- NA
    }
    query_obj@meta.data[query_cells, md_col] <- query_cells_df[query_cells, md_col]
  }
  
  return(query_obj)
}