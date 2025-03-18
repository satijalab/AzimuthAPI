#' Save a Seurat object to an RDS file
#'
#' Saves a Seurat object to an RDS file after applying "keto diet" to reduce size
#'
#' @param seurat_obj A Seurat object to save
#' @param filepath Character string specifying the output path
#'
#' @return Invisible filepath where the object was saved
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