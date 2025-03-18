#' Read a Seurat object from an RDS file
#'
#' Reads a Seurat object from an RDS file, ensuring
#' compatibility with Seurat versions >= 4.4.0
#'
#' @param filepath Character string specifying the path to the RDS file
#' @param assay_name Character string specifying the assay to use (default: "RNA")
#'
#' @return A valid Seurat object
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