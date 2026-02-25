#' Run Pan-Human Azimuth annotation on the cloud
#'
#' @param object Seurat object to annotate
#' @param assay Name of the assay to use (default: 'RNA')
#' @param ip Hostname or IP address of the cloud server (default: 'azimuthapi.satijalab.org')
#' @param port Port number for the API (default: 5000)
#' @param ... Additional arguments for the API to pass to the model
#' @return Annotated Seurat object
#' @importFrom httr POST GET upload_file content status_code
#' @importFrom RCurl url.exists
#' @importFrom SeuratObject LayerData Idents IsMatrixEmpty CreateAssay5Object CreateSeuratObject Cells
#' @export
CloudAzimuth <- function(object = object, assay = 'RNA', ip = 'azimuthapi.satijalab.org', 
                         port = 5000, ...) {
  message("Running Pan-Human Azimuth on the cloud!")
  
  layer_name <- 'data'
  tryCatch({
    data <- LayerData(object, assay = assay, layer = layer_name)
  }, warning = function(w) {
    message <- conditionMessage(w)
    if (grepl(paste0("Layer ‘", layer_name, "’ is empty"), message)) {
      stop(simpleError(
          "Please run NormalizeData on the data before running Azimuth",
          call = conditionCall(w)
      ))
    }
  })

  # check if data has been normalized, just using the first 5 cells
  # throw an error if large values, or all integer values, are detected
  data_check <- data[,1:min(5,ncol(data))]
  if ((max(data_check) > 15) || isTRUE(all.equal(data_check,floor(data_check)))) {
    stop("Please run NormalizeData on the data before running Azimuth")
  }

  feature_file <- 'https://raw.githubusercontent.com/rsatija/public_utils/refs/heads/main/features_v0.txt'
  
  if (url.exists(feature_file)) {
    features <- readLines(feature_file)
    data <- data[intersect(features, rownames(data)), ]
  }
  
  message("Uploading dataset")
  suppressWarnings(srt <- CreateSeuratObject(CreateAssay5Object(data = data)))
  
  tmpname <- tempfile()
  tmp_input <- paste0(tmpname, ".rds")
  tmp_output <- paste0(tmpname, "_ANN.rds")
  saveRDS(object = srt, file = tmp_input)
  api_base_url <- paste0('http://', ip, ":", port)
  
  process_rds_file(api_base_url, tmp_input, ...)
  srt <- readRDS(file = tmp_output)
  
  # Copy reductions
  # Match assay.used slot to the assay processed
  for (i in names(srt@reductions)) {
    object[[i]] <- srt[[i]]
    object[[i]]@assay.used <- assay
  }
  
  # Copy metadata columns except QC ones, matching by cell names
  md_cols <- grep("nCount_RNA|nFeature_RNA", 
                 colnames(srt@meta.data), 
                 invert = TRUE, 
                 value = TRUE)
  
  # Only update metadata for cells that were processed (objects can have multiple assays w/ cell IDs formatted differently)
  assay_cells <- Cells(srt)
  
  for (i in md_cols) {
    # Initialize column with NA only if it doesn't exist
    if (!i %in% colnames(object@meta.data)) {
      object@meta.data[, i] <- NA
    }
    # Update only the processed cells
    object@meta.data[assay_cells, i] <- srt@meta.data[assay_cells, i]
  }
  
  Idents(object) <- 'azimuth_label'
  return(object)
}

#' Process RDS file through the cloud API
#'
#' @param api_base_url Base URL for the API
#' @param file_path Path to input RDS file
#' @param ... Additional arguments passed to the API
#' @return NULL
#' @importFrom httr POST GET upload_file content status_code
process_rds_file <- function(api_base_url, file_path, ...) {
  progress_url <- paste0(api_base_url, "/process_rds")
  cat("Uploading file and listening for updates...\n")
  success <- listen_to_progress(progress_url, file_path, ...)
  
  if (isFALSE(success)) {
    stop("Processing failed on the server. Please check the error messages above.")
  }
  
  output_file_name <- gsub("\\.rds$", "_ANN.rds", basename(file_path))
  download_url <- paste0(api_base_url, "/download_output?output_file=/tmp/", 
                         output_file_name)
  save_path <- file.path(tempdir(), output_file_name)
  cat("Downloading the output file...\n")
  download_output(download_url, save_path)
} 

