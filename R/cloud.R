#' Run Pan-human Azimuth annotation on the cloud
#'
#' @param object Seurat object to annotate
#' @param assay Name of the assay to use (default: 'RNA')
#' @param ip Hostname or IP address of the cloud server (default: 'azimuthapi.satijalab.org')
#' @param port Port number for the API (default: 5000)
#' @param ... Additional arguments for the API to pass to the model (see ANNotate function for details)
#' @return Annotated Seurat object
#' @importFrom httr POST GET upload_file content status_code
#' @importFrom RCurl url.exists
#' @importFrom SeuratObject LayerData Idents IsMatrixEmpty CreateAssay5Object CreateSeuratObject Cells Idents<-
#' @concept annotation
#' @export
CloudAzimuth <- function(object = object, assay = 'RNA', ip = 'azimuthapi.satijalab.org', 
                         port = 5000, ...) {

  cli::cli_h1("Running Pan-human Azimuth on the cloud")

  api_base_url <- paste0('http://', ip, ":", port)
  update <- check_api_version(api_base_url)

  if (isTRUE(update)) {
    cli::cli_alert_warning("A new version of the AzimuthAPI package is available. Please update to the latest version for the best experience.")
  }

  layer_name <- 'data'
  tryCatch({
    data <- LayerData(object, assay = assay, layer = layer_name)
  }, warning = function(w) {
    warn_msg <- conditionMessage(w)
    if (grepl(paste0("Layer.*", layer_name, ".*is empty"), warn_msg)) {
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

  suppressWarnings(srt <- CreateSeuratObject(CreateAssay5Object(data = data)))

  tmpname <- tempfile()
  tmp_input <- paste0(tmpname, ".rds")
  tmp_output <- paste0(tmpname, "_ANN.rds")
  saveRDS(object = srt, file = tmp_input)

  cli::cli_alert_info("Uploading dataset...")

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
    # Convert to character if it's a factor to avoid level mismatch issues
    if (is.factor(object@meta.data[[i]])) {
      object@meta.data[[i]] <- as.character(object@meta.data[[i]])
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
#' @noRd
process_rds_file <- function(api_base_url, file_path, ...) {
  progress_url <- paste0(api_base_url, "/process_rds")
  cli::cli_alert_info("Uploading file and listening for updates...")
  result <- listen_to_progress(progress_url, file_path, ...)
  
  if (isFALSE(result$success)) {
    stop("Processing failed on the server. Please check the error messages above.")
  }

  output_file_name <- gsub("\\.rds$", "_ANN.rds", basename(file_path))
  if (!is.null(result$output_file)) {
    output_file_name <- basename(result$output_file)
  }
  download_url <- result$download_url
  if (is.null(download_url) && !is.null(result$output_file)) {
    download_url <- paste0(api_base_url, "/download_output?output_file=",
                            utils::URLencode(result$output_file, reserved = TRUE)
    )
  }
  if (is.null(download_url)) {
    stop("Server did not return output file metadata for download.")
  }
  if (!grepl("^https?://", download_url)) {
    download_url <- paste0(api_base_url, download_url)
  }
  save_path <- file.path(tempdir(), output_file_name)
  cli::cli_alert_info("Downloading the output file...")
  download_output(download_url, save_path)
  cli::cli_alert_success("Annotation complete. Output saved to: {save_path}")
} 

