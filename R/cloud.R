#' Run Pan-Human Azimuth annotation on the cloud
#'
#' @param object Seurat object to annotate
#' @param assay Name of the assay to use (default: 'RNA')
#' @param ip IP address of the cloud server (default: '34.222.135.233')
#' @param port Port number for the API (default: 5000)
#' @return Annotated Seurat object
#' @importFrom httr POST GET upload_file content status_code
#' @importFrom RCurl url.exists
#' @export
CloudAzimuth <- function(object = object, assay = 'RNA', ip = '34.222.135.233', 
                         port = 5000) {
  message("Running Pan-Human Azimuth on the cloud!")
  
  layer_name <- 'data'
  data <- LayerData(object, assay = assay, layer = layer_name)

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
  
  process_rds_file(api_base_url, tmp_input)
  srt <- readRDS(file = tmp_output)
  
  # Copy reductions
  for (i in names(srt@reductions)) {
    object[[i]] <- srt[[i]]
  }
  
  # Copy metadata columns except QC ones
  md_cols <- grep("nCount_RNA|nFeature_RNA", 
                 colnames(srt@meta.data), 
                 invert = TRUE, 
                 value = TRUE)
  for (i in md_cols) {
    object@meta.data[, i] <- srt@meta.data[, i]
  }
  
  Idents(object) <- 'azimuth_label'
  return(object)
}

#' Process RDS file through the cloud API
#'
#' @param api_base_url Base URL for the API
#' @param input_file Path to input RDS file
#' @return NULL
#' @importFrom httr POST GET upload_file content status_code
process_rds_file <- function(api_base_url, file_path) {
  progress_url <- paste0(api_base_url, "/process_rds")
  cat("Uploading file and listening for updates...\n")
  listen_to_progress(progress_url, file_path)
  output_file_name <- gsub("\\.rds$", "_ANN.rds", basename(file_path))
  download_url <- paste0(api_base_url, "/download_output?output_file=/tmp/", 
                         output_file_name)
  save_path <- file.path(tempdir(), output_file_name)
  cat("Downloading the output file...\n")
  download_output(download_url, save_path)
} 