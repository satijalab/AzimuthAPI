require('httr')
require('jsonlite')
require('curl')
library(Seurat)
# Check Seurat version
if (!requireNamespace("Seurat", quietly = TRUE)) {
  stop("Seurat is not installed. Please install it to proceed.")
} else {
  seurat_version <- as.numeric(substr(as.character(packageVersion("Seurat")),1,1))
  if (seurat_version < 5) {
    stop("Seurat version must be greater than 5. Please update Seurat.")
  }
}

# Function to handle progress updates from the API
listen_to_progress <- function(url, file_path) {
  # Create a multipart form for the upload
  form <- list(file = upload_file(file_path))
  
  # Open a connection to the API using curl
  handle <- new_handle()
  handle_setform(handle, .list = form)
  
  # Process the SSE stream in real time
  curl_fetch_stream(
    url,
    handle = handle,
    fun = function(data) {
      lines <- strsplit(rawToChar(data), "\n")[[1]]
      for (line in lines) {
        # SSE messages start with "data: "
        if (startsWith(line, "data: ")) {
          json_message <- sub("^data: ", "", line)
          tryCatch({
            parsed <- fromJSON(json_message)
            
            # Print messages based on content
            if ("message" %in% names(parsed)) {
              cat(parsed$message, "\n")
            }
            if ("output" %in% names(parsed)) {
              cat("Output:", parsed$output, "\n")
            }
            if ("error" %in% names(parsed)) {
              cat("Output:", parsed$error, "\n")
            }
          }, error = function(e) {
            # Handle JSON parsing errors gracefully
            #cat("Warning: Failed to parse message, skipping. Error:", e$message, "\n")
          })
        }
      }
    }
  )
}


# Function to download the output file
download_output <- function(output_url, save_as) {
  response <- GET(output_url)
  
  if (response$status_code == 200) {
    writeBin(content(response, "raw"), save_as)
    cat("Output file saved as:", save_as, "\n")
  } else {
    cat("Failed to download output file. Status code:", response$status_code, "\n")
  }
}

# Main function to upload and process RDS file
process_rds_file <- function(api_url, file_path) {
  # Step 1: Upload and listen for progress updates
  progress_url <- paste0(api_url, "/process_rds")
  cat("Uploading file and listening for updates...\n")
  listen_to_progress(progress_url, file_path)
  
  # Step 2: Download the output file
  output_file_name <- gsub("\\.rds$", "_ANN.rds", basename(file_path))
  download_url <- paste0(api_url, "/download_output?output_file=/tmp/", output_file_name)
  save_path <- file.path(tempdir(), output_file_name)
  
  cat("Downloading the output file...\n")
  download_output(download_url, save_path)
}


PrepLabel <- function(object, label_id = 'final_level_label', newid = 'PrepLabel', cutid = 'Other', cutoff=10) {
  rejected_names <- names(which(table(object@meta.data[,label_id])<cutoff))
  object@meta.data[,newid]=as.character(object@meta.data[,label_id])
  rejected_cells <- which(object@meta.data[,label_id]%in%rejected_names)
  object@meta.data[rejected_cells,newid]=cutid
  return(object)
}

CloudANNotate <- function(object = object, assay='RNA',ip='10.4.120.13',port=5000, normalize=FALSE) {
  
  # normalize data to be safe
  if (normalize) {
    object <- NormalizeData(object)
  }
  layer_name='data'
  data <- LayerData(object,assay = assay,layer = layer_name)
  feature_file <- '/brahms/shared/AzimuthAPI/22tissue_nsforestfeatures.txt'
  if (file.exists(feature_file)) {
    features <- readLines(feature_file)
    data <- data[intersect(features,rownames(data)),]
  }
  suppressWarnings(srt <- CreateSeuratObject(CreateAssay5Object(data=data)))  
  tmpname <- tempfile()
  tmp_input <- paste0(tmpname,".rds")
  tmp_output <- paste0(tmpname,"_ANN.rds")
  saveRDS(object = srt,file = tmp_input)
  api_base_url <- paste0('http://',ip,":",port)
  
  process_rds_file(api_base_url, tmp_input)
  srt <- readRDS(file = tmp_output)
  srt <- PrepLabel(srt,label_id = 'final_level_label',cutoff = max(5,0.001*ncol(srt)),cutid = 'Other',newid = 'azimuth_label')
  for(i in names(srt@reductions)) {
    object[[i]] <- srt[[i]]
  }
  
  # we want all metadata coluns except for the QC ones
  md_cols <- grep("nCount_RNA|nFeature_RNA",colnames(srt@meta.data),invert = TRUE,value = TRUE)
  for(i in md_cols) {
    object@meta.data[,i] <- srt@meta.data[,i]
  }
  Idents(object) <- 'azimuth_label'
  return(object)
}
