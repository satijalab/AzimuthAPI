#' Listen to progress updates from the API
#'
#' @param url API endpoint URL
#' @param file_path Path to the file being processed
#' @return NULL
#' @importFrom curl new_handle handle_setform curl_fetch_stream
#' @importFrom jsonlite fromJSON
#' @export
listen_to_progress <- function(url, file_path) {
  # Create a multipart form for the upload
  form <- list(file = upload_file(file_path))
  
  # Open a connection to the API using curl
  handle <- new_handle()
  handle_setform(handle, .list = form)
  
  # Buffer to store partial messages
  buffer <- ""
  
  # Process the SSE stream in real time
  curl_fetch_stream(
    url,
    handle = handle,
    fun = function(data) {
      # Convert raw data to character and append to buffer
      buffer <<- paste0(buffer, rawToChar(data))
      
      # Split buffer into lines
      lines <- strsplit(buffer, "\n")[[1]]
      
      # Keep the last incomplete line in buffer
      if (length(lines) > 0) {
        buffer <<- lines[length(lines)]
        lines <- lines[-length(lines)]
      }
      
      for (line in lines) {
        # Skip empty lines and ":" prefixed lines
        if (nchar(line) == 0 || startsWith(line, ":")) {
          next
        }
        
        # SSE messages start with "data: "
        if (startsWith(line, "data: ")) {
          json_message <- sub("^data: ", "", line)
          tryCatch({
            parsed <- fromJSON(json_message)
            
            # Print messages based on content
            if (!is.null(parsed$message)) {
              cat(parsed$message, "\n")
            }
            if (!is.null(parsed$output)) {
              cat(parsed$output, "\n")
            }
            if (!is.null(parsed$error)) {
              cat(parsed$error, "\n")
            }
            if (!is.null(parsed$progress)) {
              cat("Progress:", parsed$progress, "%\n")
            }
            if (!is.null(parsed$status)) {
              cat("Status:", parsed$status, "\n")
            }
          }, error = function(e) {
            # Only print actual parsing errors, not empty or malformed messages
            if (!grepl("unexpected end", e$message)) {
              cat("Warning: Failed to parse message:", e$message, "\n")
              cat("Raw message:", json_message, "\n")
            }
          })
        }
      }
    }
  )
}

#' Download output file from the API
#'
#' @param url Download URL
#' @param save_path Path to save the downloaded file
#' @return NULL
#' @importFrom httr GET write_disk
download_output <- function(url, save_path) {
  GET(url, write_disk(save_path, overwrite = TRUE))
}

#' Run the PanAzimuth API analysis
#'
#' @param input_data Input data for analysis
#' @param api_endpoint API endpoint URL
#' @return Analysis results
#' @export
run_azimuth_api <- function(input_data, api_endpoint) {
  # Implementation of the API analysis
  # This is a placeholder - implement the actual API interaction logic
  stop("Not implemented yet")
}

#' Process results from PanAzimuth API
#'
#' @param results Raw results from the API
#' @return Processed results in Seurat object format
#' @export
process_azimuth_results <- function(results) {
  # Implementation of results processing
  # This is a placeholder - implement the actual results processing logic
  stop("Not implemented yet")
} 