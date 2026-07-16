#' Helper function for `listen_to_progress`
#' Safely process SSE stream with error handling
#'
#' @param url API endpoint URL
#' @param handle Curl handle with the appropriate form data set
#' @return Logical indicating success (TRUE) or failure (FALSE)
#' @importFrom curl curl_fetch_stream
#' @importFrom jsonlite fromJSON
#' @noRd
safe_progress_stream <- function(url, handle) {
  tryCatch({
    result <- list(success = NULL, output_file = NULL, download_url = NULL, job_id = NULL)
    progress_id <- NULL

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
              
              # Check for success flag
              if (!is.null(parsed$success)) {
                result$success <<- parsed$success
              }
              if (!is.null(parsed$output_file)) {
                result$output_file <<- parsed$output_file
              }
              if (!is.null(parsed$download_url)) {
                result$download_url <<- parsed$download_url
              }
              if (!is.null(parsed$job_id)) {
                result$job_id <<- parsed$job_id
              }
              
              # Print messages based on content
              if (!is.null(parsed$message)) {
                if (isTRUE(parsed$queued)) {
                  cli::cli_alert_warning("[Queued] {parsed$message}")
                } else if (grepl("successfully", parsed$message, ignore.case = TRUE)) {
                  cli::cli_alert_success("{parsed$message}")
                } else {
                  cli::cli_alert_info("{parsed$message}")
                }
              }

              h2_strings <- c("panhumanpy version", "reference model and parameters", "running model")
              h2_strings <- paste(h2_strings, collapse = "|")

              if (!is.null(parsed$output)) {
                if (grepl(h2_strings, parsed$output, ignore.case = TRUE)) {
                  cli::cli_h2("{parsed$output}")
                } else if (grepl("successfully", parsed$output, ignore.case = TRUE)) {
                  cli::cli_alert_success("{parsed$output}")
                } else {
                  cli::cli_text("{parsed$output}")
                }
              }
              if (!is.null(parsed$error)) {
                cli::cli_alert_danger("{parsed$error}")
              }
              if (!is.null(parsed$progress)) {
                if (is.null(progress_id)) {
                  progress_id <<- cli::cli_progress_bar(
                    name = "Processing",
                    total = 100,
                    clear = FALSE,
                    format = "{name} [{bar}] {percent}%"
                  )
                }
                cli::cli_progress_update(id = progress_id, set = as.numeric(parsed$progress))
              }
              if (!is.null(parsed$status)) {
                cli::cli_alert_info("Status: {parsed$status}")
              }
            }, error = function(e) {
              # Only print actual parsing errors, not empty or malformed messages
              if (!grepl("unexpected end", e$message)) {
                cli::cli_warn("Failed to parse message: {e$message}")
                cli::cli_text("Raw message: {json_message}")
              }
            })
          }
        }
    })
    if (!is.null(progress_id)) {
      cli::cli_progress_done(id = progress_id)
    }
    if (is.null(result$success)) {
      cli::cli_alert_danger("No completion signal received from the server.")
      result$success <- FALSE
      return(result)
    }
    return(result)
  }, error = function(e) {
    cli::cli_alert_danger("Error during progress streaming: {conditionMessage(e)}")
    return(list(success = FALSE, output_file = NULL, download_url = NULL, job_id = NULL))
  })
}


#' Listen to progress updates from the API
#'
#' @param url API endpoint URL
#' @param file_path Path to the file being processed
#' @param ... Additional arguments passed to the API
#' @return Logical indicating success (TRUE) or failure (FALSE)
#' @importFrom curl new_handle handle_setform handle_setopt curl_fetch_stream form_file
#' @importFrom jsonlite fromJSON
#' @keywords internal
#' @noRd
listen_to_progress <- function(url, file_path, ...) {
  # Create a multipart form for the upload
  additional_args <- list(...)
  additional_args <- lapply(additional_args, as.character)

  form <- c(
    list(file = form_file(file_path)),
    additional_args
  )

  # Open a connection to the API using curl
  handle <- new_handle()
  handle_setopt(handle, tcp_keepalive = 1L)
  handle_setform(handle, .list = form)

  # Track whether processing succeeded
  result <- safe_progress_stream(url, handle)

  return(result)
}

#' Download output file from the API
#'
#' @param url Download URL
#' @param save_path Path to save the downloaded file
#' @return NULL
#' @importFrom httr GET write_disk status_code
#' @noRd
download_output <- function(url, save_path) {
  response <- GET(url, write_disk(save_path, overwrite = TRUE))
  if (status_code(response) >= 300) {
    stop("Failed to download output file from the server.")
  }
}

#' Check local vs latest API version / ensure connection can be established to the server
#'
#' @param api_base_url Base URL for the API
#'
#' @importFrom httr GET content status_code
#' @importFrom utils packageVersion
#'
#' @return NULL
#' @keywords internal
#' @noRd
check_api_version <- function(api_base_url) {
  version_url <- paste0(api_base_url, "/version")
  response <- tryCatch({
    GET(version_url, timeout(5))
  }, error = function(e) {
    if (grepl("Could not connect to server", e$message, fixed = TRUE)) {
      stop(simpleError(
          "Connection refused: server not running or port closed.\nPlease report at https://github.com/satijalab/AzimuthAPI/issues.",
          call = conditionCall(e)
      ))
    }
    stop(simpleError(paste0("Error connecting to the server: ", conditionMessage(e)), call = conditionCall(e)))
  })
  
  if (status_code(response) != 200) {
    stop(simpleError(paste0("Failed to retrieve version information from the server (HTTP ", status_code(response), ")."), call = NULL))
  }

  latest_version <- content(response)$version
  current_version <- as.character(packageVersion("AzimuthAPI"))
  return(utils::compareVersion(current_version, latest_version) < 0)
}