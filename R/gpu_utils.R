#' Check for NVIDIA GPU availability
#'
#' @return Logical indicating whether an NVIDIA GPU is available
#' @export
if_gpu <- function() {
  if (.Platform$OS.type == "windows") {
    # Windows-specific check
    gpu_check <- tryCatch(
      system("where nvidia-smi", intern = TRUE, ignore.stderr = TRUE),
      error = function(e) character(0)
    )
  } else {
    # Unix-like systems (Linux, macOS)
    gpu_check <- tryCatch(
      system("which nvidia-smi", intern = TRUE, ignore.stderr = TRUE),
      error = function(e) character(0)
    )
  }
  
  return(length(gpu_check) > 0)
}