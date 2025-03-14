#' Check for NVIDIA GPU availability
#'
#' @return Logical indicating whether an NVIDIA GPU is available
#' @export
if_gpu <- function() {
  gpu_check <- system("which nvidia-smi", intern = TRUE, ignore.stderr = TRUE)
  
  if (length(gpu_check) == 0) {
    return(FALSE)
  } else{
    return(TRUE)
  }
} 