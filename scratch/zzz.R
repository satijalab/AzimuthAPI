oldonLoad <- function(libname, pkgname) {
  # Check for required packages
  if (!requireNamespace("reticulate", quietly = TRUE) || 
      packageVersion("reticulate") < "1.40.0") {
    message(
      paste0("Installing or updating reticulate to version 1.40.0 alongwith ",
      "dependencies Rcpp and RcppTOML...")
      )
    install.packages("Rcpp")
    install.packages("RcppTOML")
    install.packages("reticulate")
    
    if (!requireNamespace("reticulate", quietly = TRUE) || 
        packageVersion("reticulate") < "1.40.0") {
      stop(
        paste0("Failed to install or update reticulate to version 1.40.0. ",
        "Please resolve this issue externally.")
        )
    } else {
      message("reticulate successfully updated to version 1.40.0.")
    }
  }
  
  if (!requireNamespace("Seurat", quietly = TRUE)) {
    stop("Seurat is not installed. Please install it to proceed.")
  } else {
    seurat_version <- as.numeric(substr(
      as.character(packageVersion("Seurat")),
      1,
      1
      ))
    
  }
} 
