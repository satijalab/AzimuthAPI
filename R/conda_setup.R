#' Set up conda environment for PanAzimuth
#'
#' @param yml_file Path to conda environment YAML file
#' @param requirements_file Path to requirements.txt file
#' @param conda_path Path to conda executable (optional)
#' @param force_create Logical indicating whether to force environment creation
#' @return NULL
#' @export
setup_conda_env <- function(yml_file, requirements_file, conda_path=NULL, force_create=FALSE) {
  
  if (is.null(conda_path)) {
    conda_path <- reticulate::conda_binary()
  }
  if (is.null(conda_path)) {
    stop("Conda not found. Please install conda or miniconda and try again.")
  }
  
  env_name <- "AzimuthNN_min" 
  
  env_check_command <- sprintf('%s env list | grep "%s"', conda_path, env_name)
  env_check <- system(env_check_command, intern = TRUE)
  
  if (length(env_check) == 0 || force_create) {
    message(sprintf("Creating conda environment '%s'...", env_name))
    reticulate::conda_create(env_name, yml_file)
  } else {
    message(sprintf("Conda environment '%s' already exists.", env_name))
  }
  
  reticulate::use_condaenv(env_name)
  
  # Install additional requirements if specified
  if (!is.null(requirements_file) && file.exists(requirements_file)) {
    reticulate::py_install(requirements_file, envname = env_name)
  }
} 