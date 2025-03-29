#' Set up conda environment for PanAzimuth
#'
#' @param yml_file Path to conda environment YAML file
#' @param requirements_file Path to requirements.txt file
#' @param conda_path Path to conda executable (optional)
#' @param force_create Logical indicating whether to force environment creation
#' @param tensorflow_gpu Tensorflow version for gpu based usage
#' @param tensorflow_cpu Tensorflow version for cpu based usage
#' @return NULL
#' @export
setup_conda_env <- function(
  yml_file, 
  requirements_file, 
  conda_path = NULL, 
  force_create = FALSE,
  tensorflow_gpu = "tensorflow[and-cuda]==2.17",
  tensorflow_cpu = "tensorflow-cpu==2.17.0"
) {
  # Validate conda installation
  if (is.null(conda_path)) {
    conda_path <- reticulate::conda_binary()
  }
  if (is.null(conda_path)) {
    stop("Conda not found. Please install conda or miniconda and try again.")
  }
  
  # Parse and validate YAML environment file
  library(yaml)
  if (!file.exists(yml_file)) {
    stop(sprintf("YML file does not exist at the path: %s", yml_file))
  }
  yml_data <- tryCatch(
    yaml::read_yaml(yml_file),
    error = function(e) {
      stop(sprintf("Error reading YAML file: %s", e$message))
    }
  )
  env_name <- yml_data$name
  
  # Check if environment exists
  env_exists <- env_name %in% reticulate::conda_list()$name
  
  # Handle force_create flag
  if (force_create && env_exists) {
    message(sprintf(
      "Conda environment '%s' exists, but force_create is TRUE. Deleting it...",
      env_name
    ))
    system2(conda_path, c("env", "remove", "--name", env_name, "--yes"))
    message(sprintf("Environment '%s' deleted successfully.", env_name))
    env_exists <- FALSE
  }
  
  # Create environment if it doesn't exist
  if (!env_exists) {
    message(sprintf(
      "Creating conda environment '%s' from '%s'...", 
      env_name, 
      yml_file
    ))
    system2(conda_path, c("env", "create", "-f", yml_file))
    message(sprintf("Environment '%s' created.", env_name))
  }
  
  # Activate the environment
  message(sprintf("Activating environment '%s'...", env_name))
  reticulate::use_condaenv(
    condaenv = env_name,
    conda = conda_path,
    required = TRUE
  )
  message(sprintf("Environment '%s' is now active.", env_name))
  
  # Verify environment connection
  python_binary <- ifelse(
    .Platform$OS.type == "windows",
    file.path("Scripts", "python.exe"),
    "bin/python"
  )
  
  active_env <- reticulate::py_config()$python
  expected_env <- file.path(
    dirname(dirname(conda_path)), 
    "envs", 
    env_name, 
    python_binary
  )
  
  verify_environment <- function() {
    tryCatch({
      if (normalizePath(active_env) != normalizePath(expected_env)) {
        stop(sprintf(
          paste0(
            "Error: Reticulate is not connected to the expected conda ",
            "environment properly.\nExpected: %s\nActive: %s"
          ),
          expected_env, 
          active_env
        ))
      }
      message(sprintf("Reticulate is connected to: %s", active_env))
    }, error = function(e) {
      warning(sprintf(
        paste0(
          "Could not verify environment paths due to: %s\n",
          "Raw paths - Expected: %s, Active: %s"
        ), 
        e$message, expected_env, active_env
      ))
    })
  }
  
  verify_environment()
  
  # Install additional requirements if environment was just created
  if (!env_exists && 
      !is.null(requirements_file) && 
      file.exists(requirements_file)) {
    
    message(sprintf("Installing dependencies from %s...", requirements_file))
    
    # Install requirements file
    system2(conda_path, c(
      "run", "-n", env_name, "pip", "install", "-r", requirements_file
    ))
    
    # Upgrade pip
    system2(conda_path, c(
      "run", "-n", env_name, "pip", "install", "--upgrade", "pip"
    ))
    
    # Upgrade keras
    system2(conda_path, c(
      "run", "-n", env_name, "pip", "install", "--upgrade", "keras"
    ))
    
    # Install appropriate TensorFlow version
    if (if_gpu()) {
      message("GPU detected, installing GPU version of TensorFlow...")
      system2(conda_path, c(
        "run", "-n", env_name, "pip", "install", tensorflow_gpu
      ))
    } else {
      message("No GPU detected, installing CPU version of TensorFlow...")
      system2(conda_path, c(
        "run", "-n", env_name, "pip", "install", tensorflow_cpu
      ))
    }
    
    message("All dependencies installed through pip.")
  } else if (!env_exists && 
             (is.null(requirements_file) || !file.exists(requirements_file))) {
    message("Requirements file not provided or doesn't exist.")
  }
  
  # Print final Python configuration
  print(reticulate::py_config())
}