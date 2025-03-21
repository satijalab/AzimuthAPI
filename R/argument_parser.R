#' Parse command line arguments for the ANNotate function
#'
#' @return A list of parsed arguments
#'
#' @examples
#' \dontrun{
#' args <- parse_annotate_args()
#' }
parse_annotate_args <- function() {
  # Check if argparse is available
  if (!requireNamespace("argparse", quietly = TRUE)) {
    install.packages("argparse")
    library(argparse)
  } else {
    library(argparse)
  }
  
  message("Parsing command line arguments...")
  
  # Create argument parser
  parser <- ArgumentParser(description = "Command line interface for ANNotate function")
  
  # Required arguments
  parser$add_argument(
    "filepath",
    help = "Path to the Seurat object RDS file",
    type = "character"
  )
  
  # Optional arguments 
  parser$add_argument(
    "-f", "--feature_names_col",
    default = NULL,
    help = "Column name containing feature names (default: NULL)",
    type = "character"
  )

  parser$add_argument(
    "-p", "--annotation_pipeline",
    default = "supervised",
    help = "Annotation pipeline to use (default: 'supervised')",
    type = "character"
  )  
  
  parser$add_argument(
    "-ebs", "--eval_batch_size",
    default = 40960,
    help = "Batch size for evaluation (default: 40960)",
    type = "integer"
  )
  
  parser$add_argument(
    "-no", "--normalization_override",
    action = "store_true",
    default = FALSE,
    help = "Override normalization (default: FALSE)"
  )
  
  parser$add_argument(
    "-ncb", "--norm_check_batch_size",
    default = 1000,
    help = "Batch size for normalization check (default: 1000)",
    type = "integer"
  )
  
  parser$add_argument(
    "-om", "--output_mode",
    default = "minimal",
    help = "Output mode (default: 'minimal')",
    type = "character"
  )
  
  parser$add_argument(
    "-rl", "--refine_labels",
    action = "store_true",
    default = TRUE,
    help = "Refine labels (default: TRUE)"
  )  
  
  parser$add_argument(
    "-ee", "--extract_embeddings",
    action = "store_true",
    default = TRUE,
    help = "Extract embeddings (default: TRUE)"
  )
  
  parser$add_argument(
    "-ue", "--umap_embeddings",
    action = "store_true",
    default = TRUE,
    help = "Generate UMAP embeddings (default: TRUE)"
  )
  
  parser$add_argument(
    "-nn", "--n_neighbors",
    default = 30,
    help = "Number of neighbors for UMAP (default: 30)",
    type = "integer"
  )
  
  parser$add_argument(
    "-nc", "--n_components",
    default = 2,
    help = "Number of components for UMAP (default: 2)",
    type = "integer"
  )
  
  parser$add_argument(
    "-m", "--metric",
    default = "cosine",
    help = "Distance metric for UMAP (default: 'cosine')",
    type = "character"
  )
  
  parser$add_argument(
    "-md", "--min_dist",
    default = 0.3,
    help = "Minimum distance for UMAP (default: 0.3)",
    type = "double"
  )
  
  parser$add_argument(
    "-ulr", "--umap_lr",
    default = 1.0,
    help = "Learning rate for UMAP (default: 1.0)",
    type = "double"
  )
  
  parser$add_argument(
    "-us", "--umap_seed",
    default = 42,
    help = "Random seed for UMAP (default: 42)",
    type = "integer"
  )
  
  parser$add_argument(
    "-sp", "--spread",
    default = 1.0,
    help = "Spread parameter for UMAP (default: 1.0)",
    type = "double"
  )
  
  parser$add_argument(
    "-v", "--verbose",
    action = "store_true",
    default = TRUE,
    help = "Verbose output (default: TRUE)"
  )
  
  parser$add_argument(
    "-i", "--init",
    default = "spectral",
    help = "Initialization method for UMAP (default: 'spectral')",
    type = "character"
  )
  
  parser$add_argument(
    "-po", "--process_obj",
    action = "store_true",
    default = TRUE,
    help = "Process object with PrepLabel (default: TRUE)"
  )
  
  parser$add_argument(
    "-ca", "--cutoff_abs",
    default = 5,
    help = "Absolute cutoff for PrepLabel (default: 5)",
    type = "integer"
  )
  
  parser$add_argument(
    "-cf", "--cutoff_frac",
    default = 0.001,
    help = "Fractional cutoff for PrepLabel (default: 0.001)",
    type = "double"
  )
  
  args <- parser$parse_args()
  
  return(args)
}



#' Format parsed arguments for ANNotate function
#'
#' @param args List of parsed arguments from parse_annotate_args()
#' @return Named list of formatted arguments for ANNotate function
#'
#' @examples
#' \dontrun{
#' args <- parse_annotate_args()
#' formatted_args <- format_annotate_args(args)
#' }
format_annotate_args <- function(args) {
  message("Formatting arguments for ANNotate function...")
  
  # Convert arguments to properly named list for ANNotate
  formatted_args <- list(
    # Required arguments
    filepath = args$filepath,
    
    # Optional arguments with parameter names matching ANNotate function
    feature_names_col = args$feature_names_col,
    annotation_pipeline = args$annotation_pipeline,
    eval_batch_size = args$eval_batch_size,
    normalization_override = args$normalization_override,
    norm_check_batch_size = args$norm_check_batch_size,
    output_mode = args$output_mode,
    refine_labels = args$refine_labels,
    extract_embeddings = args$extract_embeddings,
    umap_embeddings = args$umap_embeddings,
    n_neighbors = args$n_neighbors,
    n_components = args$n_components,
    metric = args$metric,
    min_dist = args$min_dist,
    umap_lr = args$umap_lr,
    umap_seed = args$umap_seed,
    spread = args$spread,
    verbose = args$verbose,
    init = args$init,
    process_obj = args$process_obj,
    cutoff_abs = args$cutoff_abs,
    cutoff_frac = args$cutoff_frac,
    
    )
  
  # Inform about parsed arguments
  if (args$verbose) {
    message("Arguments for ANNotate function:")
    for (name in names(formatted_args)) {
      value <- formatted_args[[name]]
      if (is.null(value)) {
        message(sprintf("  %s: NULL", name))
      } else if (is.logical(value) || is.numeric(value)) {
        message(sprintf("  %s: %s", name, value))
      } else {
        message(sprintf("  %s: '%s'", name, value))
      }
    }
  }
  
  return(formatted_args)
}