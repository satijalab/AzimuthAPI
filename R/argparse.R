#' Parse command line arguments for PanAzimuth
#'
#' @return Parsed arguments
#' @export
arg_parse_in_R <- function() {
  cat("Capturing arguments passed to R script...")
  cat("\n")
  
  parser <- ArgumentParser(description = "Argument parser for the script")
  
  parser$add_argument(
    "filepath",
    help = "Enter abs file path to the query. Query should be in h5ad format.",
    type = "character"
  )
  
  parser$add_argument(
    "-fn", "--feature_names_col",
    default = NULL,
    help = "Enter the column name where the feature names are stored in query.var where query is the anndata object read from the h5ad.",
    type = "character"
  )
  
  parser$add_argument(
    "-sdd", "--source_data_dir",
    default = paste0(python_module_path, "data/kfold_data"),
    help = "Source data directory",
    type = "character"
  )
  
  parser$add_argument(
    "-ft", "--features_txt",
    default = "features.txt",
    help = "Features text file",
    type = "character"
  )
  
  parser$add_argument(
    "-m", "--mode",
    default = "cumulative",
    help = "Enter the label mode: cumulative or independent",
    type = "character"
  )
  
  parser$add_argument(
    "-md", "--model",
    default = "M0.2",
    help = "Enter the model name",
    type = "character"
  )
  
  parser$add_argument(
    "-l", "--loss",
    default = "level_wt_focal_loss",
    help = "Enter the loss function used for optimization",
    type = "character"
  )
  
  parser$add_argument(
    "--epochs",
    default = 55,
    help = "Enter the number of epochs the model has been trained for",
    type = "integer"
  )
  
  parser$add_argument(
    "-ts", "--train_seed",
    default = 100,
    help = "Enter the training seed used",
    type = "integer"
  )
  
  parser$add_argument(
    "-ds", "--data_seed",
    default = 414,
    help = "Enter the data prep seed that was used",
    type = "integer"
  )
  
  parser$add_argument(
    "-dso", "--data_source",
    default = "data/kfold_data/datasets/fold10_02_26_2025_17_53_139",
    help = "Enter the source dataset with no / at either end",
    type = "character"
  )
  
  parser$add_argument(
    "-dsp", "--data_split",
    nargs = 3,
    default = c(7, 1, 2),
    help = "Enter the train:valid:test split as ints separated by a space",
    type = "integer"
  )
  
  parser$add_argument(
    "-ms", "--mask_seed",
    default = NULL,
    help = "Enter the seed used in the masking processes",
    type = "integer"
  )
  
  return(parser$parse_args())
} 