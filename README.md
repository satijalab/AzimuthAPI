# PanAzimuth

An R package providing an interface to the Pan-Azimuth Web API for single-cell RNA sequencing analysis.

## Installation

You can install the package using devtools:

```r
if (!requireNamespace("devtools", quietly = TRUE)) {
  install.packages("devtools")
}
devtools::install_github("yourusername/PanAzimuth")
```

## Dependencies

The package requires:
- R >= 4.0.0
- reticulate >= 1.40.0
- Seurat >= 5.0.0
- conda or miniconda installed on your system

## Usage

```r
library(PanAzimuth)

# Check for GPU availability
if_gpu()

# Set up conda environment
setup_conda_env(
  yml_file = "path/to/environment.yml",
  requirements_file = "path/to/requirements.txt"
)

# Run analysis
results <- run_azimuth_api(
  input_data = your_data,
  api_endpoint = "your_api_endpoint"
)

# Process results
processed_results <- process_azimuth_results(results)
```

## License

This package is licensed under the MIT License - see the LICENSE file for details. 