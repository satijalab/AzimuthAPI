# AzimuthAPI

[![Better Stack Badge](https://uptime.betterstack.com/status-badges/v1/monitor/2k7ub.svg)](https://uptime.betterstack.com/?utm_source=status_badge)

## Version 1.0.0 

An R package providing an interface to the Pan-human Azimuth neural network, enabling users to run cell type annotation on single-cell and spatial transcriptomics data.

> [!TIP]
> For more details, including an introductory vignette & function reference, visit https://satijalab.org/pan_human_azimuth.

Two options for annotation are available:

- `CloudAzimuth` - computation occurs on the cloud; requires no additional setup.
- `ANNotate` - computation occurs entirely locally in R via `reticulate`; requires a conda environment with [`panhumanpy`](https://pypi.org/project/panhumanpy/) installed.


## Installation

```r
if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes") 
}

# Install AzimuthAPI from GitHub
remotes::install_github("satijalab/AzimuthAPI")
```
