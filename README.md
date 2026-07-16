# AzimuthAPI

## Version 0.9.0

An R package providing an interface to the Pan-human Azimuth neural network, enabling users to run cell type annotation on single-cell and spatial transcriptomics data.

> [!TIP]
> For more details, including an introductory vignette & function reference, visit https://satijalab.org/pan_human_azimuth.

Two options for annotation are available:

- `CloudAzimuth`: computation occurs on the cloud; requires no additional setup.
- `ANNotate`: computation occurs entirely locally; requires a working Python installation with [`panhumanpy`](https://pypi.org/project/panhumanpy/) installed.

## Installation

```r
# Install devtools if not already installed
install.packages("remotes")

# Install AzimuthAPI from GitHub
remotes::install_github("satijalab/AzimuthAPI")
```

