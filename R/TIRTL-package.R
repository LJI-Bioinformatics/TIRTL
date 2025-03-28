#' @section Installation:
#' You can install the development version from GitHub with:
#' ```r
#' devtools::install_github("LJI-Bioinformatics/TIRTL")
#' ```
#' **Requirements**
#' This package will launch several Python scripts that have additional requirements.  The
#' requirements vary with the chosen backend and can be installed in a virtual environment
#' (or the global environemnt).
#'
#' *Numpy backend requirements*
#' ```pip install numpy pandas```
#'
#' *CuPy (CUDA) backend requirements*
#' ```pip install cupy```
#'
#' *Mac Graphics (MLX) backend requirements*
#' ```pip install mlx```
#'
#' @section Usage:
#' The Python requirements listed above must be available in you environment when running
#' the code in this package.  This can be accomplished by activating your virtual environment
#' prior to launching R.
#'
#' @section Main Functions:
#' - `run_single_point_analysis_sub_gpu()`: Run a full analysis on a single timepoint
#'
#' @keywords internal
"_PACKAGE"
