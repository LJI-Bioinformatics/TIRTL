% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/TIRTL-package.R
\docType{package}
\name{TIRTL-package}
\alias{TIRTL}
\alias{TIRTL-package}
\title{TIRTL: MAD-HYPE and TSHELL analysis of TIRTL-seq clones}
\description{
This package provides functionality for analyzying sets of TCR clone files output by MIXCR v4.6 as part of a TIRTL-seq experiment.
}
\section{Installation}{

You can install the development version from GitHub with:

\if{html}{\out{<div class="sourceCode r">}}\preformatted{devtools::install_github("LJI-Bioinformatics/TIRTL")
}\if{html}{\out{</div>}}

\strong{Requirements}
This package will launch several Python scripts that have additional requirements.  The
requirements vary with the chosen backend and can be installed in a virtual environment
(or the global environemnt).

\emph{Numpy backend requirements}
\verb{pip install numpy pandas}

\emph{CuPy (CUDA) backend requirements}
\verb{pip install cupy}

\emph{Mac Graphics (MLX) backend requirements}
\verb{pip install mlx}
}

\section{Usage}{

The Python requirements listed above must be available in you environment when running
the code in this package.  This can be accomplished by activating your virtual environment
prior to launching R.
}

\section{Main Functions}{

\itemize{
\item \code{run_single_point_analysis_sub_gpu()}: Run a full analysis on a single timepoint
}
}

\author{
\strong{Maintainer}: Mikhail Pogorelyy \email{Mikhail.Pogorelyy@stjude.org} (\href{https://orcid.org/0000-0003-0773-1204}{ORCID})

Authors:
\itemize{
  \item Jason Greenbaum \email{jgbaum@gmail.com} (\href{https://orcid.org/0000-0002-1381-0390}{ORCID})
}

}
\keyword{internal}
