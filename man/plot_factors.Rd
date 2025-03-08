% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plots.R
\name{plot_factors}
\alias{plot_factors}
\title{Plot GAUDI Factors}
\usage{
plot_factors(
  object,
  label_size = 0,
  draw_lines = FALSE,
  ad_hoc_label = NULL,
  palette = "magma"
)
}
\arguments{
\item{object}{A \code{GAUDIObject}.}

\item{label_size}{Numeric value specifying the size of labels on the plot.
A value of 0 means labels are not drawn.}

\item{draw_lines}{Logical value indicating whether to draw dashed lines at the median of UMAP coordinates.}

\item{ad_hoc_label}{Optional vector of labels to be used instead of the default cluster labels.}

\item{palette}{Character string specifying the color palette for the plot.
Supported palettes include 'inferno', 'plasma', and others from the viridis package.}
}
\value{
A ggplot object representing the factor plot.
}
\description{
This function creates a plot of factors derived from GAUDI.
It provides options to adjust label sizes, draw lines, and use ad-hoc labels for clustering.
}
