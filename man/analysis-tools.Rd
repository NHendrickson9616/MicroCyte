\name{analysis-tools}
\alias{cluster_count}
\alias{color_density}
\alias{compensate}
\alias{cull}
\alias{customize}
\alias{ergodic_slicing}
\alias{explore}
\alias{foci_caller}
\alias{foci_compact}
\alias{full_image_correlate}
\alias{get_density}
\alias{grapho}
\alias{heat_map_image}
\alias{heatmap_it_all}
\alias{icellate}
\alias{joiner}
\alias{line_analysis}
\alias{modthequad}
\alias{nearest_n}
\alias{normalizer}
\alias{plaque_data_collect}
\alias{plumeria}
\alias{pull_n_zip}
\alias{re_map}
\alias{run_ca}
\alias{save_png}
\alias{sirmixaplot}
\alias{smoke_stack}
\alias{ternaplot}
\title{MicroCyte Analysis Tools}
\description{
Functions for image-derived microcytometry analysis. Each function documents
its contract through explicit arguments and returns its result rather than
requiring caller-global objects.
}
\details{
Image-oriented functions require the optional \pkg{imager} and \pkg{magick}
packages. Plotting functions return ggplot-compatible plot objects or
structured lists containing plots and transformed data.
}
\keyword{utilities}
