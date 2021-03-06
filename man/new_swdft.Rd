% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/swdft.R
\name{new_swdft}
\alias{new_swdft}
\title{Constructor function for class 'swdft'}
\usage{
new_swdft(a, x, n, type, pad, taper_type, taper, p, smooth, m, num_convs)
}
\arguments{
\item{a}{2D complex array of SWDFT coefficients. If there is smoothing, then
this represents the smoothed squared modulus coefficients.}

\item{x}{numeric input signal}

\item{n}{window size}

\item{type}{'fftw' or 'fft'}

\item{pad}{whether or not it was padded}

\item{taper_type}{type of taper}

\item{taper}{numeric values of the taper}

\item{p}{of cosine taper (if used)}

\item{smooth}{type of smoother}

\item{m}{width of kernel for smoothing (optional)}

\item{num_convs}{number of kernel convolutions (optional)}
}
\value{
list w/ the same elements as the arguments, an S3 object of class 'swdft'
}
\description{
Constructor function for class 'swdft'
}
