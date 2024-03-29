% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sleuth_pre_processing.R
\name{sleuth_test_lrt}
\alias{sleuth_test_lrt}
\title{Automated function to run all possible Likewise Ratio tests on a given Sleuth object with models derived from the original formula used to run sleuth_fit()}
\usage{
sleuth_test_lrt(sleuth_obj)
}
\arguments{
\item{sleuth_obj}{An existing Sleuth object as generated by sleuth_prep() and fit by sleuth_fit()}
}
\value{
Results of the Likewise Ratio test assigned to the original Sleuth object input
}
\description{
Automated function to run all possible Likewise Ratio tests on a given Sleuth object with models derived from the original formula used to run sleuth_fit()
}
\examples{
# Given a Sleuth object, run a Likewise Ratio test on all possible models derived from the original formula used to run sleuth_prep()
sleuth_test_lrt(so)
}
