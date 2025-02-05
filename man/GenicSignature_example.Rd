% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/GenicSignature_example.R
\name{GenicSignature_example}
\alias{GenicSignature_example}
\title{GenicSignature_example}
\usage{
GenicSignature_example(data = "platelet_signature")
}
\arguments{
\item{data}{Name of the available signature that the user wants to test.}
}
\description{
Get path to molecular signature files used as example to illustrate how scSIGTest works.

scSIGTest comes bundled with some example files in its `inst/extdata`
directory. This function allows to get the path of files containing molecular signatures that can be used as input of the scSIGTest method.
}
\examples{

GenicSignature_example("platelet_signature")

}
