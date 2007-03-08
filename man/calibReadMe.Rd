\name{calibReadMe}
\alias{calibReadMe}
\title{ View CALIB readme file}
\description{
 Finds the location of the CALIB readme file and optionally opens it.
}
\usage{
calibReadMe(view = TRUE)
}
\arguments{
  \item{view}{ logical, \code{TRUE} means open the readme file and
        \code{FALSE} means finds out the location of the file only.}
}
\details{
 The function \code{vignette("limma")} will find the short CALIB 
 vignette which describes how to obtain the CALIB readme file. 
 The readme file is not a true vignette because it is not automatically
 generated using 'Sweave' during the package build process. This
 means that it cannot be found using \code{vignette}, hence the 
 need for this special function.

 If the operating system is other than Windows, then the PDF
 viewer used is the one given by \code{Sys.getenv("R\_PDFVIEWER")}. 
 The PDF viewer can be changed using \code{Sys.putenv(R_\PDFVIEWER=)}.
}
\value{
Character string giving the file location.
}
\references{\code{\link[limma]{limmaUsersGuide}} in the limma package}
\author{ Hui Zhao }
\examples{
calibReadMe(view=FALSE)
}
\keyword{ documentation }