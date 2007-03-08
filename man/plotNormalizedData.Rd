\name{plotNormalizedData}
\alias{plotNormalizedData}
\title{plot estimated absolute levels of two conditions}

\description{
  plot estimated absolute levels of any two user-specified conditions. The values in the plot
  are in log scale.
}

\usage{
plotNormalizedData(data, condition = c(1, 2), xlab = NULL, ylab = NULL,
                   main = NULL,xlim = NULL, ylim = NULL, pch = 19,
                   cex = 0.2,col = "black", diag=TRUE, diagcol="blue",
                   diaglwd=1.5, ...)
}

\arguments{
  \item{data}{ matrix containing estimated absolute levels. columns are conditons and
        rows are genes.}

  \item{condition}{ integer vector giving the two conditions to be plotted.}

  \item{xlab}{ a title for the x axis.}

  \item{ylab}{ a title for the y axis.}

  \item{main}{ an overall title for the plot. }

  \item{xlim}{ the x limits (min,max) of the plot.}

  \item{ylim}{ the y limits of the plot.}

  \item{pch}{ an integer code for one of a set of plotting characters 
        or symbols for the spike data set. Default is 19.}

  \item{cex}{ a numerical value giving the amount by which points
        should be scaled relative to the default. Default is 0.2.}

  \item{col}{ the color of the points. Default is black.}

  \item{diag}{ a logical value. Add diagonal on the plot if it is \code{TRUE}. 
        Default is \code{TRUE}.}

  \item{diagcol}{ the color of the diagonal. Default is blue.}

  \item{diaglwd}{ the width of the diagonal. Default is 1.5.}

  \item{\dots}{ other graphical parameters can be used in function
        \code{\link[graphics]{plot}}. }
}

\details{
  The function polts estimated absolute expression levels of two conditions. It accepts
  expression levels from the argument 'data', which should have the same data format as 
  the output value of the function \code{\link[CALIB:normalizedata]{NormalizeData}}.

  The two conditions to be plotted should be specified by the argument \code{condition}. 
  The \code{condition} should be a numeric vector with length two and it should be subset
  of condition vector of the design matrix. see function \code{NormalizeData}.

  see other graphic functions for the other arguments. 
}

\value{
   A plot is created on the current graphics device.
}

\author{ Hui Zhao}

\examples{
# load data: normalized data
data(normdata)
	
# specify the two conditions to be plotted.
cond <- c(1,2)
	
# use the default values for other parameters.
plotNormalizedData(normdata,condition = cond)
}

\keyword{ hplot }
