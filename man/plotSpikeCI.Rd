\name{plotSpikeCI}
\alias{plotSpikeCI}
\title{ plot spike concentration vs measured intensity}

\description{
  plot spike known concentration and measured intensity of one array. 
}

\usage{
plotSpikeCI(spike, parameter,array = 1, bc = FALSE, area = TRUE, 
            meanpoint = TRUE,xlab = "log(Concentration)",
            ylab = "log(Intensity)", main = colnames(spike$R)[array],
            onlycalib = TRUE, xlim = NULL, ylim = NULL, pch = 19,
            cex = 0.2, meanpch = 21, meancex = 1, lwd = 1.5, 
            cy5col = "red", cy3col = "green", ...)
}

\arguments{
  \item{spike}{ a \code{SpikeList} object.}

  \item{parameter}{a \code{ParameterList} object.}

  \item{array}{ integer giving the array to be plotted. }

  \item{bc}{a logical value. \code{TRUE} means background corrected measured 
        intensities are used. Default is FALSE.}

  \item{area}{  a logical value. \code{TRUE} means spot area is used to calculate 
        measured intensities. Namly, measured intensities are calculated by 
        foreground intensities(or background corrected intensities, if bc is 
        \code{TRUE}) multiply spot area. \code{FALSE} means spot area is not used.
        Default is \code{TRUE}.}

  \item{meanpoint}{ a logical value. \code{TRUE} is to show meanpoint of measured
        intensities with the same concentration on the plot. \code{FALSE} means not show. }

  \item{xlab}{ a title for the x axis.}

  \item{ylab}{ a title for the y axis.}

  \item{main}{ an overall title for the plot.}

  \item{onlycalib}{ a logical value. \code{TRUE} means only the calibration controls are on 
        the plot. \code{FALSE} means to plot all the spikes}

  \item{xlim}{ the x limits (min,max) of the plot.}

  \item{ylim}{ the y limits of the plot.}

  \item{pch}{ a integer code for one of plotting characters or symbols
        for the spike data set. Default is 21.}

  \item{cex}{ a numerical value giving the amount by which the points which
        indicate spike data set should be scaled relative to the default. Default is 0.4.}

  \item{meanpch}{ a integer code for one of plotting characters or symbols
        for the meanpoints. Default is 21.}

  \item{meancex}{ a numerical value giving the amount by which the meanpoints
        should be scaled relative to the default value. Default is 1.}

  \item{lwd}{ width of the model curves. Default is 1.5.}

  \item{cy5col}{ color of all symbols for cy5. Default is red.}

  \item{cy3col}{ color of all symbols for cy3. Default is green.}

  \item{\dots}{ other graphical parameters can be used in function \code{\link[graphics]{plot}}.}
}

\details{
  The function plots spike concentration and measured intensity of one array. 
  array number is specified by the argument \code{array}. It accepts the concentration
  of given array from the agrument \code{spike}, which is a \code{\link[CALIB:spike]{SpikeList}}
  object. The measured intensities are calculated from \code{spike}. Four different ways can
  be used to calculate the measured intensities. Arguments \code{bc} and \code{area} are 
  logical and their combinations are used for specifying the four differents ways. 
  \code{bc} indicates using background correction or not. \code{area} indicates multipling 
  spot area or not. The default value of these two arguments are \code{bc} = FALSE and 
  \code{area} = TRUE. 

  In order to help data visualization, meanpoints and model curve can be added
  on the plot. And the arguments \code{meanpoint} and \code{parameter} are correspond to these.
  The meadians of every group of measured intensities which have the same concentration 
  are shown on the polt if \code{meanpoint} is true. Model curves of both dye are 
  shown if the arguement \code{parameter} is specified after parameter estimation. 
}

\value{
  A plot is created on the current graphics device.
}

\author{ Hui Zhao}

\seealso{ see graphic functions \code{\link[graphics]{plot}}, \code{\link[graphics]{par}}}}

\examples{
# load data: spike
data(spike)

# specify the array to be plotted.
array <- 1

# use the default values for other parameters.
plotSpikeCI(spike,array=array)

# after parameter estimation, the model curves can be shown on the plot.
data(parameter)
plotSpikeCI(spike,parameter,array=array)
}

\keyword{ hplot}
