\name{rast.grad}
\alias{rast.grad}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Creates gradient rasters from a raster object.%%  ~~function to do ... ~~
}
\description{
  This function takes a raster stack or raster object and creates two
  matrices for each raster layer, one which contains the x coordinates
  of the gradient of the raster layer and one which contains the y
  coordinates of the gradient of the raster layer.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
rast.grad(rasterstack)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{rasterstack}{A raster layer or raster stack from package "raster".
%%     ~~Describe \code{rasterstack} here~~
}
}
\details{
  The gradient is computed using the "terrain" function in raster.
  %%  ~~ If necessary, more details than the description above ~~
}
\value{
  \item{xy}{A matrix of x and y coordinates of each cell in the raster
    stack or raster layer.  The order is the order of the cells in the
    raster object.}
  \item{grad.x}{a matrix where each column is the x-coordinates of the
    gradient for one raster layer}
  \item{grad.y}{a matrix where each column is the y-coordinates of the
    gradient for one raster layer}
  \item{rast.grad.x}{A raster stack where each raster layer is the
    x-coordinates of the gradient for one covariate}
  \item{rast.grad.y}{A raster stack where each raster layer is the
    x-coordinates of the gradient for one covariate}
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
Hanks, E. M.; Hooten, M. B. & Alldredge, M. W. Continuous-time
Discrete-space Models for Animal Movement The Annals of Applied
Statistics, 2015, 9, 145-165
%% ~put references to the literature/web site here ~
}
\author{
  Ephraim M. Hanks
%%  ~~who you are~~
}
%% \note{
%% %%  ~~further notes~~
%% }

%% ~Make other sections like Warning with \section{Warning }{....} ~

%% \seealso{
%% %% ~~objects to See Also as \code{\link{help}}, ~~~
%% }
\examples{
## For example code, do
##
## > help(ctmcMove)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
