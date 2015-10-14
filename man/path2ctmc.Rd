\name{path2ctmc}
\alias{path2ctmc}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Function to turn a discrete-time continuous-space path into a CTMC path.
  %%  ~~function to do ... ~~
}
\description{
  This function takes a movement path defined by xyt values (not
  necessarily equally spaced in time), and converts it into a CTMC path
  (a continuous-time discrete-space path on grid cells in a raster).
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
path2ctmc(xy, t, rast)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{xy}{
    A matrix of x,y locations at T time points.
    %%     ~~Describe \code{xy} here~~
}
\item{t}{
  A vector of T times associated with the T locations in "xy".
%%     ~~Describe \code{t} here~~
}
\item{rast}{
  A raster object or raster stack object that will define the
  discrete-space grid cells for the CTMC movement path.
%%     ~~Describe \code{rast} here~~
}
}
\details{
  This takes a xyt path and turns it into a list of the embedded chain
  and residence times of a continuous time Markov chain walk on the
  graph.  "Diagonal" transitions are fixed with linear interpolation.
  %%  ~~ If necessary, more details than the description above ~~
}
\value{
  \item{ec}{A vector of the sequential grid cells (the embedded chain)
    in the CTMC movement path}
  \item{rt}{A vector of residence times in each sequential grid cell in
    the CTMC movement path}
  \item{trans.times}{A vector of times in which the movement path enters
    the grid cell in "ec".}
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
