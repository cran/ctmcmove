\name{get.crawl.path}
\alias{get.crawl.path}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Simulate a CRAWL path within a raster object.
}
\description{
This function simulates a quasi-continuous CRAWL path that does not
leave the borders of a raster object.  Output includes the crawl xyt
locations and the discretized path on the raster grid. 
}
\usage{
get.crawl.path(sim.obj, raster, mintime, maxtime, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{sim.obj}{
A CRAWL sim.obj.
  }
  \item{raster}{
A raster or raster stack object.  The simulated CRAWL path will be
required to stay within the confines of the raster object.
}
  \item{mintime}{
    Time to begin simulating the CRAWL path.
  }
  \item{maxtime}{
    Time to end simulating the CRAWL path.
    %%     ~~Describe \code{maxtime} here~~
}
  \item{\dots}{
Legacy.  Not used.
}
}
\details{
  This is a wrapper function for CRAWL simulation functions and depends
  on a valid crawl sim.obj.
  %%  ~~ If necessary, more details than the description above ~~
}
\value{
  \item{t}{A vector of times at which the CTCRW path is sampled}
  \item{xy}{A 2-column matrix of xy locations of the CTCRW path}
  \item{cells.path}{A vector of the raster cell of the CTCRW location at each time in "t"}%%  ~Describe the value returned
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
