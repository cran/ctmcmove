\name{get.UD}
\alias{get.UD}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Find the stationary distribution of the CTMC.%%  ~~function to do ... ~~
}
\description{
  Finds the stationary distribution (proportional utilization
  distribution) implied by a CTMC movement model with a given rate
  matrix.  
}
\usage{
get.UD(R,zero.idx=integer())
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{R}{
    Rate matrix with R[i,j] equal to the CTMC rate of movement from
    raster cell i to neighboring raster cell j.  R[i,j]=0 implies that
    cells i and j are not first order neighbors.
  }
  \item{zero.idx}{
    Integer vector of the indices of raster cells that are not passable and should be
    excluded.  These are cells where movement should be impossible.  Default is zero.idx=integer().
    }
}
\details{
  This calculates the stationary distribution of the CTMC using the
  method on pg. 455 of Harrod and Plemmons (1984).  If R is a sparse
  Matrix object, then sparse matrix methods are used, making this
  calculation extremely efficient.
  %%  ~~ If necessary, more details than the description above ~~
}
\value{
  Vector of the stationary distribution at each raster grid cell 
  
  %%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
  Harrod, W. J. & Plemmons, R. J.  Comparison of some direct methods for
  computing stationary distributions of Markov chains.  SIAM Journal on
  Scientific and Statistical Computing, 1984, 5, 453-469
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
