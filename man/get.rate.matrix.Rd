\name{get.rate.matrix}
\alias{get.rate.matrix}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Create a CTMC rate matrix from rasters and parameter estimates.%%  ~~function to do ... ~~
}
\description{
  Creates a CTMC rate matrix from rasters and parameter estimates
  (perhaps from a GLM analysis). 
}
\usage{
get.rate.matrix(object, stack.static, stack.grad,
    normalize.gradients = FALSE, grad.point.decreasing = TRUE,
    directions=4, zero.idx=integer(), coef)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
A fitted GLM or GAM object used to fit the CTMC movement model
  }
  \item{stack.static}{
A rasterStack object, where each layer in the stack is a location-based covariate.
  }
  \item{stack.grad}{
    A rasterStack object, where each layer in the stack is a directional
    gradient-based covariate
}
  \item{normalize.gradients}{
    Logical.  Default is FALSE.  If TRUE, then all gradient
         covariates are normalized by dividing by the length of the gradient
	vector at each point.}
  \item{grad.point.decreasing}{
Logical.  If TRUE, then the gradient covariates are positive in the
direction of decreasing values of the covariate.  If FALSE, then the
gradient covariates are positive in the direction of increasing values
of the covariate (like a true gradient).
  }
  \item{directions}{
    Integer.  Either 4 (indicating a "Rook's neighborhood" of 4
    neighboring grid cells) or 8 (indicating a
    "King's neighborhood" of 8 neighboring grid cells).
  }
  \item{zero.idx}{
    Integer vector of the indices of raster cells that are not passable and should be
    excluded.  These are cells where movement should be impossible.  Default is zero.idx=integer().
    }
    \item{coef}{A vector of coefficents to use in place of those in 'object'}
}
\details{
  This function takes the covariate rasters in stack.static (motility
  covariates) and stack.grad (gradient covariates) and creates a CTMC
  rate matrix defining movement between all neighboring raster grid
  cells.  It is NOT possible to include an autocovariate here ("crw" in
  ctmc2glm).  If such was included in the original fitted model, then the crw term is set equal to zero.
%%  ~~ If necessary, more details than the description above ~~
}
\value{
  An n-by-n Matrix of CTMC rate values.
  %%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
Hanks, E. M.; Hooten, M. B. & Alldredge, M. W. Continuous-time
Discrete-space Models for Animal Movement.  The Annals of Applied
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
