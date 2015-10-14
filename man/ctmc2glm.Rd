\name{ctmc2glm}
\alias{ctmc2glm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Convert a "ctmc" object into Poisson glm format.%%  ~~function to do ... ~~
}
\description{
  Transforms a "ctmc" object and covariate rasters into data suitable for analysis using
  Poisson GLM software (like glm in R). 
}
\usage{
ctmc2glm(ctmc, stack.static, stack.grad, crw = TRUE,
    normalize.gradients = FALSE, grad.point.decreasing = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ctmc}{
A "ctmc" object (output from "path2ctmc").
  }
  \item{stack.static}{
A rasterStack object, where each layer in the stack is a location-based covariate.
  }
  \item{stack.grad}{
    A rasterStack object, where each layer in the stack is a directional
    gradient-based covariate
}
\item{crw}{
  Logical.  If TRUE (default), an autocovariate is created for each cell
  visited in the CTMC movement path.  The autocovariate is a
  unit-length directional vector pointing from the center of the most
  recent grid cell to the center of the current grid cell. 
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
}
\details{
  This code creates one data row for each possible transition from each
  grid cell visited by the CTMC path.
%%  ~~ If necessary, more details than the description above ~~
}
\value{
  \item{z}{Response variable (either zero or 1) for analysis using GLM
    software.}
  \item{X}{Matrix of predictor variables for analysis using GLM
    software.  Created from the location-based
    and gradient-based covariates.}
  \item{tau}{Offset for each row in a Poisson GLM with log link.}
  \item{t}{Vector of the time each raster grid cell was entered} 
  
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
