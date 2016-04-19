\name{mcmc.fmove}
\alias{mcmc.fmove}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Fit continuous-time functional movement model to telemetry data.%%  ~~function to do ... ~~
}
\description{
  Fits a functional movement model to telemetry data following Buderman
  et al., 2015.    
}
\usage{
mcmc.fmove(xy,t,fdabasis,tpred=t,QQ="CAR2",a=1,b=1,r=1,q=1,
           n.mcmc=100,num.paths.save=10,sigma.fixed=NA)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{xy}{
    A two-column matrix with each row corresponding to the x,y locations
    of a telemetry location.
  }
  \item{t}{
    A numeric vector of length = nrow(xy), with the i-th entry corresponding
    to the time of the i-th telemetry location in xy.
  }
  \item{fdabasis}{
    A "basisfd" object, typically resulting from a call to
    "create.bspline.basis" in the fda package.  Other basis functions
    can be used.
  }
  \item{tpred}{
    Numeric vector of times to impute the quasi-continuous path.
  }
  \item{QQ}{
    The precision matrix of the fda basis coefficients.  This can either
    be a string, taking on values of "CAR1" or "CAR2", or can be a user
    specified matrix (or sparse matrix using the Matrix package) of
    dimension equal to the number of basis functions in fdabasis.
    Defaults to "CAR2".  "CAR1" will result in less-smooth paths.
  }
  \item{a}{
    The shape parameter of the inverse gamma prior on the observation
    variance.
  }
  \item{b}{
    The scale parameter of the inverse gamma prior on the observation
    variance.
  }
  \item{r}{
    The shape parameter of the inverse gamma prior on the partial sill
    parameter of the spline basis coefficients.
  }
  \item{q}{
    The scale parameter of the inverse gamma prior on the partial sill
    parameter of the spline basis coefficients.
  }
  \item{n.mcmc}{
    Number of mcmc iterations to run.
  }
  \item{num.paths.save}{
    Number of quasi-continuous path realizations to save.  Defaults to
 10.
 }
 \item{sigma.fixed}{
   Numeric value (or the default NA).  If NA, then the observation
   variance sigma^2 is estimated using MCMC.  If a numeric value, this
   is the fixed standard deviation of the observation error.
   }
}
\details{
  Fits the functional movement model of Buderman et al., 2015, and
  outputs quasi-continuous paths that stochastically interpolate between
  telemetry locations.  The model fit is as follows (written out for 1-D):

  y_t = observed location at time t
  z_t = Sum_k beta_k*phi_k(t) = true location at time t, expressed using
  a linear combination of spline basis functions phi_k(t).

  y_t ~ N( z_t , sigma^2 )
  beta ~ N( 0 , tau^2 * QQ^{-1} )
  sigma^2 ~ IG(a,b)
  tau^2 ~ IG(r,q)
  
%%  ~~ If necessary, more details than the description above ~~
}
\value{
  A list with the following entries.
  1. s2.save = numeric vector of the values of sigma^2 at each mcmc
  iteration
  2. tau2.save = numeric vector of the values of tau^2 at each mcmc
  iteration
  3. pathlist = a list of length num.paths.save, with each item itself
  being a list with two entries:
    xy = a matrix with rows corresponding to x,y locations of the
    quasi-continuous path imputation
    t = a vector with entries corresponding to the times at which the
    quasi-continuous path was imputed
}
\references{
  Buderman, F.E.; Hooten, M. B.; Ivan, J. S. and Shenk, T. M.  A
  functional model for characterizing long-distance movement behavior.
  Methods in Ecology and Evolution, 2016, 7, 264-273. 
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
