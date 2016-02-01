\name{Pctmc}
\alias{Pctmc}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Transition Matrix of a CTMC.%%  ~~function to do ... ~~
}
\description{
  Computes the transition matrix P(t) of a CTMC with given rate matrix (Q) and time (t).  
}
\usage{
Pctmc(Q,t)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Q}{
A square matrix.  Either a rate matrix or the infinitessimal generator of the CTMC.
  }
  \item{t}{
A numeric value - the time step.
  }
}
\details{
  This code creates one data row for each possible transition from each
  grid cell visited by the CTMC path.
%%  ~~ If necessary, more details than the description above ~~
}
\value{
A square matrix P with entries P[i,j]=Prob(X(t)=j|X(0)=i)
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
