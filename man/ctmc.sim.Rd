\name{ctmc.sim}
\alias{ctmc.sim}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Code to simulate a continuous-time Markov chain.%%  ~~function to do ... ~~
}
\description{
  Simulates a CTMC with given rate matrix (Q) for a time (T), or until
  it reaches a final absorbing state.  
}
\usage{
ctmc.sim(Q,start.state=1,T=1,final.state=NA)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Q}{
A square matrix.  Either a rate matrix or the infinitessimal generator of the CTMC.
  }
  \item{start.state}{
    An integer - the starting state for the simulation.
  }
  \item{T}{
A numeric value greater than zero.  The time window for simulating the
CTMC will be [0,T].
  }
  \item{final.state}{
    Either NA or an integer.  If an integer, the chain will be simulated
    until it enters the "final.state", at which time the simulation will
    be terminated.
  }
}
\details{
  This code uses the Gillespie algorithm to simulate a CTMC path in
  continuous time.
%%  ~~ If necessary, more details than the description above ~~
}
\value{

  \item{ec}{A vector of the sequential grid cells (the embedded chain)
    in the CTMC movement path}
    \item{rt}{A vector of residence times in each sequential grid cell
    in the CTMC movement path}
  }
\references{
None%% ~put references to the literature/web site here ~
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
