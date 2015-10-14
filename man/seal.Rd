\name{seal}
\alias{seal}
\docType{data}
\title{
Data for one foraging trip by a male northern fur seal (NFS).%%   ~~ data name/kind ... ~~
}
\description{
  seal$locs contains xyt locations for argos fixes on the seal's
  location in the "datetime", "latitude", and "longitude" columns.

  seal$cov.df contains a data.frame of spatial covariate values for sea surface
  temperature (sst), chlorophyll A levels (chA) and net primary
  production (npp).

  seal$sim.obj is a fitted "sim.obj" obtained using the "crawl" package
  in R.
%%  ~~ A concise (1-5 lines) description of the dataset. ~~
}
\usage{data("seal")}
\format{
  The format is:
List of 10
 $ tripID    : num 156
 $ locs      :'data.frame':	163 obs. of  6 variables:
  ..$ datetime  : num [1:163] 36741 36741 36741 36742 36742 ...
  ..$ latitude  : num [1:163] 57.2 57.3 57.3 57.2 57.5 ...
  ..$ longitude : num [1:163] 190 190 190 190 190 ...
  ..$ landseamig: int [1:163] 0 1 1 1 1 1 1 1 1 1 ...
  ..$ lqadjust  : int [1:163] 5 1 0 -2 -2 1 -2 -2 -2 -2 ...
  ..$ lq        : Factor w/ 8 levels "0","1","2","3",..: 5 2 1 7 7 2 7 7 7 7 ...
 $ year      : int 2000
 $ rookery   : Factor w/ 2 levels "Reef","Vostochni": 2
 $ sex       : Factor w/ 2 levels "female","male": 2
 $ start.date: Factor w/ 131 levels "10/10/1999","10/11/1999",..: 80
 $ end.date  : Factor w/ 131 levels "10/10/1999","10/11/1999",..: 66
 $ cov.df    :List of 4
  ..$ X       :'data.frame':	10000 obs. of  5 variables:
  .. ..$ x  : num [1:10000] 184 184 184 184 184 ...
  .. ..$ y  : num [1:10000] 56.7 56.7 56.7 56.7 56.7 ...
  .. ..$ chA: num [1:10000] 1.19 0.924 0.744 0.709 0.733 ...
  .. ..$ sst: num [1:10000] 9.07 10.35 10.27 10.43 9.98 ...
  .. ..$ pro: num [1:10000] 853 821 823 849 886 ...
  ..$ X.grad.x: num [1:10000, 1:5] 0.0729 0.0729 0.0729 0.0729 0.0729 ...
  .. ..- attr(*, "dimnames")=List of 2
  .. .. ..$ : NULL
  .. .. ..$ : chr [1:5] "x" "y" "chA" "sst" ...
  ..$ X.grad.y: num [1:10000, 1:5] 0 0 0 0 0 0 0 0 0 0 ...
  .. ..- attr(*, "dimnames")=List of 2
  .. .. ..$ : NULL
  .. .. ..$ : chr [1:5] "x" "y" "chA" "sst" ...
  ..$ xy      : num [1:10000, 1:2] 184 184 184 184 184 ...
 $ sim.obj   :List of 31
  ..$ x            : num [1:5420] 190 9999 9999 9999 9999 ...
  ..$ y            : num [1:5420] 57.2 9999 9999 9999 9999 ...
  ..$ loctype      : num [1:5420] 0 1 1 1 1 1 1 1 1 1 ...
  ..$ P1.y         : num [1:2, 1:2] 1 0 0 1
  ..$ P1.x         : num [1:2, 1:2] 1 0 0 1
  ..$ a1.x         : num [1:2] 190 0
  ..$ a1.y         : num [1:2] 57.2 0
  ..$ n.errX       : int 2
  ..$ n.errY       : int 2
  ..$ n.mov        : int 1
  ..$ delta        : num [1:5420] 0.00417 0.00417 0.00417 0.00417 0.00417 ...
  ..$ lonAdj       : num [1:5420] 0.545 0.545 0.545 0.545 0.545 ...
  ..$ driftMod     : logi FALSE
  ..$ stopMod      : logi FALSE
  ..$ stop.mf      : NULL
  ..$ err.mfX      : num [1:5420, 1:2] 1 1 1 1 1 1 1 1 1 1 ...
  .. ..- attr(*, "dimnames")=List of 2
  .. .. ..$ : chr [1:5420] "1" "2" "3" "4" ...
  .. .. ..$ : chr [1:2] "(Intercept)" "errX"
  ..$ err.mfY      : num [1:5420, 1:2] 1 1 1 1 1 1 1 1 1 1 ...
  .. ..- attr(*, "dimnames")=List of 2
  .. .. ..$ : chr [1:5420] "1" "2" "3" "4" ...
  .. .. ..$ : chr [1:2] "(Intercept)" "errY"
  ..$ mov.mf       : num [1:5420, 1] 1 1 1 1 1 1 1 1 1 1 ...
  ..$ fixPar       : logi [1:6] NA NA NA NA NA NA
  ..$ Cmat         : num [1:6, 1:6] 0.116 -0.04442 0.0079 -0.00218 -0.00964 ...
  ..$ predType     : chr [1:5420] "p" "p" "p" "p" ...
  ..$ par          : num [1:6] -4.8167 0.5685 -6.9152 1.2077 0.0307 ...
  ..$ nms          : chr [1:6] "ln tau.x (Intercept)" "ln tau.x errX" "ln tau.y (Intercept)" "ln tau.y errY" ...
  ..$ N            : int 5420
  ..$ lower        : num -Inf
  ..$ upper        : num Inf
  ..$ loglik       : num 665
  ..$ datetime     : num [1:5420] 36741 36741 36741 36741 36741 ...
  ..$ coord        : chr [1:2] "longitude" "latitude"
  ..$ Time.name    : chr "datetime"
  ..$ thetaSampList:List of 1
  .. ..$ : num [1:100, 1:9] 0.00841 0.01441 0.02038 0.01096 0.00801 ...
  .. .. ..- attr(*, "dimnames")=List of 2
  .. .. .. ..$ : NULL
  .. .. .. ..$ : chr [1:9] "w" "lik" "prop.lik" "ln tau.x (Intercept)" ...
  .. .. ..- attr(*, "effSamp")= num 87
  .. .. ..- attr(*, "method")= chr "IS"
  ..- attr(*, "class")= chr "crwSimulator"
 $ pred.times: num [1:5258] 36741 36741 36741 36741 36741 ...
}
\details{
CRAWL sim.obj and xyt locations from ARGOS telemetry data for one NFS
  near the Pribilof islands.
  %%  ~~ If necessary, more details than the __description__ above ~~
}
\source{
Hanks, E.; Hooten, M.; Johnson, D. & Sterling, J. Velocity-Based
Movement Modeling for Individual and Population Level Inference PLoS
ONE, Public Library of Science, 2011, 6, e22795

  %%  ~~ reference to a publication or URL from which the data were obtained ~~
}
%% \references{
%% Hanks, E.; Hooten, M.; Johnson, D. & Sterling, J. Velocity-Based
%% Movement Modeling for Individual and Population Level Inference PLoS
%% ONE, Public Library of Science, 2011, 6, e22795
%% %%  ~~ possibly secondary sources and usages ~~
%% }
\examples{
## For example code, do
##
## > help(ctmcMove)
}
\keyword{datasets}
