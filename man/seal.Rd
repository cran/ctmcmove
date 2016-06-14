\name{seal}
\alias{seal}
\docType{data}
\title{
Data for one foraging trip by a male northern fur seal (NFS).%%   ~~ data name/kind ... ~~
}
\description{
  seal$locs contains xyt locations for ARGOS fixes on the seal's
  location in the "datetime", "latitude", and "longitude" columns.

  seal$cov.df contains a data.frame of spatial covariate values for sea surface
  temperature (sst), chlorophyll A levels (chA) and net primary
  production (npp).

%%  ~~ A concise (1-5 lines) description of the dataset. ~~
}
\usage{data("seal")}
\format{
  The format is:
  
 $ locs      :'data.frame':	163 obs. of  6 variables:
 
  ..$ datetime  : num [1:163] 36741 36741 36741 36742 36742 ...
  
  ..$ latitude  : num [1:163] 57.2 57.3 57.3 57.2 57.5 ...
  
  ..$ longitude : num [1:163] 190 190 190 190 190 ...
  
  ..$ landseamig: int [1:163] 0 1 1 1 1 1 1 1 1 1 ...
  
  ..$ lqadjust  : int [1:163] 5 1 0 -2 -2 1 -2 -2 -2 -2 ...
  
  ..$ lq        : Factor w/ 8 levels "0","1","2","3",..: 5 2 1 7 7 2 7 7 7 7 ...
 
 $ sex       : Factor w/ 2 levels "female","male": 2
 
  $ cov.df    :List of 4
 
  ..$ X       :'data.frame':	10000 obs. of  5 variables:
 
  .. ..$ x  : num [1:10000] 184 184 184 184 184 ...
 
  .. ..$ y  : num [1:10000] 56.7 56.7 56.7 56.7 56.7 ...
 
  .. ..$ chA: num [1:10000] 1.19 0.924 0.744 0.709 0.733 ...
 
  .. ..$ sst: num [1:10000] 9.07 10.35 10.27 10.43 9.98 ...
 
  .. ..$ pro: num [1:10000] 853 821 823 849 886 ...
}
\details{
Covariate Rasters and ARGOS telemetry data for one NFS
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
## > help(ctmcmove)
}
\keyword{datasets}
