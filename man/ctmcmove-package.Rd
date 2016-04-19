\name{ctmcmove-package}
\alias{ctmcmove-package}
\alias{ctmcmove}
\docType{package}
\title{
ctmcmove
}
\description{
  Software to facilitates taking movement data in xyt format and pairing it with raster covariates within a continuous time Markov chain (CTMC) framework.  As described in Hanks et al. (2015) <DOI:10.1214/14-AOAS803> , this allows flexible modeling of movement in response to covariates (or covariate gradients) with model fitting possible within a Poisson GLM framework.
}
\details{

Typical work flow for analysis of telemetry / GPS movement data:

  1. Fit a quasi-continuous path model to telemetry xyt data.  The
  ctmcmove package facilitates this through the "mcmc.fmove" function.
  
  2. Create or import raster layers (from package "raster") for each
  covariate.
  
  3. Impute a quasi-continuous path (done jointly with model fitting in
  the "mcmc.fmove" function.
  
  4. Turn the quasi-continuous path into a CTMC discrete-space
  path using the "path2ctmc" command.

  5. Turn discrete-space path into Poisson GLM format using the
  "ctmc2glm" command.

  6. Repeat #3 - #5 multiple times (M times).  Stack together the response "z",
  model matrix "X", and offset "tau" elements from each imputed path. 

  7. Fit a Poisson GLM model to the stacked data with response "z",
  model matrix "X", offset "log(tau)", and weights for each row equal to "1/M".

  7 (alternate).  Alternately, multiple imputation could be used, as
  described in Hanks et al., (2015).
}
\author{
  Ephraim M. Hanks

Maintainer: Ephraim M. Hanks
}
\references{
Hanks, E. M.; Hooten, M. B. & Alldredge, M. W. Continuous-time
Discrete-space Models for Animal Movement The Annals of Applied
Statistics, 2015, 9, 145-165

Hanks, E.; Hooten, M.; Johnson, D. & Sterling, J. Velocity-Based
Movement Modeling for Individual and Population Level Inference PLoS
ONE, Public Library of Science, 2011, 6, e22795

Hooten, M. B.; Johnson, D. S.; Hanks, E. M. & Lowry, J. H. Agent-Based
Inference for Animal Movement and Selection Journal of Agricultural,
Biological, and Environmental Statistics, 2010, 15, 523-538

}
\keyword{ Animal Movement }
%% \seealso{
%% }
\examples{
\dontrun{

##
## Example of using a CTMC model for movement
##
## Steps:
##  1. Fit Quasi-Continuous Path Model to telemetry data (done using Buderman et al 2015)
##  2. Create covariate raster objects (the CTMC will be on the raster
##     grid cells)
##  3. Impute a quasi-continuous path
##  4. Turn quasi-continuous path into a CTMC discrete-space path
##  5. Turn discrete-space path into latent Poisson GLM format
##  6. Fit a Poisson GLM model to the data
##

library(ctmcmove)
data(seal)
xyt=seal$locs[,3:1]
head(xyt)
plot(xyt[,1:2],type="b")
xy=xyt[,-3]
x=xyt[,1]
y=xyt[,2]
t=xyt[,3]


########################  


##########################################################################
##
## 1. Fit functional movement model to telemetry data
##
##########################################################################

library(fda)

## Define the knots of the spline expansion.
##
## Problems with fitting the functional movement model can often be fixed by
## varying the spacing of the knots.
knots = seq(min(t),max(t),by=1/6)
## create B-spline basis vectors used to approximate the path
b=create.bspline.basis(c(min(t),max(t)),breaks=knots,norder=3)
## define the sequence of times on which to sample the imputed path
tpred=seq(min(t),max(t),by=1/24/60)



## Fit latent Gaussian model using MCMC
out=mcmc.fmove(xy,t,b,tpred,QQ="CAR2",n.mcmc=400,a=1,r=1,num.paths.save=20)
str(out)

## plot 3 imputed paths
plot(xy,type="b")
points(out$pathlist[[1]]$xy,col="red",type="l")
points(out$pathlist[[2]]$xy,col="blue",type="l")
points(out$pathlist[[3]]$xy,col="green",type="l")


##########################################################################
##
## 2. Creating rasters
##
##########################################################################

cov.df=seal$cov.df
str(cov.df)

NN=sqrt(nrow(cov.df$X))
sst=matrix(seal$cov.df$X$sst,NN,byrow=TRUE)
sst=sst[NN:1,]
sst=raster(sst,xmn=min(seal$cov.df$X$x),xmx=max(seal$cov.df$X$x),
           ymn=min(seal$cov.df$X$y),ymx=max(seal$cov.df$X$y))


crs(sst)="+proj=longlat +datum=WGS84"
plot(sst)

chA=matrix(seal$cov.df$X$chA,NN,byrow=TRUE)
chA=chA[NN:1,]
chA=raster(chA,xmn=min(seal$cov.df$X$x),xmx=max(seal$cov.df$X$x),
           ymn=min(seal$cov.df$X$y),ymx=max(seal$cov.df$X$y))
crs(chA)="+proj=longlat +datum=WGS84"

pro=matrix(seal$cov.df$X$pro,NN,byrow=TRUE)
pro=pro[NN:1,]
npp=raster(pro,xmn=min(seal$cov.df$X$x),xmx=max(seal$cov.df$X$x),
           ymn=min(seal$cov.df$X$y),ymx=max(seal$cov.df$X$y))
crs(npp)="+proj=longlat +datum=WGS84"


int=sst
values(int) <- 1

d2r=int
rookery.cell=cellFromXY(int,xyt[1,1:2])
values(d2r)=NA
values(d2r)[rookery.cell]=0
d2r=distance(d2r)

grad.stack=stack(sst,chA,npp,d2r)
names(grad.stack) <- c("sst","cha","npp","d2r")

plot(sst)
points(xyt[,1:2],type="b")

plot(grad.stack)


##########################################################################
##
## 3 Impute Quasi-Continuous Paths
##
##########################################################################

plot(sst)
for(i in 1:10){
    points(out$pathlist[[i]]$xy,col=i,type="l")
}
points(xyt[,1:2],type="b",pch=20,cex=2,lwd=2)

##########################################################################
##
## 4. Turn continuous space path into a CTMC discrete space path
##
##########################################################################

path=out$pathlist[[1]]
ctmc=path2ctmc(path$xy,path$t,int)

str(ctmc)

##########################################################################
##
## 5. Turn CTMC discrete path into latent Poisson GLM data
##
##########################################################################

loc.stack=stack(int,sst)
names(loc.stack) <- c("Intercept","sst.loc")

glm.data=ctmc2glm(ctmc,loc.stack,grad.stack)

str(glm.data)
summary(glm.data)

## now repeat 3-5 for 10 imputations
for(i in 2:10){
    path=out$pathlist[[i]]
    ctmc=path2ctmc(path$xy,path$t,int)
    glm.data=rbind(glm.data,ctmc2glm(ctmc,loc.stack,grad.stack))
}

str(glm.data)

## remove transitions that are nearly instantaneous

summary(glm.data)
idx.0=which(glm.data$tau<10^-5)
idx.0
if(length(idx.0)>0){
    glm.data=glm.data[-idx.0,]
}

summary(glm.data)

##########################################################################
##
## 6. Fit Poisson GLM and Poisson GAM with time-varying coefficients
##    (here we are fitting all "M" paths simultaneously,
##     giving each one a weight of "1/M")
##
##########################################################################



fit=glm(z~cha+npp+sst+crw+d2r+sst.loc,
        weights=rep(1/10,nrow(glm.data)),family="poisson",offset=log(tau),data=glm.data)
summary(fit)

##
## 6. (Alternate) We can use any software which fits Poisson glm data.
##    The following uses "gam" in package "mgcv" to fit a time-varying
##    effect of "d2r" using penalized regression splines.  The result
##    is similar to that found in:
##
##    Hanks, E.; Hooten, M.; Johnson, D. & Sterling, J. Velocity-Based
##    Movement Modeling for Individual and Population Level Inference
##    PLoS ONE, Public Library of Science, 2011, 6, e22795
##

library(mgcv)

fit=gam(z~cha+npp+crw+s(t,by=-d2r),
        weights=rep(1/10,nrow(glm.data)),family="poisson",offset=log(tau),data=glm.data)
summary(fit)

plot(fit)
abline(h=0,col="red")

}

}
