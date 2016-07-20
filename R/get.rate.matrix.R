get.rate.matrix <-
function(beta.static,beta.grad,stack.static,stack.grad,normalize.gradients=FALSE,grad.point.decreasing=TRUE,directions=4,zero.idx=integer()){

    ##
    ## Inputs:
    ##
    ##  stack.static - a raster stack or raster layer of "location-based" covariates
    ##  stack.grad - a raster stack or raster layer of "gradient based covariates
    ##  normalize.gradients - logical.  If TRUE, then normalize all gradient
    ##        covariates by dividing by the length of the gradient vector at each point
    ##  grad.point.decreasing - logical.  If TRUE, then the gradient covariates are positive
    ##        in the direction of decreasing values of the covariate.
    ##
  p.static=nlayers(stack.static)
  p.crw=0
  if(class(stack.grad)=="RasterLayer" | class(stack.grad)=="RasterStack"){
    p.grad=nlayers(stack.grad)
    stack.gradient=rast.grad(stack.grad)

    if(normalize.gradients){
      lengths=sqrt(stack.gradient$grad.x^2+stack.gradient$grad.y^2)
      stack.gradient$grad.x <- stack.gradient$grad.x/lengths
      stack.gradient$grad.y <- stack.gradient$grad.y/lengths
    }
  }else{
    p.grad=0
  }
  p=p.static+p.crw+p.grad



  if(class(stack.static)=="RasterStack"){
        examplerast=stack.static[[1]]
  }
  if(class(stack.static)=="RasterLayer"){
    examplerast=stack.static
  }

  locs=1:ncell(examplerast)

    ##browser()
  ##
  ## Make X matrix
    ##

    nn=ncell(examplerast)
    logR=Matrix(0,nrow=nn,ncol=nn,sparse=TRUE)
    A=Matrix(0,nrow=nn,ncol=nn,sparse=TRUE)
    B=Matrix(0,nrow=nn,ncol=nn,sparse=TRUE)


    adj=adjacent(examplerast,locs,pairs=TRUE,sorted=TRUE,id=TRUE,directions=directions)
    idx.mot=adj[,2:3]

    if(p.static>1){
        logR[idx.mot]=values(stack.static)[idx.mot[,1],]%*%beta.static
    }
    if(p.static==1){
        logR[idx.mot]=values(stack.static)[idx.mot[,1]]*beta.static
    }

  ##
  ## Get x values for gradiant covariates
    ##

    start.cells=idx.mot[,1]
    adj.cells=idx.mot[,2]

  xy.cell=xyFromCell(examplerast,start.cells)
  xy.adj=xyFromCell(examplerast,adj.cells)
  ## Find normalized vectors to adjacent cells
  v.adj=(xy.adj-xy.cell)/sqrt(apply((xy.cell-xy.adj)^2,1,sum))
  ## dot product for gradient covariates
  if(p.grad>0){
    X.grad=v.adj[,1]*stack.gradient$grad.x[start.cells,]+v.adj[,2]*stack.gradient$grad.y[start.cells,]
    ## Make gradient vectors point toward LOWER values (if specified)
    if(grad.point.decreasing==TRUE){
      X.grad=-X.grad
    }
  }

    if(p.grad>1){
        logR[idx.mot]=logR[idx.mot]+X.grad%*%beta.grad
    }
    if(p.grad==1){
        logR[idx.mot]=logR[idx.mot]+X.grad*beta.grad
    }
    ##
    ## Make Rate Matrix
    ##

    R=logR
    R[idx.mot]=exp(R[idx.mot])

    ## exclude "zero" cells
    if(length(zero.idx)>0){
        R[,zero.idx]=0
    }

    R
}
