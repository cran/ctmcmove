get.rate.matrix <-
  function(object, stack.static, stack.grad, normalize.gradients=FALSE, grad.point.decreasing=TRUE, directions=4, zero.idx=integer(), coef){
    
    ##
    ## Inputs:
    ##
    ##  object - A fitted GLM or GAM (from mgcv) object
    ##  coef - a coefficient vector to be used other than what is provided in object
    ##  stack.static - a raster stack or raster layer of "location-based" covariates
    ##  stack.grad - a raster stack or raster layer of "gradient based covariates
    ##  normalize.gradients - logical.  If TRUE, then normalize all gradient
    ##        covariates by dividing by the length of the gradient vector at each point
    ##  grad.point.decreasing - logical.  If TRUE, then the gradient covariates are positive
    ##        in the direction of decreasing values of the covariate.
    ##
    if(inherits(stack.static, "Raster")){
      p.static=nlayers(stack.static)
    } else p.static = 0
    p.crw=0
    if(inherits(stack.grad, "Raster")){
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
    nn=ncell(examplerast)
    R=Matrix(0,nrow=nn,ncol=nn,sparse=TRUE)
    adj=adjacent(examplerast,locs,pairs=TRUE,sorted=TRUE,id=TRUE,directions=directions)
    idx.mot=adj[,2:3]
    
    if(p.static>0){
      X = data.frame(values(stack.static)[idx.mot[,1],])
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
      X = cbind(X,X.grad)
    }
    
    X$tau = 1
      X$crw=0
    
    if(missing(coef)){
      R[idx.mot] = as.vector(predict(object, newdata=X, type="response"))
    } else{
      if(length(coef) != length(coefficients(object))) stop("'coef' vector is not the correct length!")
      R[idx.mot] = exp(X%*%coef)
    }
    
    ## exclude "zero" cells
    if(length(zero.idx)>0){
      R[zero.idx,zero.idx]=0
    }
    
    R
  }
