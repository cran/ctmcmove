ctmc2glm <-
function(ctmc,stack.static,stack.grad,crw=TRUE,normalize.gradients=FALSE,grad.point.decreasing=TRUE,include.cell.locations=TRUE,directions=4){

    ## Function to take a CTMC path and covariate rasters and return data
    ##  that can be fit as a Poisson GLM
    ##
    ## Inputs:
    ##
    ##  ctmc - a "ctmc" object (output from "path2ctmc.r")
    ##  stack.static - a raster stack or raster layer of "location-based" covariates
    ##  stack.grad - a raster stack or raster layer of "gradient based covariates
    ##  crw - logical.  If TRUE, then an autocovariate is computed defined as a
    ##        directional covariate pointed in the direction of the most recent
    ##        transition in the CTMC
    ##  normalize.gradients - logical.  If TRUE, then normalize all gradient
    ##        covariates by dividing by the length of the gradient vector at each point
    ##  grad.point.decreasing - logical.  If TRUE, then the gradient covariates are positive
    ##        in the direction of decreasing values of the covariate.
    ##
  p.static=nlayers(stack.static)
  p.crw=0
  if(crw){
    p.crw=1
  }
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

  locs=ctmc$ec
  wait.times=ctmc$rt
    
  ##
  ## Make X matrix 
  ##

  ## sort.idx=sort(locs,index.return=TRUE)$ix
  ## This is for a rook's neighborhood
##  n.nbrs=4
  ## sort.idx=rep(sort.idx,each=n.nbrs)
  adj=adjacent(examplerast,locs,pairs=TRUE,sorted=TRUE,id=TRUE,directions=directions)
  adj.cells=adj[,3]
  rr=rle(adj[,1])
    time.idx=rep(rr$values,times=rr$lengths)
  start.cells=adj[,2]
  z=rep(0,length(start.cells))
  idx.move=rep(0,length(z))
  diag.move=rep(0,length(locs))
  for(i in 1:(length(locs))){
    idx.t=which(time.idx==i)
    idx.m=which(adj.cells[idx.t]==locs[i+1])
    z[idx.t[idx.m]] <- 1
    if(length(idx.m)==0){
      diag.move[i]=1
    }
  }

  #browser()

  
  ## Tau
  tau=rep(wait.times,times=rr$lengths)
  ##
  t=rep(ctmc$trans.times,times=rr$lengths)
      
  ##
  ## Get x values for static covariates
  ##

  if(nlayers(stack.static)>1){
    X.static=values(stack.static)[start.cells,]
  }else{
    X.static=matrix(values(stack.static)[start.cells],ncol=1)
  }
  #colnames(X.static) <- layerNames(stack.static)
  colnames(X.static) <- names(stack.static)
  
  ##
  ## Get x values for gradiant covariates
  ##

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
    colnames(X.grad) <- colnames(stack.gradient$grad.x)
  }
  

  ##
  ## Get crw covariate
  ##

  #browser()

  ##
  ## should probably figure a better way to handle diagonals.
  ##

    ## update 20160519
        ##     idx.move=rep(,length(locs))
        ## idx.move[-which(diag.move==1)] <- which(z==1)
    idx.move=which(z==1)
    ## last point
    idx.move=c(idx.move,length(z))

    
  v.moves=v.adj[rep(idx.move[1:(length(rr$lengths)-1)],times=rr$lengths[-1]),]
  ## shift v.moves to be a lag 1 direction
  v.moves=rbind(matrix(0,ncol=2,nrow=rr$lengths[1]),v.moves)
  ## dot product with vectors to adjacent cells
  #browser()
  X.crw=apply(v.moves*v.adj,1,sum)

  #browser()
 
  
  ## Compiling Matrices
  if(crw==FALSE & p.grad>0){
    X=cbind(X.static,X.grad)
  }
  if(crw==TRUE & p.grad>0){
    X=cbind(X.static,X.grad,X.crw)
    colnames(X)[ncol(X)]="crw"
  }
  if(crw==FALSE & p.grad==0){
    X=cbind(X.static)
  }
  if(crw==TRUE & p.grad==0){
    X=cbind(X.static,X.crw)
    colnames(X)[ncol(X)]="crw"
  }
  
  ## browser()

    if(include.cell.locations){
        xys=cbind(xy.cell,xy.adj)
        colnames(xys)=c("x.current","y.current","x.adj","y.adj")
        X=cbind(X,xys)
    }
                

  T=length(wait.times)
  p=ncol(X)
  
  #browser()
    
  

  out=data.frame(z=z,X,tau=tau,t=t)
  ## remove last time step
    T=nrow(out)
    out=out[-((T-3):T),]
    out
}
