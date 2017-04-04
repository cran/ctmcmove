path2ctmc <- function(xy,t,rast,directions=4,zero.idx=integer(),print.iter=FALSE,method="ShortestPath"){

    ##
    ## Function to turn a discrete-time continuous-space path into a CTDS path
    ##
    ## xy - Tx2 matrix of x,y locations at T time points
    ## t - time points of the observations
    ## rast - a raster object defining the nodes (grid cells) of the CTDS path
    ## method = either "ShortestPath" or "LinearInterp"


    if(class(rast)=="RasterStack"){
        rast=rast[[1]]
    }
    
    values(rast) <- 1
    values(rast)[zero.idx] <- 0
    trans=transition(rast,prod,directions=directions)


    
    ncell=ncell(rast)

    if(method=="LinearInterp"){
        ## adjacency matrix
        A=Matrix(0,nrow=ncell,ncol=ncell,sparse=TRUE)
        adj=adjacent(rast,1:ncell)
        A[adj] <- 1
    }
    

    
    ## path should be a Tx3 matrix with columns: x,y,t
    path=cbind(xy,t)

    ## check to make sure t is in time order
    tidx=sort(t,index.return=T)$ix
    path=path[tidx,]

    
    T=nrow(path)

    ec.all=cellFromXY(rast,xy)

    ##head(cbind(path,ec.all))
    
    ec=ec.all[1]
    current.cell=ec
    rt=integer()
    current.rt=0
    if(print.iter){
        cat("Total locations =",T,"\n")
    }
    for(i in 2:(T)){
        if(print.iter){
            cat(i," ")
        }
        if(ec.all[i]==current.cell){
            current.rt=current.rt+(t[i]-t[i-1])
        }else{
            if(method=="ShortestPath"){
                sl=shortestPath(trans,as.numeric(xy[i-1,]),as.numeric(xy[i,]),"SpatialLines")
                slc=coordinates(sl)[[1]][[1]]
                ## get cells between start and end cell.
                ## this includes both the end cell, and the first cell
                sl.cells=cellFromXY(rast,slc)
                ## evenly divide time spent between x-y locations into all cells visited
                t.in.each.cell=(t[i]-t[i-1])/length(sl.cells)
                ## add time spent in current cell
                current.rt=current.rt+(t.in.each.cell)
                rt=c(rt,current.rt)
                ## add residence times in new cells visited
                current.rt=t.in.each.cell
                if(length(sl.cells)>2){
                    rt=c(rt,rep(t.in.each.cell,length(sl.cells)-2))
                }
                ec=c(ec,sl.cells[-1])
                current.cell=ec[length(ec)]
            }
            if(method=="LinearInterp"){
                if(A[current.cell,ec.all[i]]==1){
                    ##    cat("trans ")
                    rt=c(rt,current.rt+(t[i]-t[i-1]))
                    current.rt=0
                    ec=c(ec,ec.all[i])
                    current.cell=ec.all[i]
                }else{
                    ##    cat("impute ")
                    xyt.1=path[i-1,]
                    xyt.2=path[i,]
                    d=sqrt(sum((xyt.1[-3]-xyt.2[-3])^2))
                    rast.res=res(rast)[1]
                    ## linearly interpolate
                    xapprox=approx(c(xyt.1[3],xyt.2[3]),c(xyt.1[1],xyt.2[1]),n=max(100,round(d/rast.res*100)))
                    yapprox=approx(c(xyt.1[3],xyt.2[3]),c(xyt.1[2],xyt.2[2]),n=max(100,round(d/rast.res*100)))
                    tapprox=xapprox$x
                    xapprox=xapprox$y
                    yapprox=yapprox$y
                    ##
                    xycell.approx=cellFromXY(rast,cbind(xapprox,yapprox))
                    rle.out=rle(xycell.approx)
                    ec.approx=rle.out$values
                    transition.idx.approx=cumsum(rle.out$lengths)
                    transition.times.approx=tapprox[transition.idx.approx]
                    rt.approx=transition.times.approx[-1]-transition.times.approx[-length(transition.times.approx)]
                    rt=c(rt,current.rt+transition.times.approx[1]-as.numeric(xyt.1[3]),rt.approx[-length(rt.approx)])
                    current.rt=rt.approx[length(rt.approx)]
                    ec=c(ec,ec.approx[-1])
                    current.cell=ec[length(ec)]
                }
            }
        }
    }
    list(ec=ec,rt=c(rt,current.rt),trans.times=c(t[1]+cumsum(rt),cumsum(rt)[length(rt)]+current.rt))
}
