path2ctmc <- function(xy,t,rast,directions=4,zero.idx=integer(),print.iter=FALSE){

    ##
    ## Function to turn a discrete-time continuous-space path into a CTDS path
    ##
    ## xy - Tx2 matrix of x,y locations at T time points
    ## t - time points of the observations
    ## rast - a raster object defining the nodes (grid cells) of the CTDS path
    ##


    if(class(rast)=="RasterStack"){
        rast=rast[[1]]
    }
    
    values(rast) <- 1
    values(rast)[zero.idx] <- 0
    trans=transition(rast,prod,directions=directions)


    
    ncell=ncell(rast)
    

    
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
    }
        
    list(ec=ec,rt=c(rt,current.rt),trans.times=c(t[1]+cumsum(rt),cumsum(rt)[length(rt)]+current.rt))
}
