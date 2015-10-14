path2ctmc <-
function(xy,t,rast){

    ##
    ## Function to turn a discrete-time continuous-space path into a CTDS path
    ##
    ## xy - Tx2 matrix of x,y locations at T time points
    ## t - time points of the observations
    ## rast - a raster object defining the nodes (grid cells) of the CTDS path
    ##
    
    ## path should be a Tx3 matrix with columns: x,y,t
    path=cbind(xy,t)
    
    delta.t=t[2]-t[1]

    ##
    ## get a first pass at the ctmc path
    ##  - transition times are assumed to be at the first observed time
    ##    in the new cell
    ##  
    xycell=cellFromXY(rast,xy)
    rle.out=rle(xycell)
    ec=rle.out$values
    transition.idx=cumsum(rle.out$lengths)
    transition.times=c(t[1],t[transition.idx])
    rt=transition.times[-1]-transition.times[-length(transition.times)]
    
    ##
    ## finding skipped transitions (diagonals)
    ##
    
    ##    browser()
    ncell=ncell(rast)
    ## adjacency matrix
    A=Matrix(0,nrow=ncell,ncol=ncell,sparse=TRUE)
    adj=adjacent(rast,1:ncell)
    A[adj] <- 1
        
    trans=cbind(ec[-length(ec)],ec[-1])
    probs=A[trans]
    idx=which(probs==0)
    n.0=length(idx)
    cat("\n","Fixing ",n.0," diagonal transitions with linear interpolation","\n")
    ec.full=ec[1:idx[1]]
    rt.full=rt[1:idx[1]]
    rt.adjust=rt
    for(i in 1:n.0){
        cat(i," ")
        idx.current=transition.idx[idx[i]]
        xyt.1=path[idx.current,]
        xyt.2=path[idx.current+1,]
        midpoint=1/2*(xyt.1+xyt.2)
        corner=apply(xyFromCell(rast,ec[idx[i]+0:1]),2,mean)
        which.min.1=which.min(abs(corner[1:2]-xyt.1[1:2]))
        which.min.2=which.min(abs(xyt.2[1:2]-corner[1:2]))
        prop.cell.1=(corner[which.min.1]-xyt.1[which.min.1])/(xyt.2[which.min.1]-xyt.1[which.min.1])
        prop.cell.2=(xyt.2[which.min.2]-corner[which.min.2])/(xyt.2[which.min.2]-xyt.1[which.min.2])
        prop.diag=1-prop.cell.1-prop.cell.2
        ## ## plotting what is going on here:
        ## plot(rbind(xyt.1,xyt.2)[,-3],pch=20,type="b")
        ## points(midpoint[1],midpoint[2])
        ## abline(v=corner[1])
        ## abline(h=corner[2])
        
        ## randomly assign diagonal cell if it is too close
        rc.12=rowColFromCell(rast,ec[idx[i]+0:1])
        cell.diag=cellFromRowCol(rast,rc.12[1,1],rc.12[2,2])
        
        ##
        ## fix embedded chain and residence times
        ##
        if(i!=n.0){
            ec.full=c(ec.full,cell.diag,ec[(idx[i]+1):idx[i+1]])
            rt.full[length(rt.full)]=rt.full[length(rt.full)]-delta.t*prop.diag/2
            rt.adjust[idx[i]+1]=rt.adjust[idx[i]+1]-delta.t*prop.diag/2
            rt.full=c(rt.full,delta.t*prop.diag,rt[(idx[i]+1):idx[i+1]])
        }
    }
    ec.full=c(ec.full,cell.diag,ec[(idx[n.0]+1):length(ec)])
    rt.full[length(rt.full)]=rt.full[length(rt.full)]-delta.t*prop.diag/2
    rt.adjust[idx[i]+1]=rt.adjust[idx[i]+1]-delta.t*prop.diag/2
    rt.full=c(rt.full,delta.t*prop.diag,rt[(idx[i]+1):length(rt)])
    
    list(ec=ec.full,rt=rt.full,trans.times=t[1]+cumsum(rt.full))
}
