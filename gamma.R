### runs simulations for a bunch of alternate gamma dists (true being gamma (1,1)) ###
### All GAMMA distributions here use the second param as the mean ####
source("percentile_residuals.r")

shapevec=meanvec=c(0.2,0.5,1,2,5)
l=length(shapevec)
settings = as.matrix(expand.grid(shapevec,meanvec))

dev.new()
par(mfrow=c(l,l))
s=apply(settings,1,function(par){
    if(par[1]==1 & par[2]==1){
        plot(1,col="white",xaxt='n',yaxt='n',ann=FALSE,bty='n')
        l.txt <- c(as.expression(bquote(R^~"*")),
                   as.expression(bquote(R^~"\u2020")),"N(0,1)");
        l.col <- c("red", "blue","black");
        legend("center", bty = "n", legend = l.txt, col = l.col, lty = 1,cex=2,lwd=3)
        }else{
        s=simulations("gamma",c(1,1),"gamma",c(par[1],par[2]),onlyroc=T,newimg=F)
        }
    })


### same runs with ROC curve focussing on only upto 10% FPR ###
dev.new()
par(mfrow=c(l,l))
s10=apply(settings,1,function(par){
    if(par[1]==1 & par[2]==1){
        plot(1,col="white",xaxt='n',yaxt='n',ann=FALSE,bty='n')
        l.txt <- c(as.expression(bquote(R^~"*")),
                   as.expression(bquote(R^~"\u2020")),"N(0,1)");
        l.col <- c("red", "blue","black");
        legend("center", bty = "n", legend = l.txt, col = l.col, lty = 1,cex=2,lwd=3)
    }else{
        s=simulations("gamma",c(1,1),"gamma",c(par[1],par[2]),onlyroc=T,newimg=F,xlim=c(0,0.1))
    }
})


### runs with right sided alternate hypotheses ###
dev.new()
par(mfrow=c(l,l))
sright=apply(settings,1,function(par){
    if(par[1]==1 & par[2]==1){
        plot(1,col="white",xaxt='n',yaxt='n',ann=FALSE,bty='n')
        l.txt <- c(as.expression(bquote(R^~"*")),
                   as.expression(bquote(R^~"\u2020")),"N(0,1)");
        l.col <- c("red", "blue","black");
        legend("center", bty = "n", legend = l.txt, col = l.col, lty = 1,cex=2,lwd=3)
    }else{
        simulations("gamma",c(1,1),"gamma",c(par[1],par[2]),onlyroc=T,newimg=F,xlim=c(0,0.1),alt="right")
    }
})


### runs with left sided alternate hypotheses ###
dev.new()
par(mfrow=c(l,l))
sleft=apply(settings,1,function(par){
    if(par[1]==1 & par[2]==1){
        plot(1,col="white",xaxt='n',yaxt='n',ann=FALSE,bty='n')
        l.txt <- c(as.expression(bquote(R^~"*")),
                   as.expression(bquote(R^~"\u2020")),"N(0,1)");
        l.col <- c("red", "blue","black");
        legend("center", bty = "n", legend = l.txt, col = l.col, lty = 1,cex=2,lwd=3)
    }else{
        simulations("gamma",c(1,1),"gamma",c(par[1],par[2]),onlyroc=T,newimg=F,xlim=c(0,0.1),alt="left")
    }
})
