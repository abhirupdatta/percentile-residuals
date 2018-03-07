## gets traditional residuals R^* using samples from F_w and F
## datavec: ndata x 1 vector of data points
## predmat ndata x nsamples matrix with i^th row giving 
## predictive distribution for the i^th data point
get.rstar <- function(datavec, predmat) {
    til.mu  <- rowMeans(predmat);
    til.sig <- apply(predmat,1,sd)
    rst     <- (datavec - til.mu)/til.sig;
}

## gets percentile based residuals R-dag
## datavec: ndata x 1 vector of data points
## predmat ndata x nsamples matrix with i^th row giving 
## predictive distribution for the i^th data point
get.rdagger <- function(datavec, predmat){
    #eps=1e-32
    #qnorm(ecdf(c(rep(smp.Fw,1),-1/eps,1/eps))(smp.F))
    n=length(datavec)
    sapply(1:n,function(i) qnorm(ecdf(predmat[i,])(datavec[i])))
}

### roc curve ###
## datamat: ndata x nsamples vector of data points
## truemat ndata x nsamples matrix with i^th row giving 
## true distribution for the i^th data point
## workmat ndata x nsamples matrix with i^th row giving 
## working distribution for the i^th data point
### txt: name of the file where the plots are ging to be stored
### save: flag whether to save plot or display it
### path: folder for saving the plot
### onlyroc: flag, it true plots only the roc, otherwise plots both null and alternate distributions and roc
### newimg: flag, opens a new image window
### xlim: specifies the x-axis limits for the roc
### alt: type of alternate hypotheses -- "both", "left" or "right"
roc=function(datamat, truemat, workmat,txt,save,path,onlyroc,newimg,xlim,alt){
  
    ##under H0, Fw = F
    rstar0   <- apply(datamat,2,get.rstar, truemat);
    rdagger0 <- apply(datamat,2,get.rdagger, truemat);
    
    ##under H1, Fw <> F
    rstar1   <- apply(datamat,2,get.rstar, workmat);
    rdagger1 <- apply(datamat,2,get.rdagger, workmat);
    
    ### average over the replicate datasets 
    fprstar=rowMeans(apply(rstar0,2,power,alt))
    powstar=rowMeans(apply(rstar1,2,power,alt))

    fprdagger=rowMeans(apply(rdagger0,2,power,alt))
    powdagger=rowMeans(apply(rdagger1,2,power,alt))
        
    dstar0=density(as.vector(rstar0))
    dstar1=density(as.vector(rstar1))
    ddag0=density(as.vector(rdagger0))
    ddag1=density(as.vector(rdagger1))
    
    if(save=="y"){
    fname=paste0(path,txt)
    #pdf(fname,encoding="MacRoman")
    png(fname,height=1024,width=768)
    }else{
        if(newimg) dev.new()
        }
    
    if(onlyroc){
        title=paste(unlist(strsplit(txt,".png")))
        plot(fprstar,powstar,xlab="FPR",ylab="Power",type="l",col="red",
             lwd=3,cex.lab=1.5,main=title,cex.main=1,xlim=xlim)
        lines(fprdagger,powdagger,col="blue",lwd=3)
        lines(0:100/100,0:100/100,lty=2,col="grey")
    }
    else{
        par(mfrow=c(2,2))
        
        plot(dstar0,xlab='x',ylab="density",xlim=range(-3.5,dstar0$x,ddag0$x,3.5),main="Null",
             col="red",type='l',lwd=3,ylim=range(dnorm(0),dstar0$y,ddag0$y),cex.lab=1.5,cex.main=2)
        lines(ddag0,col="blue",lwd=3)
        lines(seq(-3.5,3.5,0.01),dnorm(seq(-3.5,3.5,0.01)),lwd=2)
        
        plot(dstar1,xlab='x',ylab="density",xlim=range(-3.5,dstar1$x,ddag1$x,3.5),main="Alt.",
             col="red",type='l',lwd=3,ylim=range(dnorm(0),dstar1$y,ddag1$y),cex.lab=1.5,cex.main=2)
        lines(ddag1,col="blue",lwd=3)
        lines(seq(-3.5,3.5,0.01),dnorm(seq(-3.5,3.5,0.01)),lwd=2)
        
        plot(fprstar,powstar,xlab="FPR",ylab="Power",type="l",col="red",lwd=3,cex.lab=1.5,main="ROC",
             cex.main=2,xlim=xlim)
        lines(fprdagger,powdagger,col="blue",lwd=3)
        lines(0:100/100,0:100/100,lty=2,col="grey")
        
        plot(1,col="white",xaxt='n',yaxt='n',ann=FALSE,bty='n')
        l.txt <- c(as.expression(bquote(R^~"*")),
                   as.expression(bquote(R^~"\u2020")),"N(0,1)");
        l.col <- c("red", "blue","black");
        legend("center", bty = "n", legend = l.txt, col = l.col, lty = 1,cex=3,lwd=3)
        }
    
    if(save=="y")  dev.off()
    cbind(rstar0,rstar1,rdagger0,rdagger1)
    
}

power=function(samples,alt){
  z=qnorm(seq(0.5,1,length=1000))
  if(alt=="both"){
    
    1-ecdf(abs(samples))(z)
    #powstar=1-ecdf(abs(rstar1))(z)
    
    #fprdagger=1-ecdf(abs(rdagger0))(z)
    #powdagger=1-ecdf(abs(rdagger1))(z)
    }else if(alt=="right"){
      1-ecdf(samples)(z)
      }else{
        ecdf(samples)(-z)
        }
  }

### generates random samples ###
### type = current distributions allowed: normal, gamma and mixture normal
### for normal and gamma parameter argument is a 2x1 vector
### for a k-component mixnorm it is a matrix kx3 matrix where first column gives weights, second gives means and third gives sds
randomgen=function(type,par,N,seed,partype){
  set.seed(seed)
  if(type=="normal"){
    rnorm(N,par[1],par[2])
  }else  if(type=="gamma"){
    if(partype==0){
        rgamma(N,par[1],par[2])
    }else{
        beta=par[1]/par[2]
        #print(beta)
        rgamma(N,par[1],beta)
        }
  }else  if(type=="mixnorm"){
    components <- sample(1:nrow(par),prob=par[,1],size=N,replace=TRUE)
    rnorm(N,mean=par[components,2],sd=par[components,3])    
  }
}

#### MAIN FUNCTION ####
#### plots densties of R* and Rdag under null and alt, also plots ROC curve
#### Ftype=true distribution (currently supports "normal","gamma","mixnorm")
#### par= parameters of the true dist. 
#### par for Ftype normal Ftype should be a 2x1 vector with mean and sd
#### par for Ftype gamma should be a 2x1 vector with shape and rate
### partype additional optional parameter for gamma distribution
### partype=1 implies the gamma distribution is parametrized in terms of its shape and mean
### partype=0 implies it is parametrized in terms of shape and rate
#### par for Ftype mixnorm Ftype should be a px3 matrix with 
# the columns corresponding to mixture weights, mean and sd respectively
#### Fwtype=working distribution (currently supports "normal","gamma","mixnorm")
####  parw=parameters of the working dist. (same format as par)
#### N,Nw are sample sizes for null and alt dist. (defaults 10000)
#### seed=random seed (default 1)
#### save=saves the plot (if "y") or shows the plot in a popup window (if "n")
#### path=folder path if save=="y" (defaults to save in the working directory)
### xlim: specifies the x-axis limits for the roc
### alt: type of alternate hypotheses -- "both", "left" or "right"
simulations=function(Ftype,par,Fwtype,parw,ndata=100,nsamples=500,nseed=100,txt=NULL,
                     save="n",path="",partype=1,onlyroc=F,newimg=T,xlim=c(0,1),alt="both"){
  txt=paste0(paste(Ftype,paste(round(par,2),collapse="_"),
                   Fwtype,paste(round(parw,2),collapse="_"),alt,sep="_"),".png")
  print(txt)
  #print(seed)
  #print(save)
  datamat=sapply(1:nseed,function(i) randomgen(Ftype,par,ndata,i,partype))
  truemat=sapply(nseed+(1:nsamples),function(i) randomgen(Ftype,par,ndata,i,partype))
  predmat=sapply(nseed+(1:nsamples),function(i) randomgen(Fwtype,parw,ndata,i,partype))
  roc(datamat,truemat,predmat,txt,save,path,onlyroc,newimg,xlim,alt)
}


