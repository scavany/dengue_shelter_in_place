library(data.table)
library(mgcv)
library(grDevices)
library(MASS)

##load data and put in data frame
set.seed(10)
load("./shelter_analysis_sweep_mozprob.RData",verbose=T)
parms.shelter <- fread("march17_sweep_mozprob_parameters.csv")
parms.none <- fread("none_sweep_mozprob_parameters.csv")
output.shelter <- data.table(parms = parms.shelter$MosquitoMoveProbability,
                             output=unlist(sweep.mozprob.output[1]))
output.none <- data.table(parms = parms.none$MosquitoMoveProbability,
                          output=unlist(sweep.mozprob.output[2]))
output.none$ld = as.factor(0)
output.shelter$ld = as.factor(1)
output = rbind(output.none,output.shelter)

## generate GAM fits
## knots=6
gam.shelter <- gam(output~s(parms),data=output.shelter)
gam.none <- gam(output~s(parms),data=output.none)
g = gam(output~s(parms,by=ld),data=output)
# https://stats.stackexchange.com/questions/190348/can-i-use-bootstrapping-to-estimate-the-uncertainty-in-a-maximum-value-of-a-gam
beta = coef(g)
Vb = vcov(g)
mrand = mvrnorm(1e4,beta,Vb)
parms = seq(min(output$parms),max(output$parms),length.out=1e3)
Xp.ld0 = predict(
    g,newdata=data.frame(
          parms=parms,
          ld=as.factor(0)),
    type='lpmatrix')
Xp.ld1 = predict(
    g,newdata=data.frame(
          parms=parms,
          ld=as.factor(1)),
    type='lpmatrix')
pred.ld0 = mrand %*% t(Xp.ld0)
pred.ld1 = mrand %*% t(Xp.ld1)
ratio = pred.ld1 / pred.ld0
ratio.median = apply(ratio,2,median)
ratio.lower = apply(ratio,2,function(x)quantile(x,0.025))
ratio.upper = apply(ratio,2,function(x)quantile(x,0.975))

## Plot
axis.size=1
label.size=1
tiff("../figures/Fig5.tif", res=600,
     width=4152, height=4152/2, compression="lzw",pointsize=10)
par(mfrow=c(1,2))
plot(gam.shelter, residuals=T, rug=F, shift=coef(gam.shelter)[1],
     se=F, cex.axis=axis.size, cex.lab=label.size, select=1, bty="n",
     ylab="Cumulative infections", xlab="Daily mosquito move probability",ylim=c(0,2e5),
     xaxs='i',yaxs='i',xlim=c(0,0.5),lwd=2)
par(new = TRUE)
plot(gam.none, residuals=T, rug=F, shift=coef(gam.none)[1],
     se=F, cex.axis=axis.size, cex.lab=label.size, select=1,
     ylim=c(0,2e5),col="red",xaxt="n",yaxt="n",xlab="",ylab="", bty="n",
     xaxs='i',yaxs='i',xlim=c(0,0.5),lwd=2)
abline(v=0.3, lty="dashed")
mtext("A",side=3,line=0, 
       at=par("usr")[1]+0.05*diff(par("usr")[1:2]),
      cex=1.2)
legend(0.03,50000,c("Lockdown","No lockdown"),col=c("black","red"),lty=1,pch=20,bty="n")
plot(parms,ratio.median,lwd=2,type='l',ylim=c(0.5,1.2),bty="n", xlim=c(0,0.5),
     xlab="Daily mosquito move probability",ylab="Ratio of infections",xaxs='i')
polygon(c(parms, rev(parms)),
        c(ratio.lower,
          rev(ratio.upper)),
        col=adjustcolor("red",0.25),border=FALSE)
abline(h=1.0,lty="dashed")
abline(v=0.3,lty="dashed")
mtext("B",side=3,line=0, 
       at=par("usr")[1]+0.05*diff(par("usr")[1:2]),
      cex=1.2)
dev.off()

