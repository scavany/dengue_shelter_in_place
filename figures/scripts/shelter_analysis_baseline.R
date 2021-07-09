rm #======================================================================================
# Author: Sean Cavany
# project: Outbreak response
# Year: 2019
# 
# Code to analyze the output of the agent based model following space spraying
#
# Files needed:
#======================================================================================
library(lattice)
library(stringr)
library(vioplot)
library(viridis)
library(gridExtra)
library(fields)
library(grDevices)
directory <- ".."
foi.outname <- "foi"
none.folder <- 'shelter_none'

n.sims <- 400
n.outputs <- 2

##analysis - load output
setwd(directory)
load("scripts/shelter_analysis_baseline.RData",verbose=T)

## process timeseries
med.none <- lapply(none.baseline,function(x) apply(x, 1, median))
low.none <- lapply(none.baseline, function(y) apply(y,1,function(x) quantile(x, 0.05)))
upp.none <- lapply(none.baseline, function(y) apply(y,1,function(x) quantile(x, 0.95)))

tiff("figures/Fig1.tif", width = 6*600, height =3*600,
     compression="lzw", res=600,pointsize=10)
par(mar=c(4.1,4.1,2.1,4.1))
start <- as.numeric(as.Date("2000-07-01") - as.Date("1999-12-31"))
end <- as.numeric(as.Date("2005-06-30") - as.Date("1999-12-31"))
tvec <- seq(as.Date("2000-07-01"),as.Date("2005-06-30"),"1 day")
plot(tvec,med.none[[1]][start:end]+med.none[[3]][start:end],type='l',xlab="Season",ylab="Median infections",yaxt="n",lwd=2,
     xaxt="n",bty="n",xaxs="i",yaxs="i",las=1)
axis(1,at=seq(as.Date("2001-01-01"),as.Date("2010-01-01"),"1 year")-0.5,
     labels=paste0("0",0:9,"-",str_pad(1:10,width=2,side="left",pad="0")),
     tick=FALSE)
axis(1,at=seq(as.Date("2000-07-01"),as.Date("2010-07-01"),"1 year")-0.5,
     labels=FALSE)
axis(2)
axis(3,at=seq(as.Date("2001-01-01"),as.Date("2010-01-01"),"1 year")-0.5,
     labels=c("High","Sero. inv.",rep("",2),"Low",rep("",5)),
     tick=FALSE)
par(new=TRUE)
plot(tvec,med.none[[2]][start:end],type='l',xlab="",ylab="",
     yaxt="n",xaxt="n",col="red",lty="dashed",bty="n",lwd=2,xaxs="i",yaxs="i",las=1)
axis(4,col="red",col.axis="red")
mtext(expression(italic("Aedes aegypti")),4,3,col="red")
abline(v=seq(as.Date("2000-07-01"),as.Date("2010-07-01"),"1 year")-0.5,lwd=1.5)
year1 <- seq(as.Date("2000-07-01"),as.Date("2001-06-30"),by="1 day")
year2 <- seq(as.Date("2001-07-01"),as.Date("2002-06-30"),by="1 day")
year3 <- seq(as.Date("2004-07-01"),as.Date("2005-06-30"),by="1 day")
polygon(c(year1,rev(year1)),c(rep(0,length(year1)),rep(1e6,length(year1))),
        col=adjustcolor("green",alpha.f=0.25),bty="n")
polygon(c(year2,rev(year2)),c(rep(0,length(year2)),rep(1e6,length(year2))),
        col=adjustcolor("green",alpha.f=0.25),bty="n")
polygon(c(year3,rev(year3)),c(rep(0,length(year3)),rep(1e6,length(year3))),
        col=adjustcolor("green",alpha.f=0.25),bty="n")
dev.off()
