#======================================================================================
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
directory <- ".."
foi.outname <- "foi"
folder <- "shelter_monthly_none_vc"
## none.folder <- 'shelter_none'
## analysis.folders <- c("shelter_march17")
## monthly.folders <- c("shelter_01low","shelter_01high")
## all.folder <- "shelter_all"
## suffixes <- c("/", "_2/")

n.sims <- 400
n.batches <- 1
## yearly.cats <- 8
## none.cats <- 4
monthly.cats <- 12
## n.cats <- c(yearly.cats)
## names(n.cats) <- analysis.folders
n.outputs <- 1

##analysis - load output
setwd(directory)
load("scripts/shelter_analysis_ts.RData",verbose=T)
load("scripts/shelter_analysis_extra_ts.RData",verbose=T)
monthly.folders <- names(mon.extra.ts)

## process timeseries
med.mon.ts <- lapply(mon.ts,
                     function(x) apply(x, c(1,3), median))
low.mon.ts <- lapply(mon.ts,
                     function(x) apply(x, c(1,3),
                                       function(y) quantile(y, 0.05,na.rm=TRUE)))
upp.mon.ts <- lapply(mon.ts,
                     function(x) apply(x, c(1,3),
                                       function(y) quantile(y, 0.95,na.rm=TRUE)))

med.mon.extra.ts <- lapply(mon.extra.ts,
                     function(x) apply(x, c(1,3), median))
low.mon.extra.ts <- lapply(mon.extra.ts,
                     function(x) apply(x, c(1,3),
                                       function(y) quantile(y, 0.05)))
upp.mon.extra.ts <- lapply(mon.extra.ts,
                     function(x) apply(x, c(1,3),
                                       function(y) quantile(y, 0.95)))


month.startday.vec <- c(c(182,213,244,274,305,335),365+c(1,32,60,91,121,152))-182

##proportional increase for Vector control vs lockdown
start <- 182 + 365*1
end <- 181 + 365*3

ratios <- matrix(rep(NA,1000*12),nrow=1000)
for (i in 1:1000) {
    print(i)
    indices <- sample(1:400,100,replace=TRUE)
    meanA <- apply(apply(mon.ts[[2]][start:end,indices,],c(2,3),sum),2,mean)
    meanB <- apply(apply(mon.extra.ts[[folder]][start:end,indices,],c(2,3),sum),2,mean)
    ratios[i,] <- meanA/meanB
}

## Plot vector control as ratio.
tiff("figures/Fig4.tif",
     width = 6*600, height =7.4*600,
     compression="lzw", res=600)
par(mfrow=c(4,3),
    mar = 0.1+c(0,0,0,0),
    oma = c(4,4,1,1))
cols <- viridis(5)
start <- 182 + 365*1
end <- 181 + 365*3
tvec <- seq(as.Date("2001-07-01"),by="1 day",length.out=end-start+1)
cols.alpha <- viridis(5, alpha=c(0.3, 0.3))
panels <- c("Jul","Aug","Sep","Oct","Nov","Dec",
            "Jan","Feb","Mar","Apr","May","Jun")
par(mar = 0.1+c(1,1,0,0))
for (i in 1:12){
    shelter.start <- 1 + month.startday.vec[i]
    shelter.end <- shelter.start + floor(3*365/12 + 0.5)
    vector.control <- rollmean(med.mon.extra.ts[[folder]][,i],31,na.pad=T,align="center",fill=0)[start:end]
    lockdown <- rollmean(med.mon.ts[[2]][,i],31,na.pad=T,align="center")[start:end]
    vector.control[vector.control < 1e-10]  <- NA
    ratio <- lockdown/vector.control
    if (i %in% c(1,4,7,10)) {
        plot(tvec,ratio,
             type='l',log="y",
             ## ylim=c(0,1.1*max(upp.mon.ts[[2]][start:end])),
             xaxt="n",xlab="",ylab="",ylim=c(0.2,20),
             lwd=1)
        mtext("Ratio of infections",2,3)
    } else {
        plot(tvec,ratio,
             type='l', log="y",ylim=c(0.2,20),
             ## ylim=c(0,1.1*max(upp.mon.ts[[2]][(start):end])),
             xaxt="n",xlab="",yaxt="n",ylab="",lwd=1)
    }
    if (i <=6) {
        text(as.Date("2003-04-01"),
             18,panels[i],
             cex=1.2)
    } else {
        text(as.Date("2001-09-01"),
             18,panels[i],
             cex=1.2)
    }
    if (i %in% 10:12) {
        mtext("Season",1,3)
        axis(1,at=seq(as.Date(paste0("2002-01-01")),by="1 year",length.out=2),
             tick=F,labels=c("Sero. inv.",expression(italic("following"))))
    }
    axis(1,at=seq(as.Date(paste0("2001-07-01")),by="1 year",length.out=3),
         labels=F)
    polygon.indexes <- !is.na(ratio)
    polygon(c(tvec[shelter.start:shelter.end],rev(tvec[shelter.start:shelter.end])),
            c(rep(0.01,length(tvec[shelter.start:shelter.end])),
              rep(100,length(tvec[shelter.start:shelter.end]))),
            col=adjustcolor("gray",0.375),border=FALSE)
    polygon(c(tvec[polygon.indexes],rev(tvec[polygon.indexes])),
            c(pmax(1,ratio[polygon.indexes]),rep(1,sum(polygon.indexes))),
            col=adjustcolor("green",0.25),border=FALSE)
    polygon(c(tvec[polygon.indexes],rev(tvec[polygon.indexes])),
            c(pmin(1,ratio[polygon.indexes]),rep(1,length(tvec[polygon.indexes]))),
            col=adjustcolor("red",0.25),border=FALSE)
    ## abline(v=tvec[shelter.start],lty="dashed")
    ## abline(v=tvec[shelter.end],lty="dotted")
    if (grepl("twice", folder)) {
            abline(v=tvec[shelter.start]+365,lty="dashed")
            abline(v=tvec[shelter.end]+365,lty="dotted")
    }
    abline(h=1)
}
dev.off()
