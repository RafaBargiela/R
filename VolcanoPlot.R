##########################################################################################
# Title: Volcano plot script for ASV abundance for Nat Med manuscript
# Author: Dr. Rafael Bargiela
# LAST UPDATE: 5th of May, 2023                                  
##########################################################################################
# Change below the working directory
setwd("")
library(colorspace)

  File<-"" # Add here the path to the file with your data
  # Data table will be loaded in the PIstats value.
  # It is set for a tab-separated file, with header and rownames on the first column (ASV IDs, I guess)
  # I assume you use "." for decimals
  # Then, you need a column for the fold-change and another for the p-values
  PIstats<-read.table(file=File,row.names=1,header=TRUE,dec=".",sep="\t",check.names=FALSE)
  OUTPUT<-"Volcano.jpg"
  ## 1. Volcano plot ####  
    ### 1.1. Setting up the variables ---------------
      FC<-log(PIstats[,1],base=2) # Fold change
      pv<--log(PIstats[,2],base=10) # p-values
      names(FC)<-rownames(PIstats)
      names(pv)<-rownames(PIstats)
      pv<-sort(pv,decreasing=FALSE)
      FC<-FC[names(pv)]
    ### 1.2. Graphical parameters ------------------
      h<-2000
      w<-2000
      r<-300
      ylim<-c(0,45)
      xlim<-c(-30,30)
      xlab<-expression(bold("Fold Change of log"["2"]~"(ASV abundance)"))
      ylab<-expression(bold("-log"["10"]~"(p-values)"))
      vertical.lines.x<-c(-2.5,2.5)
      horiz.lines.y<-1.25  
      pv.labels.threshold<-10
      FC.labels.threshold<-20
      cols<-hcl.colors(4,"plasma")
      # Values below are for the correction of the labels possition
      # All row labels (taxonomic groups names) are added, however, only those reaching some parameters will be labelled on the volcano
      LABELS<-TRUE # change to false if you don't want any labels and add them by yourself manually
      ycorr<-vector("numeric",length(FC))
      names(ycorr)<-names(FC)
      xcorr<-rep(0.15,length(FC)) # By default I set up 0.15 of separation from the dots, change as needed
      names(xcorr)<-names(FC)
      # NOTE: In case labels overlapped, uncomment lines below and modify individually each label on the x or y axis
      # Mark the specific label name on the xcorr or ycorr vectors and add additional correction
      # ycorr[c("")]<-c()
      # xcorr[""]<-c()
    ## Plotting Volcano Plot  
      jpeg(OUTPUT,height=h,width=w,res=r,quality=100)
        par(mar=c(5,4,4,2),mgp=c(2.5,1,0))
        plot(NA,NA,xlim=xlim,ylim=ylim,xlab=list(xlab,cex=1.2,font=2),
             ylab=list(ylab,font=2,cex=1.2),axes=FALSE)
          axis(1,seq(xlim[1],xlim[2],10),labels=TRUE,tick=TRUE,lwd.ticks=2,las=1,cex.axis=1.2,font.axis=2)
          yi<-ylim[2]+((ylim[2]*15)/100)
          yw<-(yi-0)/100
          rect(xlim[1],yi,0,0,border=NA,col=rgb(t(col2rgb("thistle2")),alpha=100,maxColorValue=255),xpd=TRUE)
          rect(xlim[2],yi,0,0,border=NA,col=rgb(t(col2rgb("grey95")),alpha=100,maxColorValue=255),xpd=TRUE)
        axis(2,seq(ylim[1],ylim[2],1),labels=TRUE,tick=TRUE,lwd.ticks=2,las=1,cex.axis=1.2,font.axis=2)
        arrows(vertical.lines.x,ylim[1],vertical.lines.x,ylim[2],code=0,lty=3,lwd=1,col="grey50")    
        arrows(xlim[1],horiz.lines.y,xlim[2],horiz.lines.y,code=0,lty=3,lwd=1,col="grey50")
        box(lwd=2)
        
          for(p in 1:length(FC)){
            col<-1
            if(abs(FC[p])>FC.labels.threshold){col<-col+1}
            if(abs(pv[p])>=pv.labels.threshold){col<-col+1}
            points(FC[p],pv[p],pch=19,xpd=TRUE,col=cols[col],cex=1.2) 
            if(col==3){
              print(names(FC)[p])
              if(LABELS==TRUE){
                text(FC[p]+xcorr[p],pv[p]+ycorr[p],labels=names(FC)[p],adj=c(0,0.5),cex=0.8,xpd=TRUE,font=2)
              }
            }
          }
        
        text(xlim[1],ylim[2]+((ylim[2]*10)/100),labels="Overrepresented on\nHSIL",cex=1.2,font=2,adj=c(0,0.5),xpd=TRUE,col="grey50")
        text(xlim[2],ylim[2]+((ylim[2]*10)/100),labels="Overrepresented on\nnon-HSIL",cex=1.2,font=2,adj=c(1,0.5),xpd=TRUE,col="grey50")
        
      dev.off()
