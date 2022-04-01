############################################################################################
###   Diversity indexes and rarefaction curves  ####### last update: 16th of September, 2021
############################################################################################

DiversityRichnessAnalysis<-function(M,prediction=FALSE){
  require(vegan)
  DRA<-{}
  M<-as.matrix(M)
  
  ### Diversity indexes ####
  ShannonI<-diversity(t(M),index="shannon") 
  ShannonI<-format(ShannonI,nsmall=2,digits=2)
  DRA$Shannon<-ShannonI
  # Rarefaction curves ####
  Curves<-vector("list")   # Loop were several rarefaction curves are calculated and store in a vector list
  Subsample<-vector("list")
  for(sample in 1:ncol(M)){
    v<-as.numeric(M[,sample])
    curve<-rarefy(v,sample=seq(0,sum(v),by=50))
    Curves[[colnames(M)[sample]]]<-curve[1,]
    Subsample[[colnames(M)[sample]]]<-attributes(curve)$Subsample  

  }
  DRA$rarefaction.curve<-Curves
  DRA$rarefaction.curve.subsample<-Subsample

  if(prediction==TRUE){
      require(drc)
      Predictions<-vector("list")
      max<-max(unlist(DRA$rarefaction.curve.subsample))
          for(sample in 1:length(Curves)){
              # For prediction of additional data on the curve
              if(Subsample[[sample]][length(Subsample[[sample]])] < max){
                    dat<-data.frame(curve=Curves[[sample]],subsample=Subsample[[sample]]) # creating a data.frame
                    mod<-drm(curve~subsample,data=dat,fct=W1.4())               
                    newdata<-seq(dat$subsample[length(dat$subsample)],max,by=50)
                    pre<-as.data.frame(predict(mod,newdata=data.frame(subsample=newdata),interval="predict"))   
              }else{
                newdata<-NA
                pre<-NA
              }
              prediction<-cbind(newdata,pre)
              Predictions[[colnames(M)[sample]]]<-prediction
            }
      DRA$rarefaction.curve.prediction<-Predictions
      message("\nCurve predictions have been assessed\n")
  }

  return(DRA)
}

printRarefaction<-function(Curvesx,Curvesy,xlim=c(0,15000),xlim.labs=seq(0,xlim[2],by=7500),ylim=c(0,400),cex.axis=1.2,lwd=2,cols="black",tittle=""){
  max<-max(unlist(Curvesx))
  plot(NA,NA,las=1,xlab="",ylab="",xlim=xlim,ylim=ylim,axes=FALSE)
  axis(2,seq(ylim[1],ylim[2],by=100),labels=TRUE,las=1,tick=TRUE,lwd.ticks=2,cex.axis=cex.axis,font.axis=2,lwd=2)
  axis(1,xlim.labs,labels=TRUE,las=1,tick=TRUE,lwd.ticks=2,cex.axis=cex.axis,font.axis=2,lwd=2)
  arrows(0,seq(ylim[1],ylim[2],by=100),xlim[2],seq(ylim[1],ylim[2],by=100),lwd=1,lty=3,col="grey80",code=0,xpd=TRUE)
  text(xlim[1],ylim[2]-((ylim[2]*10)/100),labels=tittle,cex=1.2,font=2,adj=c(0,0.5),xpd=TRUE)
  box(lwd=2)
  for(curve in 1:length(Curvesy)){
    if(length(lwd)>1){
      width<-lwd[curve]
    }else{
      width<-2
    }
    points(Curvesx[[curve]],Curvesy[[curve]],type="l",lwd=width+1,xpd=TRUE,col="black")        
    points(Curvesx[[curve]],Curvesy[[curve]],type="l",lwd=width,xpd=TRUE,col=cols[curve],lty=1)
     if(Curvesx[[curve]][length(Curvesx[[curve]])] < max){

      dat<-data.frame(curve=Curvesy[[curve]],subsample=Curvesx[[curve]]) # creating a data.frame
      mod<-drm(curve~subsample,data=dat,fct=W1.4())      # Calculating model         
      newdata<-seq(dat$subsample[length(dat$subsample)],max,by=50) # data frame for new data
      pre<-as.data.frame(predict(mod,newdata=data.frame(subsample=newdata),interval="predict")) # Getting predictions and coefficients
      
      polygon(x=c(newdata,rev(newdata)),y=c(pre$Lower,rev(pre$Upper)),border=rgb(t(col2rgb(cols[curve])),alpha=100,maxColorValue=255),col=rgb(t(col2rgb(cols[curve])),alpha=100,maxColorValue=255),xpd=TRUE)
      points(newdata,pre$Prediction,type="l",lty=3,lwd=2,col=cols[curve],xpd=TRUE)
    }
  }

}

# Printing rarefaction curves with Shannon indexes on legend ####
cols<-brewer.pal(8,"Spectral")
max<-max(unlist(Subsample))
jpeg("GRAPHICS/Transbaikal.Diversity.Rarefaction.jpg",height=2000,width=2000,res=300,quality=100)
par(lwd=2,mar=c(5,4,4,0.5))
plot(NA,NA,las=1,xlab=list("Sample size",cex=1.2,font=2),ylab=list("Species",cex=1.1,font=2),xlim=c(0,12000),ylim=c(0,1200),axes=FALSE)
axis(2,seq(0,1000,by=200),labels=TRUE,las=1,tick=TRUE,lwd.ticks=2,cex.axis=1.1,font.axis=2,lwd=2)
axis(1,seq(0,30000,by=5000),labels=TRUE,las=1,tick=TRUE,lwd.ticks=2,cex.axis=1.1,font.axis=2,lwd=2)
box(lwd=2)
for(c in 1:length(Curves)){
  points(Subsample[[c]],Curves[[c]],type="l",lwd=4,xpd=TRUE,col="black")        
  points(Subsample[[c]],Curves[[c]],type="l",lwd=2,xpd=TRUE,col=cols[curve])
  
  # For prediction of additional data on the curve ## See below how I got this
  if(Subsample[[c]][length(Subsample[[c]])] < max){
    dat<-data.frame(curve=Curves[[c]],subsample=Subsample[[c]]) # creating a data.frame
    mod<-drm(curve~subsample,data=dat,fct=W1.4())      # Calculating model         
    newdata<-seq(dat$subsample[length(dat$subsample)],max,by=50) # data frame for new data
    pre<-as.data.frame(predict(mod,newdata=data.frame(subsample=newdata),interval="predict")) # Getting predictions and coefficients
    
    polygon(x=c(newdata,rev(newdata)),y=c(pre$Lower,rev(pre$Upper)),border=colsP[c],col=colsP[c],xpd=TRUE)
    points(newdata,pre$Prediction,type="l",lty=3,lwd=2,col=cols[c],xpd=TRUE)
    
  }
}

arrows(0,seq(0,1000,by=200),35000,seq(0,1000,by=200),lwd=1,lty=3,col="grey80",code=0,xpd=TRUE)

#Legend
yl<-seq(300,0,length.out=ncol(TB3))
xl<-30000
rect(xl-1000,min(yl)-3,xl+2200,max(yl)+8,border="black",col="white")
for(l in 1:ncol(TB3)){
  text(xl,yl[l],labels=colnames(TB3)[l],cex=1,font=2,adj=c(1,0.5))
  arrows(xl+100,yl[l],xl+900,yl[l],code=0,lwd=5,col="black",lty=1)
  arrows(xl+100,yl[l],xl+900,yl[l],code=0,lwd=3,col=cols[l],lty=1)
  text(xl+1100,yl[l],labels=ShannonI[l],cex=1,font=2,adj=c(0,0.5),xpd=TRUE)
}
text(xl+1400,max(yl)+2.5,labels="Shannon\nindex",cex=0.8,font=2,adj=c(0.5,0),xpd=TRUE)

dev.off()



## Predicting further points to complete the curve using exponential model####

    # Making data frame with curve data
    dat<-data.frame(curve=unlist(Curves[1]),subsample=unlist(Subsample[1]))

    # Best model fitting to rarefaction curve seem Yield-loss curve,
    # which is a parametrization of Michaellis-Menten model.
    model<-drm(curve~subsample,data=dat,fct=W1.4()) # Michaellis-Menten would be MM.2 or MM.3
    
    # Predicting new data
    newdata<-seq(800,10000,by=50)
    pre<-predict(model,newdata=data.frame(subsample=newdata),interval="predict")
    
    # plot data
    x11()
    plot(NA,NA,las=1,xlab=list("Sample size",cex=1.2,font=2),ylab=list("Number of OTUs",cex=1.1,font=2),xlim=c(0,12000),ylim=c(0,85),axes=FALSE)
    axis(2,seq(0,80,by=20),labels=TRUE,las=1,tick=TRUE,lwd.ticks=2,cex.axis=1.1,font.axis=2,lwd=2)
    axis(1,seq(0,10000,by=2000),labels=TRUE,las=1,tick=TRUE,lwd.ticks=2,cex.axis=1.1,font.axis=2,lwd=2)
    points(dat$subsample,dat$curve,type="p",lwd=2,xpd=TRUE,col="red3") # Observed data
    points(dat$subsample,predict(model),lwd=1,type="l",col="red3")
    points(newdata,pre[,1],type="l",lwd=1,lty=3,col="red3")
    
    
    
    
