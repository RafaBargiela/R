##########################################################################################
###		            Genome/Proteome ring graph                                   ###########
### Author: Dr. Rafael Bargiela Bargiela                                          ########
### Resume: An example of drawing data over whole proteome/Genome                 ########
###         circular representation                                               ######## 
##########################################################################################

library(aspace)
library(seqinr)
library(TeachingDemos)
library(png)
library(RColorBrewer)
library(stringr)

#### COLORS AND GRADIENTS ####
cols<-c("indianred","yellowgreen","royalblue1") # Colors reprenseting each sample
#### LOADING DATA ####
TB<-read.table("Prueba.txt",header=TRUE,row.names=1,check.names=FALSE,dec=".",sep="\t") ## CIRCULAR PROTEOME: 3 CIRCLES (I,A,E) USING REL. ABUN. (5 JUN 2019) ##
M<-as.matrix(TB[,c(1,3,5)]) # Matrix with values
Md<-as.matrix(TB[,c(2,4,6)]) # Matrix with SD


#### SETTING PREVIOUS VALUES BEFORE PRINTING THE DIAGRAM ####
lvExp<-c(0.001,0.01,0.1,1,10,100,1000) # Exponencial expression levels (subcircles)
AI<-98 # Angle where first KO is drawn
AF<-442# Angle where last KO is drawn
AW<-(360/nrow(M))/2 # Because we are going to make bars, is divided by 2
Samples<-c("I","A","E") # Samples names

#### PRINTING THE DIAGRAM ####
jpeg("GRAPHICS/Aging.Proteome.KO.jpg",height=4000,width=4000,res=600,quality=100)
anSeq<-seq(AI+AW,AF-AW,length.out=nrow(M)) # Angles where KO are distributed over the circumference
ri<-c(0.125,0.275,0.425) # Initial radius for base circles
rw<-0.0175 # Additional wide to add to the radius for each subcircle
par(mar=c(3,3,3,3))
plot(NA,NA,xlim=c(0,1),ylim=c(0,1),axes=FALSE,ylab="",xlab="")
lvExpLabs<-c("0.001","0.01","0.1","1","10","100","1000")
lvExpCex<-c(0.5,0.6,0.7)
yl<-1 # Values for final legend
xl<--0.1
for(s in 1:ncol(M)){
  circle(0.5,0.5,r=ri[s],border="black",col="transparent",xpd=TRUE) # Base circle
  for(sub in 1:length(lvExp)){  # Drawing just subcircles
    circle(0.5,0.5,r=ri[s]+rw*sub,border="black",col="transparent",xpd=TRUE,lty=3)  
    shadowtext(0.5+(ri[s]+rw*sub)*cos(as_radians(90)),0.5+(ri[s]+rw*sub)*sin(as_radians(90)),labels=lvExpLabs[sub],cex=lvExpCex[s],bg="white",col="black",xpd=TRUE,font=2)
  }
  
  for(k in 1:nrow(M)){
    n<-as.numeric(M[k,s]) # Expression value
    nSd<-as.numeric(Md[k,s]) # SD of the expression value
    if(n>0){
      rr<-sum(lvExp<=n) # Closest value to n in lvExp (position in the vector)
      rrr<-n*1/lvExp[rr+1] # Additional proportion to add to the radius to get specific value
      rrrSd<-(n+nSd)*1/lvExp[rr+1] # Same but for SD
      rf<-ri[s]+(rw*rr)+(rw*rrr) # End radius for each KO is the Initial radius (ri) + radius of corresponding subcircle + proporcional % inside each level
      rfSd<-ri[s]+(rw*rr)+(rw*rrrSd)
      # Calculating all coordinates for each specific bar
      x<-c(0.5+ri[s]*cos(as_radians(seq(anSeq[k]-AW,anSeq[k]+AW,length.out=50))),rev(0.5+(rf)*cos(as_radians(seq(anSeq[k]-AW,anSeq[k]+AW,length.out=50)))))
      y<-c(0.5+ri[s]*sin(as_radians(seq(anSeq[k]-AW,anSeq[k]+AW,length.out=50))),rev(0.5+(rf)*sin(as_radians(seq(anSeq[k]-AW,anSeq[k]+AW,length.out=50)))))
      polygon(x,y,border="black",lwd=0.1,col=cols[s],xpd=TRUE)
      # Drawing Error lines
      x1Sd<-0.5+rf*cos(as_radians(anSeq[k]))
      x2Sd<-0.5+rfSd*cos(as_radians(anSeq[k]))
      y1Sd<-0.5+rf*sin(as_radians(anSeq[k]))
      y2Sd<-0.5+rfSd*sin(as_radians(anSeq[k]))
      arrows(x1Sd,y1Sd,x2Sd,y2Sd,code=3,length=0.01,xpd=TRUE,lwd=1,angle=90)
    }
  }
  rect(xl-0.02,yl-0.02,xl+0.02,yl+0.02,border="black",col=cols[s],xpd=TRUE)
  text(xl+0.04,yl,labels=Samples[s],cex=1.2,font=2,xpd=TRUE)
  yl<-yl-0.05
}
text(0.5,0.5,labels="Relative\nabundance\n(ng/µg)",cex=1.2,font=2,xpd=TRUE)
# Printing outer ring
outR<-0.565
lines(x=0.5+outR*cos(as_radians(seq(AI-AW,AF+AW,length.out=nrow(M)))),y=0.5+outR*sin(as_radians(seq(AI-AW,AF+AW,length.out=nrow(M)))),xpd=TRUE,lwd=2)
for (r in 1:nrow(M)){
  if(r%%15==0 | grepl("K01667|K01696",rownames(M)[r])==TRUE){
    lines(x=c(0.5+outR*cos(as_radians(anSeq[r])),0.5+(outR+0.01)*cos(as_radians(anSeq[r]))),y=c(0.5+outR*sin(as_radians(anSeq[r])),0.5+(outR+0.01)*sin(as_radians(anSeq[r]))),lwd=2,xpd=TRUE)
    adj<-c(0,0.5)
    aL<-anSeq[r]
    if(anSeq[r]>90 & anSeq[r]<270){
      aL<-anSeq[r]-180
      adj<-c(1,0.5)
    }
    text(0.5+(outR+0.015)*cos(as_radians(anSeq[r])),0.5+(outR+0.015)*sin(as_radians(anSeq[r])),srt=aL,labels=rownames(M)[r],cex=0.6,font=2,xpd=TRUE,adj=adj,
         col=ifelse(grepl("K01667|K01696",rownames(M)[r])==TRUE,"red4","black"))
  } 
}
dev.off()


