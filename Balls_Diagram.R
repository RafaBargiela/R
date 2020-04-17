
DrawBalls<-function(M,pal,balls=c(0.001,0.01,0.1,1,10,100),rr=seq(0.8,2.25,length.out=length(balls))){
  # M : Matrix of values to plot
  # pal : color palette (vector or matrix). Should be same length than balls and rr vectors, from light to dark.
  # balls : Vector with values determining the boundaries for each ball to which the values are adjusted
  # rr : Vector with corresponding radius values for each ball level. 
  # rr and balls must be same length
  pal<-as.matrix(pal)
  
  for(r in 1:nrow(M)){
    palette<-1
    for(c in 1:ncol(M)){
      n<-as.numeric(M[r,c])
      if(n==0){
        # print (n)           
      }else{
        if(n<min(balls)){
          points(c,(nrow(M)+1)-r,pch=4,cex=rr[1],lwd=3.5,xpd=TRUE,col="white")  
          points(c,(nrow(M)+1)-r,pch=4,cex=rr[1]-(rr[1]/4),lwd=3,xpd=TRUE,col=pal[1,palette])
        }else{
          v<-sum(balls<=n) # Fine closest minimum value
          # v<-balls[which.min(abs(balls-n))]  # Fine the closest number
          points(c,(nrow(M)+1)-r,pch=19,cex=rr[v]+(rr[v]/10),xpd=TRUE,col="white")           
          for(rad in v:1){
            points(c,(nrow(M)+1)-r,pch=19,cex=rr[rad],xpd=TRUE,col=pal[rad,palette])
          }
        }
      }
      palette<-palette+1  # Changing the column on the matrix color palette
      if(palette>ncol(pal)){ # if M has more columns than PAL, the palette is restarted
        palette<-1
      }
      
    }
  }
}