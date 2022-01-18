# load needed libraries
library(MASS)
library(clusterGeneration)
library(repmis)
library(rgl)
library(shape)
library(magick)
library(jagsUI) # to run jags
library(plotly)

#be sure to use the data source 3d from niche rover run

nGroups<-nspp
sigma<-output4#out$mean$cov
mu<-meanmu#out$mean$mu

# Call and construct plot
#for 10 groups: c("orange", "red", "blue", "purple","purple4", "plum2", "black", "green", "palegreen", "green4")
#for 6 groups:c("orange", "red", "blue", "purple","black", "green")
#for 5 groups:c("orange", "red", "purple","black", "green")
#for PECO decades: colorit<-c("black","gray29", "gray52", "gray85")
#for DCCO: colorit<-c("red3","red", "orangered", "indianred1")
#for LESP: c("purple4","purple", "plum2")
colorit<-c("orange", "red", "purple","black", "green")
open3d()
for (k in nGroups:1){
  nam <- paste("A", k, sep = "") 
  assign(nam, ellipse3d(qmesh=TRUE,sigma[,,k],centre=meanmu[k,],level=0.95, subdivide = 5, trans=diag(4)))
  wire3d(get(paste("A", k, sep="")), col=colorit[k])
  nam2 <- paste("cords", k, sep = "")
  assign(nam2,t(get(paste("A", k, sep=""))$vb)[,1:3])
}
axes3d(lwd=2, col="black")
grid3d("x+",lwd=2,col="black")
grid3d("y+",lwd=2,col="black")
grid3d("z+",lwd=2,col="black")
axes3d(c('x--','x++','y--','y++','z--','z++'),expand=2)
title3d(,,'Carbon','Nitrogen','Sulfur')

# Save a movie rotation
M <- par3d("userMatrix")
if (!rgl.useNULL())
  play3d( par3dinterp(time=(0:2)*0.75,userMatrix=list(M,
                                                      rotate3d(M, pi/2, 1, 0, 0),
                                                      rotate3d(M, pi/2, 0, 1, 0) ) ), 
          duration=10 )
#movie3d( spin3d(axis=c(0,1,0),rpm=3), duration=10, dir=getwd() )


#plotly plot 3D

fig <- plot_ly(C, x = ~isotope1, y = ~isotope2, z = ~isotope3, color = ~group)
fig <- fig %>% add_markers()
fig <- fig %>% layout(scene = list(xaxis = list(title = 'Carbon'),
                                   yaxis = list(title = 'Nitrogen'),
                                  zaxis = list(title = 'Sulfur')))
fig

# 2D Graphs (must run 3d graph first)
cords<-data.frame(matrix(NA, nrow=nrow(cords1),ncol=nGroups*3))
seq3<-seq(1,nGroups*3,by=3)
for (i  in 1:nGroups){
  cords[,seq3[i]:(seq3[i]+2)]<-get(paste("cords",i,sep=""))
}

cmin<-floor(min(cords[,seq3]))
cmax<-ceiling(max(cords[,seq3]))
nmin<-floor(min(cords[,(seq3+1)]))
nmax<-ceiling(max(cords[,(seq3+1)]))
smin<-floor(min(cords[,(seq3+2)]))-1
smax<-ceiling(max(cords[,(seq3+2)]))
carbon<-expression(paste(delta^"13"*"C (\211)"))
nitrogen<-expression(paste(delta^"15"*"N (\211)"))
sulfur<-expression(paste(delta^"34"*"S (\211)"))
#rgb<-matrix(c(1,0,0,0,1,0,0,0,1),ncol=nGroups)

matrix.col<-matrix(rep(c(1,rep(0,nGroups)),nGroups),ncol=nGroups,byrow = T)
matrix.col<-matrix.col[1:nGroups,]
rgb<-matrix.col

#Carbon Nitrogen
par(las=1)
par(mgp=c(0, 0.7, 0))
plot(NA,NA,ylim=c(nmin,nmax),xlim=c(cmin,cmax),xlab="",ylab="", bty="n", 
     xaxt="n",yaxt="n",xaxs="i",yaxs="i", main="", cex.lab=2, cex.main=2)
rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col = "lightgray", border=NA)
abline(v=seq(-100,100,by=1),col="white")
abline(h=seq(-100,100,by=1),col="white")
axis(side = 1, at = seq(-100,100,by=1), cex.axis=1.5,lwd=0, lwd.ticks=0)
axis(side = 2, at = seq(-100,100,by=1),cex.axis=1.5, lwd=0, lwd.ticks=0)
mtext(carbon, side=1, line=3.3, cex=2)
par(las=3)
mtext(nitrogen, side=2, line=1.8, cex=2)
for(i in 1:nGroups){
  polygon(hull(cords[,seq3[i]],cords[,(seq3[i]+1)]),col="lightgray",border=NA)
}
for(i in 1:nGroups){
  polygon(hull(cords[,seq3[i]],cords[,(seq3[i]+1)]),col=colorit[i],border=NA)
}

#Carbon Sulfur
par(las=1)
par(mgp=c(0, 0.7, 0))
plot(NA,NA,ylim=c(smin,smax),xlim=c(cmin,cmax),xlab="",ylab="", bty="n", 
     xaxt="n",yaxt="n",xaxs="i",yaxs="i", main="", cex.lab=2, cex.main=2)
rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col = "lightgray", border=NA)
abline(v=seq(-100,100,by=1),col="white")
abline(h=seq(-100,100,by=2),col="white")
axis(side = 1, at = seq(-100,100,by=1), cex.axis=1.5,lwd=0, lwd.ticks=0)
axis(side = 2, at = seq(-100,100,by=2),cex.axis=1.5, lwd=0, lwd.ticks=0)
mtext(carbon, side=1, line=3.3, cex=2)
par(las=3)
mtext(sulfur, side=2, line=1.8, cex=2)

for(i in 1:nGroups){
  polygon(hull(cords[,seq3[i]],cords[,(seq3[i]+2)]),col="lightgray",border=NA)
}
for(i in 1:nGroups){
  polygon(hull(cords[,seq3[i]],cords[,(seq3[i]+2)]),col=colorit[i],border=NA)
}

# Sulfur Nitrogen
par(las=1)
par(mgp=c(0, 0.7, 0))
plot(NA,NA,xlim=c(smin,smax),ylim=c(nmin,nmax),xlab="",ylab="", bty="n", 
     xaxt="n",yaxt="n",xaxs="i",yaxs="i", main="", cex.lab=2, cex.main=2)
rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col = "lightgray", border=NA)
abline(v=seq(-100,100,by=2),col="white")
abline(h=seq(-100,100,by=1),col="white")
axis(side = 1, at = seq(-100,100,by=2), cex.axis=1.5,lwd=0, lwd.ticks=0)
axis(side = 2, at = seq(-100,100,by=1),cex.axis=1.5, lwd=0, lwd.ticks=0)
mtext(sulfur, side=1, line=3.3, cex=2)
par(las=3)
mtext(nitrogen, side=2, line=1.8, cex=2)
for(i in 1:nGroups){
  polygon(hull(cords[,(seq3[i]+2)],cords[,(seq3[i]+1)]),col="lightgray",border=NA)
}
for(i in 1:nGroups){
  polygon(hull(cords[,(seq3[i]+2)],cords[,(seq3[i]+1)]),col=colorit[i],border=NA)
}
