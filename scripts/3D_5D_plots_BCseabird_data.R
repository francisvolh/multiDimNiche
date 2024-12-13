
#be sure to use the data source 3d from niche rover run
#saveRDS(meanmu, "data/meanmu_for_3Dplot.RDS")
#saveRDS(output4, "data/output4_for_3Dplot.RDS")
MM<-read.csv("C:/Users/francis van oordt/OneDrive - McGill University/Documents/McGill/00Res Prop v2/Chap 0 BC isotopes/MM.csv") # read object with original filtered dataframe to run PCA on

nspp <- length(unique(MM$Species)) # number of species

output4 <- readRDS("data/output4_for_3Dplot.RDS")# for 100000 iterations 3D 5 sp from nicherover
meanmu <- readRDS("data/meanmu_for_3Dplot.RDS")# for 100000 iterations 3D 5 spp from nicherover
nGroups<-nspp
sigma<-output4 #out$mean$cov
mu<-meanmu #out$mean$mu


# Call and construct plot
#for 10 groups: c("orange", "red", "blue", "purple","purple4", "plum2", "black", "green", "palegreen", "green4")
#for 6 groups:c("orange", "red", "blue", "purple","black", "green")
#for 5 groups:c("orange", "red", "purple","black", "green")
#for PECO decades: colorit<-c("black","gray29", "gray52", "gray85")
#for DCCO: colorit<-c("red3","red", "orangered", "indianred1")
#for LESP: c("purple4","purple", "plum2")
colorit<-c("orange", "red", "purple","black", "green")
{rgl::open3d()
rgl::par3d(cex=1.75)
for (k in nGroups:1){
  nam <- paste("A", k, sep = "") 
  assign(nam, rgl::ellipse3d(qmesh=TRUE, sigma[,,k], centre=meanmu[k,],level=0.95, subdivide = 5, trans=diag(4)))
  rgl::wire3d(get(paste("A", k, sep="")), col=colorit[k])
  nam2 <- paste("cords", k, sep = "")
  assign(nam2,t(get(paste("A", k, sep=""))$vb)[,1:3])
}
rgl::axes3d(lwd=2, col="black")
rgl::grid3d("x-",lwd=2,col="black")
rgl::grid3d("y-",lwd=2,col="black")
rgl::grid3d("z+",lwd=2,col="black")
#rgl::axes3d(c('x--','x++','y--','y++','z--','z++'),expand=2)
rgl::par3d(windowRect = c(20, 30, 1000, 800))
M <- rgl::par3d("userMatrix") 

rgl::view3d(userMatrix = rgl::rotate3d(M, 3.6, 3.8, -0.8,  0.1))

rgl::mtext3d(colnames(sigma[,,1])[3], "z-+", line = 7)
rgl::mtext3d(colnames(sigma[,,1])[2], "y-", line = 7)
rgl::mtext3d(colnames(sigma[,,1])[1], "x-", line = 7)

#rgl::title3d(,,colnames(sigma[,,1])[1],colnames(sigma[,,1])[2],colnames(sigma[,,1])[3])

}
# Plot 3d and rotate two different angles

rgl::snapshot3d( "plots/test3D1stv2.png", fmt = "png", top = TRUE , webshot = TRUE) #needs the snapshot 3d package installed


{rgl::open3d()
  rgl::par3d(cex=1.75)
  
  for (k in nGroups:1){
    nam <- paste("A", k, sep = "") 
    assign(nam, rgl::ellipse3d(qmesh=TRUE, sigma[,,k], centre=meanmu[k,],level=0.95, subdivide = 5, trans=diag(4)))
    rgl::wire3d(get(paste("A", k, sep="")), col=colorit[k])
    nam2 <- paste("cords", k, sep = "")
    assign(nam2,t(get(paste("A", k, sep=""))$vb)[,1:3])
  }
  rgl::par3d(windowRect = c(20, 30, 1200, 800))
  
  rgl::axes3d(lwd=2, col="black")
  rgl::grid3d("x+",lwd=2,col="black")
  rgl::grid3d("y+",lwd=2,col="black")
  rgl::grid3d("z+",lwd=2,col="black")
  M <- rgl::par3d("userMatrix") 
  rgl::view3d(userMatrix = rgl::rotate3d(M, 3,4, 11,  0))
  
  #rgl::axes3d(c('x--','x++','y--','y++','z--','z++'),expand=2)
  
  rgl::mtext3d(colnames(sigma[,,1])[3], "z-+", line = 6)
  rgl::mtext3d(colnames(sigma[,,1])[2], "y-+", line = 5)
  rgl::mtext3d(colnames(sigma[,,1])[1], "x-+", line = 4)
  #rgl::title3d(,,colnames(sigma[,,1])[1],colnames(sigma[,,1])[2],colnames(sigma[,,1])[3])
}


rgl::snapshot3d( "plots/test3D2ndv2.png", fmt = "png", top = TRUE , webshot = TRUE) #needs the snapshot 3d package installed



colorit<-c( "#E69F00", "#D55E00" ,  "#CC79A7","black",  "#009E73")
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
#carbon<-expression(paste(delta^"13"*"C (\u2030)"))
#nitrogen<-expression(paste(delta^"15"*"N (\u2030)"))
carbon<- "PC1Css"
nitrogen<-"PC1Nss"

sulfur<-expression(paste(delta^"34"*"S (\u2030)"))
#rgb<-matrix(c(1,0,0,0,1,0,0,0,1),ncol=nGroups)

matrix.col<-matrix(rep(c(1,rep(0,nGroups)),nGroups),ncol=nGroups,byrow = T)
matrix.col<-matrix.col[1:nGroups,]
rgb<-matrix.col


#function for convex hulls from Rossman

hull<-function (x, y) 
{
  chI <- chull(x, y)
  chI <- c(chI, chI[1])
  hullx <- x[chI]
  hully <- y[chI]
  out<- list()
  out$x<-hullx
  out$y<-hully
  out
}


#Carbon Nitrogen
par(las=1)
par(mgp=c(0, 0.7, 0))
plot(NA,NA,ylim=c(nmin,nmax),xlim=c(cmin,cmax),xlab="",ylab="", bty="n", 
     xaxt="n",yaxt="n",xaxs="i",yaxs="i", main="", cex.lab=2, cex.main=2)
rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col = "lightgray", border=NA)
abline(v=seq(-100,100,by=1),col="white")
abline(h=seq(-100,100,by=1),col="white")
axis(side = 1, at = seq(-100,100,by=2), cex.axis=1.5,lwd=0, lwd.ticks=0)
axis(side = 2, at = seq(-100,100,by=2),cex.axis=1.5, lwd=0, lwd.ticks=0)
mtext(carbon, side=1, line=3.3, cex=2)
par(las=3)
mtext(nitrogen, side=2, line=1.8, cex=2)
for(i in 1:nGroups){
  polygon(hull(cords[,seq3[i]],cords[,(seq3[i]+1)]),col="lightgray",border=NA)
}
for(i in 1:nGroups){
  polygon(hull(cords[,seq3[i]],cords[,(seq3[i]+1)]),col=colorit[i],border=NA)
}

text(-5, -4,  substitute(paste(bold('Pelagic cormorant'))),
     cex=1, pos=3,col="black") 
text(-5, 4.5,  substitute(paste(bold('Double-crested cormorant'))),
     cex=1, pos=3,col="#D55E00") 
text(3.5, 3.8,  substitute(paste(bold("Leach's storm petrel"))),
     cex=1, pos=3,col="#752E58") 
text(1.75, 0,  substitute(paste(bold("Rhinoceros auklet"))),
     cex=1, pos=3,col="#024E37") 
text(2, -7,  substitute(paste(bold("Ancient murrelet"))),
     cex=1, pos=3,col="#714C02") 

dev.copy(png, filename = "plots/plot2dims.png", width = 2200, height = 1500, 
         res = 300)
dev.off()


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
mtext(carbon
  , side=1, line=3.3, cex=2)
par(las=3)
mtext(sulfur, side=2, line=1.8, cex=2)

for(i in 1:nGroups){
  polygon(magick::hull(cords[,seq3[i]],cords[,(seq3[i]+2)]),col="lightgray",border=NA)
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
  polygon(MASS::hull(cords[,(seq3[i]+2)],cords[,(seq3[i]+1)]),col="lightgray",border=NA)
}
for(i in 1:nGroups){
  polygon(hull(cords[,(seq3[i]+2)],cords[,(seq3[i]+1)]),col=colorit[i],border=NA)
}


### PLOTLY 5d plot
#plotly plot 3D
MM<-read.csv("C:/Users/francis van oordt/OneDrive - McGill University/Documents/McGill/00Res Prop v2/Chap 0 BC isotopes/MM.csv")
set.seed(4321)
nsamples<-1000

ndim<-5

system.time({
  
          dataCNS<-MM[,c("Species","PC1Css","PC2Css",
                         "PC1Nss","PC2Nss",
                         "stdDelta.S")]#subset with 5 dimension for runs
          BCiso.par <- tapply(1:nrow(dataCNS), dataCNS$Species,
                              function(ii) nicheROVER::niw.post(nsamples = nsamples,
                                                    X = dataCNS[ii, 2:ncol(dataCNS)]))
})

ncol <- 1:nspp # seq of species numbers
nrows <- 1:nsamples # seq of sample numbers

output <- NULL
for (i in 1:nspp) {
  data <- BCiso.par[[i]]
  mu <- as.data.frame(data[1])
  output <- abind::abind(output, mu, along = 3)
  
}

output2 <- NULL
for (i in 1:nisotope) {
  output1.1 <- output[,i,]
  output2 <- abind::abind(output2, output1.1, along = 3)
}
#be sure to use 5d data source from NicheRover

spp <- sort(unique(dataCNS$Species)) # species names

nspp <- length(unique(dataCNS$Species)) # number of species

isotopes <- colnames(dataCNS)[2:ncol(dataCNS)] # isotope label names

nisotope <- length(isotopes) # numner of isotopes/dimensions used


matrix.data <- output2[sample(nrow(output2), 1000), ,]

data <-matrix.data 

dataallspp <- NULL

for (i in 1:5) {
  spone<-spp[i]
  dataspp <- data.frame(rep(spone, length(data[,1,i])), data[,i,])
  
  dataallspp <- rbind(dataallspp, dataspp)
}

colnames(dataallspp) <- c("spp",c(isotopes))


dataallspp$PC2Nss2<-dataallspp$PC2Nss+10

fig2<-  plotly::plot_ly(dataallspp,
                      x = ~PC1Css, 
                      y = ~PC1Nss, 
                      z = ~stdDelta.S, 
                      color = ~PC2Css, size = ~PC2Nss2, sizes = c(5,450))
  fig2|> plotly::layout(
    scene = list(
                 camera = list(eye = list(x = cos(0.1)*1,
                                          y = sin(0.3)*1,
                                          z=4), # multiplying by 2 as zoom
                               center = list(x = 0,
                                             y = 0,
                                             z = 0
                               )
                 )
                 ))
  
  fig2|> plotly::layout(
    scene = list(
      camera = list(eye = list(x = cos(0.1)*2,
                               y = sin(0.3)*2,
                               z = -sin(2)), # multiplying by 2 as zoom
                    center = list(x = 0,
                                  y = 0,
                                  z = 0
                    )
      )
    ))
  
  

##external saving not really working or functional

Sys.setenv("PATH" = paste(Sys.getenv("PATH"), "C:/Users/francis van oordt/AppData/Local/Programs/orca", sep = .Platform$path.sep))

plotly::orca(fig2, "test15dim.png", width = 600, height =  700)


# found this, but not applicable, not sure where to get a key or what is this function calling
plotly_IMAGE(graph,
             width = 1200,
             height = 1050,
             format = "png",
             username="xxx",
             key="xxxx",
             scale = 1,
             out_file = paste(outfile,"png", sep="."))
cam.zoom = 2
plotly::ver.angle = 0

#################
# species grouping by colour

fig<-plotly::plot_ly(dataallspp, 
                     x = ~PC1Css, 
                     y = ~PC1Nss, 
                     z = ~stdDelta.S, 
                     color = ~Species, size = 1,
                     colors = c("red", "blue", "purple","black", "green"))|>
  plotly::layout( #not reay changing the size
    xaxis = list(tickfont = list(size = 15)), 
    yaxis = list(tickfont = list(size = 5)))
fig


p1 <- cowplot::ggdraw() + cowplot::draw_image("plots/plot2dims.png", scale = 0.9)

p2 <- cowplot::ggdraw() + cowplot::draw_image("plots/test3D1stv2.png", clip= "on", scale = 1.1)
p3 <- cowplot::ggdraw() + cowplot::draw_image("plots/test3D2ndv2.png", clip= "on",  scale = 1.1)# need to crop the image ahead
settwo<-cowplot::plot_grid(p2, p3, labels = c("B"),  ncol = 2)
p4 <- cowplot::ggdraw() + cowplot::draw_image("plots/5dim_lastscreen1stang.jpg")
p5 <- cowplot::ggdraw() + cowplot::draw_image("plots/newplot2ang2.png")
setthree<-cowplot::plot_grid(p4, p5, labels = c("C"))

full.set <-cowplot::plot_grid(p1, settwo , setthree, nrow = 3, labels = c("A"))

ggplot2::ggsave(full.set,  filename = "C:/Users/francis van oordt/OneDrive - McGill University/Documents/McGill/00Res Prop v2/Chap 0 BC isotopes/multiDimNiche/plots/full.set.png",  
                units = 'in', width = 9, height = 10, bg = "white")

