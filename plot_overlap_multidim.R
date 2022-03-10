
MM<-read.csv("MM.csv") #after running PCA and cleaning up columns in datasheet

nsamples<-100

ndim<-2
ndim<-3
ndim<-5

#RUN nicherover for corresponding dimension from MM dataframe

    if (ndim == 2) {
    dataCNS<-MM[,c("Species","PC1Css","PC1Nss")]
    BCiso.par <- tapply(1:nrow(dataCNS), dataCNS$Species,
                        function(ii) niw.post(nsamples = nsamples,
                                              X = dataCNS[ii, 2:ncol(dataCNS)])) 
  }else{
    if (ndim == 3){
      dataCNS<-MM[,c("Species","PC1Css","PC1Nss","stdDelta.S")]
      BCiso.par <- tapply(1:nrow(dataCNS), dataCNS$Species,
                          function(ii) niw.post(nsamples = nsamples,
                                                X = dataCNS[ii, 2:ncol(dataCNS)]))
    }else{
      if (ndim == 5) {
        dataCNS<-MM[,c("Species","PC1Css","PC2Css",
                       "PC1Nss","PC2Nss",
                       "stdDelta.S")]#subset with 5 dimension for runs
        BCiso.par <- tapply(1:nrow(dataCNS), dataCNS$Species,
                            function(ii) niw.post(nsamples = nsamples,
                                                  X = dataCNS[ii, 2:ncol(dataCNS)]))
      }
    }
  } 
#Overlap metrics
over.stat <- overlap(BCiso.par, nreps = nsamples, nprob = 1000, alpha = c(0.95, 
                                                                          0.99))

# The mean overlap metrics calculated across iteratations for both niche
# region sizes (alpha = .95 and alpha = .99) can be calculated and displayed
# in an array.
over.mean <- apply(over.stat, c(1:2, 4), mean) * 100
round(over.mean, 2)


overmean<-as.matrix(over.mean[,,1])
overmean<-as.numeric(overmean)
Dim2<-overmean

Dim3<-overmean

Dim5<-overmean

all.dim.over<-as.data.frame(rbind(overmean2,overmean3,overmean5), na.m=TRUE)

colname <- sort(unique(MM$Species))
rowname <- sort(unique(MM$Species))

labels.1 <- NULL

for (i in 1:length(rowname)) {
  for (j in 1:length(colname)) {
    label <- paste0(colname[j], "-", rowname[i])
    
    labels.1 <- c(labels.1, label) 
    
  }
}
colnames(all.dim.over)<-labels.1

all.dim.over<-all.dim.over[, which(colSums(all.dim.over) > 0)]


par(mar = c(6,3,1,1))
plot(as.numeric(all.dim.over[1,]), pch = 19, 
     col = "red", ylim= c(0,100), axes=FALSE,
     xlab = NA, ylab=NA)
par(new = TRUE)
plot(as.numeric(all.dim.over[2,]), pch = 19, 
     col = "blue", ylim = c(0,100), axes=FALSE,
     xlab = NA, ylab=NA)
par(new = TRUE)
plot(as.numeric(all.dim.over[3,]), pch = 19, 
     col = "dark grey", ylim = c(0,100), axes=FALSE,
     xlab = NA, ylab=NA)
axis(side = 1, at = 1:ncol(all.dim.over), 
     labels = colnames(all.dim.over), las=2, cex.axis = 0.7)
axis(side = 2, at = seq(0,100,10),cex.axis = 0.7)
legend("top", legend = c("2D", "3D", "5D"), 
       pch = 19, col = c("red", "blue", "dark grey"), bty = "n", horiz = TRUE)
box()
dev.copy(png, filename = "overlapplot.png", width = 3800, height = 2000, 
         res = 300)
dev.off()
