#install.packages("fmsb")
#library(fmsb)

#load updated function from fmsb that changes the size of the points

radarchart <- function(df, axistype=0, seg=4, pty=16, pcol=1:8, plty=1:6, plwd=1,
                       pdensity=NULL, pangle=45, pfcol=NA, cglty=3, cglwd=1,
                       cglcol="navy", axislabcol="blue", title="", maxmin=TRUE,
                       na.itp=TRUE, centerzero=FALSE, vlabels=NULL, vlcex=NULL,
                       caxislabels=NULL, calcex=NULL,
                       paxislabels=NULL, palcex=NULL, ...) {
  if (!is.data.frame(df)) { cat("The data must be given as dataframe.\n"); return() }
  if ((n <- length(df))<3) { cat("The number of variables must be 3 or more.\n"); return() }
  if (maxmin==FALSE) { # when the dataframe does not include max and min as the top 2 rows.
    dfmax <- apply(df, 2, max)
    dfmin <- apply(df, 2, min)
    df <- rbind(dfmax, dfmin, df)
  }
  plot(c(-1.2, 1.2), c(-1.2, 1.2), type="n", frame.plot=FALSE, axes=FALSE, 
       xlab="", ylab="", main=title, asp=1, ...) # define x-y coordinates without any plot
  theta <- seq(90, 450, length=n+1)*pi/180
  theta <- theta[1:n]
  xx <- cos(theta)
  yy <- sin(theta)
  CGap <- ifelse(centerzero, 0, 1)
  for (i in 0:seg) { # complementary guide lines, dotted navy line by default
    polygon(xx*(i+CGap)/(seg+CGap), yy*(i+CGap)/(seg+CGap), lty=cglty, lwd=cglwd, border=cglcol)
    if (axistype==1|axistype==3) CAXISLABELS <- paste(i/seg*100,"(%)")
    if (axistype==4|axistype==5) CAXISLABELS <- sprintf("%3.2f",i/seg)
    if (!is.null(caxislabels)&(i<length(caxislabels))) CAXISLABELS <- caxislabels[i+1]
    if (axistype==1|axistype==3|axistype==4|axistype==5) {
      if (is.null(calcex)) text(-0.05, (i+CGap)/(seg+CGap), CAXISLABELS, col=axislabcol) else
        text(-0.05, (i+CGap)/(seg+CGap), CAXISLABELS, col=axislabcol, cex=calcex)
    }
  }
  if (centerzero) {
    arrows(0, 0, xx*1, yy*1, lwd=cglwd, lty=cglty, length=0, col=cglcol)
  }
  else {
    arrows(xx/(seg+CGap), yy/(seg+CGap), xx*1, yy*1, lwd=cglwd, lty=cglty, length=0, col=cglcol)
  }
  PAXISLABELS <- df[1,1:n]
  if (!is.null(paxislabels)) PAXISLABELS <- paxislabels
  if (axistype==2|axistype==3|axistype==5) {
    if (is.null(palcex)) text(xx[1:n], yy[1:n], PAXISLABELS, col=axislabcol) else
      text(xx[1:n], yy[1:n], PAXISLABELS, col=axislabcol, cex=palcex)
  }
  VLABELS <- colnames(df)
  if (!is.null(vlabels)) VLABELS <- vlabels
  if (is.null(vlcex)) text(xx*1.2, yy*1.2, VLABELS) else
    text(xx*1.2, yy*1.2, VLABELS, cex=vlcex)
  series <- length(df[[1]])
  SX <- series-2
  if (length(pty) < SX) { ptys <- rep(pty, SX) } else { ptys <- pty }
  if (length(pcol) < SX) { pcols <- rep(pcol, SX) } else { pcols <- pcol }
  if (length(plty) < SX) { pltys <- rep(plty, SX) } else { pltys <- plty }
  if (length(plwd) < SX) { plwds <- rep(plwd, SX) } else { plwds <- plwd }
  if (length(pdensity) < SX) { pdensities <- rep(pdensity, SX) } else { pdensities <- pdensity }
  if (length(pangle) < SX) { pangles <- rep(pangle, SX)} else { pangles <- pangle }
  if (length(pfcol) < SX) { pfcols <- rep(pfcol, SX) } else { pfcols <- pfcol }
  for (i in 3:series) {
    xxs <- xx
    yys <- yy
    scale <- CGap/(seg+CGap)+(df[i,]-df[2,])/(df[1,]-df[2,])*seg/(seg+CGap)
    if (sum(!is.na(df[i,]))<3) { cat(sprintf("[DATA NOT ENOUGH] at %d\n%g\n",i,df[i,])) # for too many NA's (1.2.2012)
    } else {
      for (j in 1:n) {
        if (is.na(df[i, j])) { # how to treat NA
          if (na.itp) { # treat NA using interpolation
            left <- ifelse(j>1, j-1, n)
            while (is.na(df[i, left])) {
              left <- ifelse(left>1, left-1, n)
            }
            right <- ifelse(j<n, j+1, 1)
            while (is.na(df[i, right])) {
              right <- ifelse(right<n, right+1, 1)
            }
            xxleft <- xx[left]*CGap/(seg+CGap)+xx[left]*(df[i,left]-df[2,left])/(df[1,left]-df[2,left])*seg/(seg+CGap)
            yyleft <- yy[left]*CGap/(seg+CGap)+yy[left]*(df[i,left]-df[2,left])/(df[1,left]-df[2,left])*seg/(seg+CGap)
            xxright <- xx[right]*CGap/(seg+CGap)+xx[right]*(df[i,right]-df[2,right])/(df[1,right]-df[2,right])*seg/(seg+CGap)
            yyright <- yy[right]*CGap/(seg+CGap)+yy[right]*(df[i,right]-df[2,right])/(df[1,right]-df[2,right])*seg/(seg+CGap)
            if (xxleft > xxright) {
              xxtmp <- xxleft; yytmp <- yyleft;
              xxleft <- xxright; yyleft <- yyright;
              xxright <- xxtmp; yyright <- yytmp;
            }
            xxs[j] <- xx[j]*(yyleft*xxright-yyright*xxleft)/(yy[j]*(xxright-xxleft)-xx[j]*(yyright-yyleft))
            yys[j] <- (yy[j]/xx[j])*xxs[j]
          } else { # treat NA as zero (origin)
            xxs[j] <- 0
            yys[j] <- 0
          }
        }
        else {
          xxs[j] <- xx[j]*CGap/(seg+CGap)+xx[j]*(df[i, j]-df[2, j])/(df[1, j]-df[2, j])*seg/(seg+CGap)
          yys[j] <- yy[j]*CGap/(seg+CGap)+yy[j]*(df[i, j]-df[2, j])/(df[1, j]-df[2, j])*seg/(seg+CGap)
        }
      }
      if (is.null(pdensities)) {
        polygon(xxs, yys, lty=pltys[i-2], lwd=plwds[i-2], border=pcols[i-2], col=pfcols[i-2])
      } else {
        polygon(xxs, yys, lty=pltys[i-2], lwd=plwds[i-2], border=pcols[i-2], 
                density=pdensities[i-2], angle=pangles[i-2], col=pfcols[i-2])
      }
      points(xx*scale, yy*scale, cex=2, pch=ptys[i-2], col=pcols[i-2])
    }
  }
}

par(mfrow = c(1,1), mar = c(2,2,2,2))

centroids<-read.csv("C:/Users/francis van oordt/OneDrive - McGill University/Documents/McGill/00Res Prop v2/Chap 0 BC isotopes/centroids.nicheRoverRes5spp.csv", header = TRUE)
centroids
centroids <-rbind(rep(5,6) , rep(-5,6) , centroids) #adding lines to center points correctly

rgb.valRED <- col2rgb("red")
rgb.valOrange <- col2rgb("orange")
rgb.valPURP <- col2rgb("purple")
rgb.valBLK <- col2rgb("black")
rgb.valGRE <- col2rgb("green")




rgb.valRED <- rgb(rgb.valRED[1], rgb.valRED[2], rgb.valRED[3], max = 255, alpha = 95, names = "red")
rgb.valOrange <-rgb(rgb.valOrange[1], rgb.valOrange[2], rgb.valOrange[3], max = 255, alpha = 95, names = "orange")
rgb.valPURP <- rgb(rgb.valPURP[1], rgb.valPURP[2], rgb.valPURP[3], max = 255, alpha = 95, names = "purple")
rgb.valBLK <-rgb(0, 0, 0, max = 255, alpha = 95, names = "black")
rgb.valGRE <- rgb(rgb.valGRE[1], rgb.valGRE[2], rgb.valGRE[3], max = 255, alpha = 95, names = "green")

radarchart(centroids[,2:6], axistype=1, 
            
            #custom polygon
            pcol=c(rgb.valOrange, rgb.valRED,rgb.valPURP, rgb.valBLK, rgb.valGRE), plwd=1, plty=1, 
            
            #custom the grid
            cglcol="grey", cglty=1, axislabcol="grey", caxislabels=seq(-5,5,2.5), cglwd=1.1,
            
            #custom labels
            vlcex=0.8 
)
legend(x=0.85, y=1.2, legend = c("ANMU", "DCCO","LESP", "PECO", "RHAU"), bty = "n", pch=20 , 
       col=c(rgb.valOrange, rgb.valRED,rgb.valPURP, rgb.valBLK, rgb.valGRE) , 
       text.col = "black", cex=0.9, pt.cex=1.6)

dev.copy(png, filename = "plots/spider5spp.png", width = 2000, height = 1700, 
         res = 300)
dev.off()



################## DCCO PECO Plots

rgb.valRED <- col2rgb("red")
rgb.valOrange <- col2rgb("orange")
rgb.valRED <- rgb(rgb.valRED[1], rgb.valRED[2], rgb.valRED[3], max = 255, alpha = 95, names = "red")
rgb.valOrange <-rgb(rgb.valOrange[1], rgb.valOrange[2], rgb.valOrange[3], max = 255, alpha = 95, names = "orange")


par(mfrow = c(1,2), mar = c(0.1, 0.1, 1, 0.1))
centroids<-read.csv( "C:/Users/francis van oordt/OneDrive - McGill University/Documents/McGill/00Res Prop v2/Chap 0 BC isotopes/centroidsDCCO.csv", header = TRUE)
centroids
centroids <-rbind(rep(3,6) , rep(-4,6) , centroids) #adding lines to center points correctly

radarchart(centroids[,2:6], axistype=1, 
           
           #custom polygon
           pcol=c(rgb.valOrange, rgb.valRED) , plwd=1, plty=3, 
           
           #custom the grid
           cglcol="grey", cglty=1, axislabcol="grey",caxislabels=seq(-4,3,1.75), cglwd=1.1,
           
           #custom labels
           vlcex=0.8 ,
           title = "DCCO"
)
legend(x=0.85, y=1, legend = c("1970-1989", "1990-2006"), bty = "n", pch=20 , 
       col=c(rgb.valOrange, rgb.valRED) , 
       text.col = "black", cex=0.9, pt.cex=1.6)




centroids<-read.csv( "C:/Users/francis van oordt/OneDrive - McGill University/Documents/McGill/00Res Prop v2/Chap 0 BC isotopes/centroidsPECO.csv", header = TRUE)
centroids
centroids <-rbind(rep(3,6) , rep(-4,6) , centroids) #adding lines to center points correctly

radarchart(centroids[,2:6], axistype=1, 
           
           #custom polygon
           pcol=c(rgb.valOrange, rgb.valRED) , plwd=1, plty=3, 
           
           #custom the grid
           cglcol="grey", cglty=1, axislabcol="grey", caxislabels=seq(-4,3,1.75), cglwd=1.1,
           
           #custom labels
           vlcex=0.8 ,
           title = "PECO"
)
legend(x=0.85, y=1, legend = c("1970-1989", "1990-2006"), bty = "n", pch=20 , 
       col=c(rgb.valOrange, rgb.valRED) , 
       text.col = "black", cex=0.9, pt.cex=1.6)

dev.copy(png, filename = "plots/spiderDCCOPECOv2.png", width = 4000, height = 1900, 
         res = 300)
dev.off()
  