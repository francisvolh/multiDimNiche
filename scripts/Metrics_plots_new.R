###niches
setwd("C:/Users/francis van oordt/OneDrive - McGill University/Documents/McGill/00Res Prop v2/Chap 0 BC isotopes/")
load("niches_1_100k.RData")
niches1d<-niches

load("niches_2_100k.RData")
niches2d<-niches
load("niches_3_100k.RData")
niches3d<-niches
load("niches_5_100k.RData")
niches5d<-niches

par(mfrow=c(4,1), mar = c(0, 0, 1, 0), oma = c(4,4,1,1))
boxplot(niches1d, axes = FALSE)
axis(side = 2, at = seq(0,c(range(niches1d)[2]+0.5),50))
mtext("A", side = 2, line = -1, at = trunc(c(range(niches1d)[2]+0.5)),  las= 2)#
box()

boxplot(niches2d, axes = FALSE)
axis(side = 2, at = seq(0,c(range(niches2d)[2]+0.5),50))
mtext("B", side = 2, line = -1, at = trunc(c(range(niches2d)[2]+0.5)),  las= 2)#
box()

boxplot(niches3d, axes = FALSE)
axis(side = 2, at = seq(0,c(range(niches3d)[2]+0.5),100))
mtext("C", side = 2, line = -1, at = trunc(c(range(niches3d)[2]+0.5)),  las= 2)#
box()


boxplot(niches5d, axes = FALSE)
axis(side = 2, at = seq(0,c(range(niches5d)[2]+0.5),400))
axis(side = 1, at = 1:5, labels = c("ANMU", "DCCO", "LSPE", "PECO", "RHAU"))
mtext("D", side = 2, line = -1, at = trunc(c(range(niches5d)[2]+0.5)),  las= 2)#
mtext("Niche Size", side = 1, line = 3, cex = 1, font = 1, outer = TRUE)
box()
dev.copy(jpeg, filename = "nicheSize5dims.jpg", width = 2000, height = 3000, 
         res = 300)
dev.off()
#################################################
#load CC outputs

#1d

load("output5_1_100k.RData")
ndim1cd <- as.numeric(output5)

load("output12_1_100k.RData")
ndim1 <- as.numeric(output12)


load("matrixrange_1_100k.RData")
range1d<-matrixrange

load("samplessd_1_100k.RData")
sdnnd1d<-samplessd




#2D
load("C:/Users/francis van Oordt/OneDrive - McGill University/Documents/McGill/00Res Prop v2/Chap 0 BC isotopes/nicheRover/CC outputs/new CC fixed/new run 02 10 Feb 2020/output5_2_100k.RData")
ndim2cd <- as.numeric(output5)

load("C:/Users/francis van Oordt/OneDrive - McGill University/Documents/McGill/00Res Prop v2/Chap 0 BC isotopes/nicheRover/CC outputs/new CC fixed/new run 02 10 Feb 2020/output12_2_100k.RData")
ndim2 <- as.numeric(output12)


load("C:/Users/francis van Oordt/OneDrive - McGill University/Documents/McGill/00Res Prop v2/Chap 0 BC isotopes/nicheRover/CC outputs/new CC fixed/run new metrics Apr 2021/matrixrange_2_100k.RData")
range2d<-matrixrange

load("C:/Users/francis van Oordt/OneDrive - McGill University/Documents/McGill/00Res Prop v2/Chap 0 BC isotopes/nicheRover/CC outputs/new CC fixed/run new metrics Apr 2021/samplessd_2_100k.RData")
sdnnd2d<-samplessd

#3D
load("C:/Users/francis van Oordt/OneDrive - McGill University/Documents/McGill/00Res Prop v2/Chap 0 BC isotopes/nicheRover/CC outputs/new CC fixed/new run 02 10 Feb 2020/output5_3_100k.RData")
ndim3cd <- as.numeric(output5)
load("C:/Users/francis van Oordt/OneDrive - McGill University/Documents/McGill/00Res Prop v2/Chap 0 BC isotopes/nicheRover/CC outputs/new CC fixed/new run 02 10 Feb 2020/output12_3_100k.RData")
ndim3 <- as.numeric(output12)

load("C:/Users/francis van Oordt/OneDrive - McGill University/Documents/McGill/00Res Prop v2/Chap 0 BC isotopes/nicheRover/CC outputs/new CC fixed/run new metrics Apr 2021/matrixrange_3_100k.RData")
range3d<-matrixrange
load("C:/Users/francis van Oordt/OneDrive - McGill University/Documents/McGill/00Res Prop v2/Chap 0 BC isotopes/nicheRover/CC outputs/new CC fixed/run new metrics Apr 2021/samplessd_3_100k.RData")
sdnnd3d<-samplessd

#5D
load("C:/Users/francis van Oordt/OneDrive - McGill University/Documents/McGill/00Res Prop v2/Chap 0 BC isotopes/nicheRover/CC outputs/new CC fixed/new run 02 10 Feb 2020/output5_5_100k.RData")
ndim5cd <- as.numeric(output5)
load("C:/Users/francis van Oordt/OneDrive - McGill University/Documents/McGill/00Res Prop v2/Chap 0 BC isotopes/nicheRover/CC outputs/new CC fixed/new run 02 10 Feb 2020/output12_5_100k.RData")
ndim5 <- as.numeric(output12)

load("C:/Users/francis van Oordt/OneDrive - McGill University/Documents/McGill/00Res Prop v2/Chap 0 BC isotopes/nicheRover/CC outputs/new CC fixed/run new metrics Apr 2021/matrixrange_5_100k.RData")
range5d<-matrixrange

load("C:/Users/francis van Oordt/OneDrive - McGill University/Documents/McGill/00Res Prop v2/Chap 0 BC isotopes/nicheRover/CC outputs/new CC fixed/run new metrics Apr 2021/samplessd_5_100k.RData")
sdnnd5d<-samplessd

#################################################
#matrix.ndim <- matrix(data = c(1:100), ncol = 5)
#numeric.ndim <- as.numeric(matrix.ndim)
#ndim2 <- numeric.ndim
#ndim3 <- seq(2,200,2)
#ndim5 <- c(400:499)

ndim.box <- list(ndim1, ndim2, ndim3, ndim5)
ndim.box.cd<-list(ndim1cd, ndim2cd, ndim3cd, ndim5cd)
ndim.box.sdnnd<-list(sdnnd1d, sdnnd2d, sdnnd3d, sdnnd5d)

par(mfrow=c(2,2), mar = c(4, 4, 1, 1))

#par(mar = c(6, 7, 4, 8) + 0.1)
boxplot(ndim.box.cd, axes = FALSE, ylim = c(0,c(range(ndim.box.cd)[2]+0.5)))
axis(side = 1, at = 1:4, labels = c("1d","2d", "3d", "5d"))
axis(side = 2, at = seq(0,c(range(ndim.box.cd)[2]+0.5),4))
mtext("CD", side = 2, line = 3, cex = 1, font = 1)
mtext("A", side = 2, line = -1, at = trunc(c(range(ndim.box.cd)[2]+0.5)), las= 2)
box()

#pltCD <- recordPlot()


#par(mar = c(6, 7, 4, 8) + 0.1)
boxplot(ndim.box, axes = FALSE, ylim = c(0,c(range(ndim.box)[2]+0.5)))
axis(side = 1, at = 1:4, labels = c("1d","2d", "3d", "5d"))
axis(side = 2, at = seq(0,c(range(ndim.box)[2]+0.5)))
mtext("NND", side = 2, line = 3, cex = 1, font = 1)
mtext("B", side = 2, line = -1, at = trunc(c(range(ndim.box)[2]+0.5)), las= 2)
box()
#pltNND <- recordPlot()

#par(mar = c(4, 7, 4, 8) + 0.1)
boxplot(ndim.box.sdnnd, axes = FALSE, ylim = c(0,c(range(ndim.box.sdnnd)[2]+0.5)))
axis(side = 1, at = 1:4, labels = c("1d","2d", "3d", "5d"))
axis(side = 2, at = seq(0,c(range(ndim.box.sdnnd)[2]+0.5)))
mtext("SDNND", side = 2, line = 3, cex = 1, font = 1)
mtext("C", side = 2, line = -1, at = trunc(c(range(ndim.box.sdnnd)[2]+0.5)), las= 2)
box()
#pltSNND <- recordPlot()


#ranges 5d
#par(mar = c(4, 5, 4, 2) + 0.1)
boxplot(range5d[,,1], range5d[,,2], range5d[,,3],  range5d[,,4], range5d[,,5],axes = FALSE, ylim = c(0,c(range(range5d)[2]+0.5)))
axis(side = 1, at = 1:5, labels = c("PC1Css", "PC2Css","PC1Nss", "PC2Nss","stdDSul"))
axis(side = 2, at = seq(0,c(range(range5d)[2]+0.5),4))
mtext("Dimension Range", side = 2, line = 3, cex = 1, font = 1)
mtext("D", side = 2, line = -1, at = trunc(c(range(range5d)[2]+0.5)), las= 2)
box()

dev.copy(jpeg, filename = "multigraphmetrics.jpg", width = 3500, height = 2000, 
         res = 300)
dev.off()



p <- ggplot(ndim.box.sdnnd) + 
  geom_boxplot()

pltRange <- recordPlot()

