install.packages("plotly")
library(plotly)

#be sure to use 5d data source from NicheRover

matrix.data <- output2[sample(nrow(output2), 1000), ,]

sort(unique(dataCNS$Species))
colnames(dataCNS)

df1<-as.data.frame(matrix.data[,,1])
df2<-as.data.frame(matrix.data[,,2])
df3<-as.data.frame(matrix.data[,,3])
df4<-as.data.frame(matrix.data[,,4])
df5<-as.data.frame(matrix.data[,,5])

ANMU<-df1$V1
DCCO<-df1$V2
LSPE<-df1$V3
PECO<-df1$V4
RHAU<-df1$V5
dff1<-data.frame(ANMU,DCCO,LSPE,PECO,RHAU)
sdff1 <- stack(dff1)

ANMU<-df2$V1
DCCO<-df2$V2
LSPE<-df2$V3
PECO<-df2$V4
RHAU<-df2$V5
dff2<-data.frame(ANMU,DCCO,LSPE,PECO,RHAU)
sdff2 <- stack(dff2)

ANMU<-df3$V1
DCCO<-df3$V2
LSPE<-df3$V3
PECO<-df3$V4
RHAU<-df3$V5
dff3<-data.frame(ANMU,DCCO,LSPE,PECO,RHAU)
sdff3 <- stack(dff3)

ANMU<-df4$V1
DCCO<-df4$V2
LSPE<-df4$V3
PECO<-df4$V4
RHAU<-df4$V5
dff4<-data.frame(ANMU,DCCO,LSPE,PECO,RHAU)
sdff4 <- stack(dff4)

ANMU<-df5$V1
DCCO<-df5$V2
LSPE<-df5$V3
PECO<-df5$V4
RHAU<-df5$V5
dff5<-data.frame(ANMU,DCCO,LSPE,PECO,RHAU)
sdff5 <- stack(dff5)

bounddff<-cbind(sdff1, sdff2, sdff3, sdff4, sdff5)
head(bounddff)
bounddff$ind<-NULL
bounddff$ind<-NULL
bounddff$ind<-NULL
bounddff$ind<-NULL

colnames(bounddff)<-c(colnames(dataCNS)[2:6], "Species")



bounddff$PC2Nss2<-bounddff$PC2Nss+10
fig2<-plot_ly(bounddff, 
             x = ~PC1Css, 
             y = ~PC1Nss, 
             z = ~stdDelta.S, 
             color = ~PC2Css, size = ~PC2Nss2, sizes = c(5,450))
fig2

#################

fig<-plot_ly(bounddff, 
              x = ~PC1Css, 
              y = ~PC1Nss, 
              z = ~stdDelta.S, 
              color = ~Species, size = 1,
              colors = c("red", "blue", "purple","black", "green"))
fig
