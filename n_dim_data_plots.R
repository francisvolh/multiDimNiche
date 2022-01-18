###dummy data as a matrix/ array from CSVs produced 3D

dataA <- read.csv(file = "dummydataAiii.csv", stringsAsFactors = FALSE)
dataB <- read.csv(file = "dummydataBiii.csv", stringsAsFactors = FALSE)
dataC <- read.csv(file = "dummydataCiii.csv", stringsAsFactors = FALSE)

newdatarray <- array(data = NA, dim = c(10,3,3))

newdatarray[,,1] <- as.matrix(dataA)
newdatarray[,,2] <- as.matrix(dataB)
newdatarray[,,3] <- as.matrix(dataC)

new.matrixi <- newdatarray

matrix.data <- new.matrixi

######################################################
#load dummy data as a dataframe 3D
dummydata<-read.csv("C:/Users/Francis van Oordt/OneDrive - McGill University/Documents/McGill/00Res Prop v2/Chap 0 BC isotopes/dummydata games/dummydata (1).csv", stringsAsFactors = FALSE)

dummydata<- read.csv(file = "dummydatai.csv", stringsAsFactors = FALSE)
dummydata$Species<-as.factor(dummydata$Sp)

#if using the simulator from MBH
dummydata<-dat1.7d
dummydata<-as.data.frame(dummydata)

#plotly plot 3D
library(plotly)
fig <- plot_ly(dummydata, x = ~var1  , y = ~var2 , z = ~var3       , color = ~Sp)
fig <- fig %>% add_markers()
fig <- fig %>% layout(scene = list(xaxis = list(title = 'Var1'),
                                   yaxis = list(title = 'Var2'),
                                   zaxis = list(title = 'Var3')))
fig

############ 3d Graph with RGL

mycolors <- c('royalblue1', 'darkcyan', 'oldlace', 'green', 'red')
mycolors <- mycolors[ as.numeric(factor(dummydata$Sp)) ]

open3d()

rgl::plot3d(dummydata$var1, dummydata$var2, dummydata$var3, 
            col = mycolors, type = 's', radius = 1,)
rgl::aspect3d(1, 1, 1)
rgl::legend3d("topright", 
              legend = paste('Species', sort(unique(dummydata$Sp))), 
              pch = 16, col = sort(unique(mycolors)), 
              cex=1, inset=c(0.02))

mean(dummydata$var1)
mean(dummydata$var2)
mean(dummydata$var3)
mean(dummydata$var4)



###########plotly##########

fig <- plot_ly(dummydata, x = ~var2, y = ~var3, z = ~var4, color = ~Sp)
fig <- fig %>% add_markers()
fig <- fig %>% layout(scene = list(xaxis = list(title = 'Var2'),
                                   yaxis = list(title = 'Var3'),
                                   zaxis = list(title = 'Var4')))
fig


fig <- plot_ly(dummydata, x = ~var1, y = ~var3, z = ~var4, color = ~Sp)
fig <- fig %>% add_markers()
fig <- fig %>% layout(scene = list(xaxis = list(title = 'Var1'),
                                   yaxis = list(title = 'Var3'),
                                   zaxis = list(title = 'Var4')))
fig

fig <- plot_ly(dummydata, x = ~var1, y = ~var2, z = ~var4, color = ~Sp)
fig <- fig %>% add_markers()
fig <- fig %>% layout(scene = list(xaxis = list(title = 'Var1'),
                                   yaxis = list(title = 'Var2'),
                                   zaxis = list(title = 'Var4')))
fig

fig <- plot_ly(dummydata, x = ~var5, y = ~var6, z = ~var7, color = ~Sp)
fig <- fig %>% add_markers()
fig <- fig %>% layout(scene = list(xaxis = list(title = 'Var1'),
                                   yaxis = list(title = 'Var2'),
                                   zaxis = list(title = 'Var4')))
fig
###########if coming from simulation of MBH and made dataframe##########
fig <- plot_ly(dummydata, x = ~data.V3, y = ~data.V4, z = ~data.V5, color = ~data.Group)
fig <- fig %>% add_markers()
fig <- fig %>% layout(scene = list(xaxis = list(title = 'Var1'),
                                   yaxis = list(title = 'Var2'),
                                   zaxis = list(title = 'Var4')))
fig
#########################################################################

###########if coming from 7D df ########################################
fig <- plot_ly(C, x = ~Delta.S, y = ~dNGluPhe, z = ~PC1ess, color = ~Species)
fig <- fig %>% add_markers()
fig <- fig %>% layout(scene = list(xaxis = list(title = 'Var1'),
                                   yaxis = list(title = 'Var2'),
                                   zaxis = list(title = 'Var4')))
fig

fig <- plot_ly(C, x = ~dNsource, y = ~PC1nones, z = ~Delta.S, color = ~Species)
fig <- fig %>% add_markers()
fig <- fig %>% layout(scene = list(xaxis = list(title = 'Var1'),
                                   yaxis = list(title = 'Var2'),
                                   zaxis = list(title = 'Var4')))
fig
#########################################################################

dummydata<- read.csv(file = "dummydataUneven1.csv", stringsAsFactors = FALSE)

hv2dsimU <- fitMBH(dummydataUN1, groups = "Sp", vars = c("var1", "var2"))
hv3dsimU <- fitMBH(dummydataUN2, groups = "Sp", vars = c("var1" ,"var2", "var3"))
hv7dsimU <- fitMBH(dummydataUN3, groups = "Sp", vars = c("var1" ,"var2" ,"var3", "var4", "var5" ,"var6","var7"))
