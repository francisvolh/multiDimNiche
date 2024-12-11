library(nicheROVER)
set.seed(4321)
MM<-read.csv("MM.csv") #after running PCA and cleaning up columns in datasheet

nsamples<-100000
ndim.vector<-c("1d", "2d", "3d", "5d")

system.time({
  all.dim.over<-NULL
  for (i in ndim.vector) {
    print(i)
    
    if (i == "1d") {
      dataCNS<-MM[,c("Species","PC1Nss")]
      BCiso.par <- tapply(1:nrow(dataCNS), dataCNS$Species,
                          function(ii) niw.post(nsamples = nsamples,
                                                X = dataCNS[ii, 2:ncol(dataCNS)])) 
      over.stat <- overlap(BCiso.par, nreps = nsamples, nprob = 1000, alpha = c(0.95, 
                                                                                0.99))
      # The mean overlap metrics calculated across iteratations for both niche
      # region sizes (alpha = .95 and alpha = .99) can be calculated and displayed
      # in an array.
      over.mean1 <- apply(over.stat, c(1:2, 4), mean) * 100
      round(over.mean1, 2)
      over.mean<-as.matrix(over.mean1[,,1])
      overmean<-as.numeric(overmean)
      Dim1<-overmean
      
      #print(Dim1)
      all.dim.over <- cbind(all.dim.over, Dim1)
      
    }else{
      
      if (i == "2d") {
        dataCNS<-MM[,c("Species","PC1Css","PC1Nss")]
        BCiso.par <- tapply(1:nrow(dataCNS), dataCNS$Species,
                            function(ii) niw.post(nsamples = nsamples,
                                                  X = dataCNS[ii, 2:ncol(dataCNS)])) 
        over.stat <- overlap(BCiso.par, nreps = nsamples, nprob = 1000, alpha = c(0.95, 
                                                                                  0.99))
        # The mean overlap metrics calculated across iteratations for both niche
        # region sizes (alpha = .95 and alpha = .99) can be calculated and displayed
        # in an array.
        over.mean2 <- apply(over.stat, c(1:2, 4), mean) * 100
        round(over.mean2, 2)
        overmean<-as.matrix(over.mean2[,,1])
        overmean<-as.numeric(overmean)
        Dim2<-overmean
        #print(Dim2)
        all.dim.over <- cbind(all.dim.over, Dim2)
        
        
      }else{
        if (i == "3d"){
          dataCNS<-MM[,c("Species","PC1Css","PC1Nss","stdDelta.S")]
          BCiso.par <- tapply(1:nrow(dataCNS), dataCNS$Species,
                              function(ii) niw.post(nsamples = nsamples,
                                                    X = dataCNS[ii, 2:ncol(dataCNS)]))
          over.stat <- overlap(BCiso.par, nreps = nsamples, nprob = 1000, alpha = c(0.95, 
                                                                                    0.99))
          # The mean overlap metrics calculated across iteratations for both niche
          # region sizes (alpha = .95 and alpha = .99) can be calculated and displayed
          # in an array.
          over.mean3 <- apply(over.stat, c(1:2, 4), mean) * 100
          round(over.mean3, 2)
          overmean<-as.matrix(over.mean3[,,1])
          overmean<-as.numeric(overmean)
          Dim3<-overmean
          #print(Dim3)
          all.dim.over <- cbind(all.dim.over, Dim3)
          
          
        }else{
          if (i == "5d") {
            dataCNS<-MM[,c("Species","PC1Css","PC2Css",
                           "PC1Nss","PC2Nss",
                           "stdDelta.S")]#subset with 5 dimension for runs
            BCiso.par <- tapply(1:nrow(dataCNS), dataCNS$Species,
                                function(ii) niw.post(nsamples = nsamples,
                                                      X = dataCNS[ii, 2:ncol(dataCNS)]))
            over.stat <- overlap(BCiso.par, nreps = nsamples, nprob = 1000, alpha = c(0.95, 
                                                                                      0.99))
            # The mean overlap metrics calculated across iteratations for both niche
            # region sizes (alpha = .95 and alpha = .99) can be calculated and displayed
            # in an array.
            over.mean5 <- apply(over.stat, c(1:2, 4), mean) * 100
            round(over.mean5, 2)
            overmean<-as.matrix(over.mean5[,,1])
            overmean<-as.numeric(overmean)
            Dim5<-overmean
            #print(Dim5)
            all.dim.over <- cbind(all.dim.over, Dim5)
            
            
          }
        }
      } 
    }
    
  }
})

#saveRDS(all.dim.over, "all.dim.over.1-2-3-5dims.RDS")
all.dim.over <- readRDS("C:/Users/francis van oordt/OneDrive - McGill University/Documents/McGill/00Res Prop v2/Chap 0 BC isotopes/all.dim.over.1-2-3-5dims.RDS")
df <- as.data.frame(all.dim.over)
df$samples <- row.names(df)


colname <- sort(unique(MM$Species))
rowname <- sort(unique(MM$Species))

labels.1 <- NULL
for (i in 1:length(rowname)) {
  for (j in 1:length(colname)) {
    label <- paste0(colname[j], "-", rowname[i])
    
    labels.1 <- c(labels.1, label) 
    
  }
}


df$labels1 <- labels.1

df$labels1 <- factor(df$labels1, levels = df$labels1) #lock in the factor order!!!!

df<-df|>
  dplyr::filter(!is.na(Dim1))

#without a legent because I dont have a grouping var (I think)

plot<-df|>
  tidyr::pivot_longer(cols = c(Dim1, Dim2, Dim3, Dim5), names_to = "Dimensions")|>
  ggplot2::ggplot(ggplot2::aes(color = Dimensions))+
  ggplot2::geom_point(ggplot2::aes(x=labels1, y = value), cex = 3, alpha = 0.5)+
  ggplot2::theme_bw()+
  ggplot2::xlab(NULL)+
  ggplot2::ylab("Percentage of directional overlap")+
  ggplot2::theme(axis.text.x= ggplot2::element_text(angle=90))+
  ggplot2::scale_color_manual( labels = c("1D", "2D", "3D", "5D"), 
                               values = c("green","red", "blue", "darkgray"))+
  ggplot2::theme(legend.position = "top")

plot

ggplot2::ggsave(plot, filename = "plots/plot_overlap_1235dims.png", width=10, height = 5, units = "in")


     ##### old version attempt
#df|>
#  ggplot2::ggplot()+
#  ggplot2::geom_point(ggplot2::aes(x=labels1, y = Dim1), cex = 3, alpha = 0.5, color = "green")+
 # ggplot2::geom_point(ggplot2::aes(x=labels1, y = Dim2), cex = 3, alpha = 0.5, color = "red")+
  #ggplot2::geom_point(ggplot2::aes(x=labels1, y = Dim3), cex = 3, alpha = 0.5, color = "blue")+
  #ggplot2::geom_point(ggplot2::aes(x=labels1, y = Dim5), cex = 3, alpha = 0.5, color = "darkgray")+
  #ggplot2::theme_bw()+
  #ggplot2::xlab(NULL)+
  #ggplot2::ylab("Percent of directional overlap")+
  #ggplot2::theme(axis.text.x= ggplot2::element_text(angle=90))+
  #ggplot2::scale_color_manual( values = "1D", "2D", "3D", "5D")

#ggsave(plot, "plot5dimsoverlap.png" )  
