MM$Glx<-MM$dN_Glu-MM$dN_Lys
MM$Gl2<-MM$dN_Glu-MM$dN_Phe


MM$Gluscaled<-MM$dN_Glu-(mean(MM$dN_Glu))/sd(MM$dN_Glu)

MM$dN_Phescaled<-MM$dN_Phe-(mean(MM$dN_Phe))/sd(MM$dN_Phe)
MM$Lysscaled<-MM$dN_Lys-(mean(MM$dN_Lys))/sd(MM$dN_Lys)

MM$Glxscaled<-MM$Glx-(mean(MM$Glx))/sd(MM$Glx)

MM$Gl2scaled<-MM$Gl2-(mean(MM$Gl2))/sd(MM$Gl2)

MM$Glxiii<-MM$Gluscaled-MM$Lysscaled
MM$Gl2iii<-MM$Gluscaled-MM$dN_Phescaled

MM$TP1 <- ((MM$dN_Glu-MM$dN_Phe -2.42)/5.63)+1

MM$TP2 <- ((((MM$dN_Glu+ MM$dN_Ala+MM$dN_Ile+ MM$dN_Leu+MM$dN_Pro+MM$dN_Val)/6)-MM$dN_Phe -2.42)/5.63)+1


cor(MM$PC1Nss, MM$Glx)

GGally::ggpairs(data.frame(MM[,c("PC1Nss","Gl2scaled", "Glxscaled", "PC1Css")]))


summary(MM)


#install.packages("SuessR")

library(SuessR)

suesscorrections<-MM|>
  dplyr::mutate(region = "Gulf of Alaska")|>
  dplyr::select(id=USOXcod, year=Year.x, dC_Ala,dC_Val,dC_Gly,dC_Ile,dC_Leu,dC_Pro,dC_Asp,dC_Phe,dC_Glu,dC_Lys, region)


caa.names <- c("dC_Ala","dC_Val","dC_Gly","dC_Ile","dC_Leu","dC_Pro","dC_Asp","dC_Phe","dC_Glu","dC_Lys")

df.dCcor <-NULL
for (i in caa.names) {
  xx<-suesscorrections|>
    dplyr::select(id, d13c=i, region, year)|>
    SuessR::SuessR()|>
    dplyr::pull(d13c.cor)
  

  df.dCcor <- cbind(df.dCcor,xx )
}

df.dCcor.df<-as.data.frame(df.dCcor)
names(df.dCcor.df) <- c("dC_Ala_cor","dC_Val_cor","dC_Gly_cor","dC_Ile_cor","dC_Leu_cor","dC_Pro_cor","dC_Asp_cor","dC_Phe_cor","dC_Glu_cor","dC_Lys_cor")

MM<-cbind(MM, df.dCcor.df)

res.pca2 <- prcomp(MM[, c("dC_Ala_cor","dC_Val_cor","dC_Gly_cor","dC_Ile_cor",
                         "dC_Leu_cor","dC_Pro_cor","dC_Asp_cor","dC_Phe_cor",
                         "dC_Glu_cor","dC_Lys_cor")], scale = TRUE)

MM$PC1CssCor <- as.numeric(scores(res.pca2, choices=c(1)))
MM$PC2CssCor <- as.numeric(scores(res.pca2, choices=c(2)))
MM$PC3CssCor <- as.numeric(scores(res.pca2, choices=c(3)))
names(MM)
GGally::ggpairs(data.frame(MM[,c("PC1Nss","PC2Nss", #"PC3Nss", 
                                 "PC1Css", "PC1CssCor",
                          "Gl2",  "Glx", "TP1", "TP2")]))
