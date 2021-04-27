
M<-merSIABAAxxx #cleaned database with AA and C N S
M$dNGluPhe<-NA
M$dNGluPhe<-M$dN_Glu-M$dN_Phe #trophic level
dNsource<-NA
M$dNsource<-(M$dN_Gly+M$dN_Phe)/2 #average source

M<-M[,c("Species","Delta.15N","d13Ccor","Delta.S", 
        "dNGluPhe", "dNsource", "decade",
        "dC_Ala","dC_Val","dC_Gly","dC_Ile","dC_Leu","dC_Pro","dC_Asp","dC_Phe",
        "dC_Glu","dC_Lys", "dN_Ala", "dN_Asp" ,"dN_Glu","dN_Gly", "dN_Ile",
        "dN_Leu","dN_Lys","dN_Phe","dN_Pro","dN_Val" , "Year")]

MM <-M[complete.cases(M[c("Delta.S")]),]#after this, the gaps in SULFUR are gone, SMALLEST dataset
MM$Species<-as.factor(MM$Species)
MM<-droplevels(MM)
table(MM$Species)
##PCA of all AASIA Cs and Ns ###################

#CHOOSE SET OF AASIA: all carbons
res.pca <- prcomp(MM[, c("dC_Ala","dC_Val","dC_Gly","dC_Ile","dC_Leu","dC_Pro",
                         "dC_Asp","dC_Phe","dC_Glu","dC_Lys")], scale = TRUE)

#Append PC1 score to data ALL Css
MM$PC1Css <- scores(res.pca, choices=c(1))
MM$PC2Css <- scores(res.pca, choices=c(2))
MM$PC3Css <- scores(res.pca, choices=c(3))

#run and append PCs for ALL Nss
res.pca <- prcomp(MM[, c("dN_Ala","dN_Asp","dN_Glu","dN_Gly",
                         "dN_Ile","dN_Leu","dN_Lys","dN_Phe",
                         "dN_Pro","dN_Val" )], scale = TRUE)

MM$PC1Nss <- scores(res.pca, choices=c(1))
MM$PC2Nss <- scores(res.pca, choices=c(2))
MM$PC3Nss <- scores(res.pca, choices=c(3))

########Standardizing variables##################################################
MM$stddNGluPhe<-NULL
MM$stddNGluPhe<-(MM$dNGluPhe-(mean(MM$dNGluPhe)))/sd(MM$dNGluPhe)

MM$stddNsource<-NULL
MM$stddNsource<-(MM$dNsource-(mean(MM$dNsource)))/sd(MM$dNsource)

MM$stdDelta.S<-NULL
MM$stdDelta.S<-(MM$Delta.S-(mean(MM$Delta.S)))/sd(MM$Delta.S)
#write.csv(MM, "MM.csv")
#MM<-read.csv("MM.csv")
