library(tidyverse)
library(MuMIn)
library(MASS)

#read dataframe with 
isotopes.dat<-read_csv("MM.csv", col_types = cols(Species = col_factor(levels = c("PECO", 
                                                                                "DCCO", "RHAU", "ANMU", "LSPE"))))

#Full model with 3 PCs for each N and C and std Sulfur

isotopes.lda.cv7<-lda( Species ~ PC1Css + PC2Css + PC3Css + PC1Nss + PC2Nss + PC3Nss + stdDelta.S , data = isotopes.dat, CV = TRUE)

isotopes.lda7<-lda( Species ~ PC1Css + PC2Css + PC3Css + PC1Nss + PC2Nss + PC3Nss + stdDelta.S , data = isotopes.dat)


plot(isotopes.lda7, dimen = 2)

ct7<-table(isotopes.dat$Species, isotopes.lda.cv7$class)

#proportion of correctly classified groups
diag(prop.table(ct7, 1))

sum(diag(prop.table(ct7)))


#model with 5 dimension 1 and 2 PCAs, and std Sulfur
isotopes.lda.cv5<-lda( Species ~ PC1Css + PC2Css + PC1Nss + PC2Nss + stdDelta.S , data = isotopes.dat, CV = TRUE)

isotopes.lda5<-lda( Species ~ PC1Css + PC2Css + PC1Nss + PC2Nss + stdDelta.S , data = isotopes.dat)


plot(isotopes.lda5, dimen = 2)

ct5<-table(isotopes.dat$Species, isotopes.lda.cv5$class)

#proportion of correctly classified groups
diag(prop.table(ct5, 1))

sum(diag(prop.table(ct5)))

#model with 3 dim PCA 1 for C and N, and  stdSulfur
isotopes.lda.cv3<-lda( Species ~ PC1Css + PC1Nss + stdDelta.S , data = isotopes.dat, CV = TRUE)

isotopes.lda3<-lda( Species ~ PC1Css + PC1Nss + stdDelta.S , data = isotopes.dat)


plot(isotopes.lda3, dimen = 2)

ct2<-table(isotopes.dat$Species, isotopes.lda.cv3$class)

#proportion of correctly classified groups
diag(prop.table(ct3, 1))

sum(diag(prop.table(ct3)))


#model with only 2 dim 1 PCA for N and C

isotopes.lda.cv2<-lda( Species ~ PC1Css + PC1Nss , data = isotopes.dat, CV = TRUE)

isotopes.lda2<-lda( Species ~ PC1Css + PC1Nss , data = isotopes.dat)


plot(isotopes.lda2, dimen = 2)

ct2<-table(isotopes.dat$Species, isotopes.lda.cv2$class)

#proportion of correctly classified groups
diag(prop.table(ct2, 1))

sum(diag(prop.table(ct2)))


