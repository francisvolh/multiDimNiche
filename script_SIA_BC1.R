#library(stringr)
#library(expss)
library(RColorBrewer)
library(Rmisc)

setwd("C:/Users/Francis van Oordt/OneDrive - McGill University/Documents/McGill/00Res Prop v2/Chap 0 BC isotopes")
#Load SIA for C, N and S, and also deh AASIA
wkSIA <- read.csv("WK file Seabird Hg.csv", header = T, na.strings=c("", "NA"))
AASIA<- read.csv("AASIA.csv", header = TRUE) #species code has been split to merge with other SIA dataframe

#produce a column to save species code again to save it for later after mergin
AASIA$code1<-AASIA$code3


#MERGING 2 dataframes of SIA bulk and AA
merSIA<-merge(wkSIA, AASIA, by.x="USOXcod", by.y="code3", all=TRUE)

#apply correction factor for dCarbon
merSIA$d13Ccor <- NA
merSIA$d13Ccor <-  merSIA$Delta.13C - 4.46 + (7.32 * log10(merSIA$C.Amount..ug./merSIA$N.Amount..ug.))

##TO FIX LOCATIONS NAMES
##change labels for locations to "Name" Island instead of "Name Is."#
levels(merSIA$Location.x)[levels(merSIA$Location.x)=="Langara Is."] <- "Langara Island"
merSIA$Location.x[merSIA$Location.x == 'Langara Is.'] <- 'Langara Island'
levels(merSIA$Location.x)[levels(merSIA$Location.x)=="Reef Is."] <- "Reef Island"#adding a level that doesnt exist
merSIA$Location.x[merSIA$Location.x == 'Reef Is.'] <- 'Reef Island'
####
#to check extra values
#merSIA[is.na(merSIA$d13Ccor),]
#
nrow(merSIA[!is.na(merSIA$Delta.S),])## to check HOW MANY EXTRA values!!!! change in which variable we want to check
##
## first visualization of data
plot(merSIA$d13Ccor,merSIA$Delta.15N)#outlier detected
plot(merSIA$Delta.S,merSIA$Delta.15N)#wont work because delta S is character, solve later
plot(merSIA$Delta.S,merSIA$d13Ccor)#wont work because dS is character, sol later
##

##
###to delete outlier, first check for the value and position to delete row
#max(merSIA$Delta.15N, na.rm = T)
merSIA[which(merSIA$Delta.15N>= 64.04 & merSIA$Delta.15N<= 64.04),]     
#merSIA[204,]
#deleting row CONFIRM ROW NUMBER!!!
merSIA <- merSIA[-c(204),] 
####

######FIX/correct character and factors into Numeric when needed
merSIA$Year.x<-as.factor(as.character(merSIA$Year.x))
merSIA$Year.y<-as.factor(as.character(merSIA$Year.y))

#Fix character value error in dS and d34S(AA-SIA) and replace "-" for NA for 
merSIA$Delta.S<-as.numeric(as.character(merSIA$Delta.S))#this we WILL USE
####
summary(merSIA$Delta.S)
###
merSIA$d34S<-as.numeric(as.character(merSIA$d34S))#we WONT USE THIS ONE ANYWAY
summary(merSIA$d34S)

##check if there are LOCATION values with NA for certain variables that were duplicated when merging
nrow(merSIA[is.na(merSIA$Species),])
nrow(merSIA[is.na(merSIA$Location.x),])#switch to different duplicated variables
merSIA[is.na(merSIA$Location.x),]
#replacing values from duplicated columns specially SPECIES and LOCATIONS from the dif dbs
merSIA$Species[is.na(merSIA$Species)] <- merSIA$SP[is.na(merSIA$Species)]
merSIA$Location.x[is.na(merSIA$Location.x)] <- merSIA$Location.y[is.na(merSIA$Location.x)]
merSIA$Year.x[is.na(merSIA$Year.x)] <- merSIA$Year.y[is.na(merSIA$Year.x)]

#produce a Decade variable
merSIA$Year<-NA
merSIA$Year<-as.numeric(as.character(merSIA$Year.x))
merSIA$decade<-NA 
merSIA$decade[merSIA$Year >= 1950 & merSIA$Year<1990] <- "1970-1989"
merSIA$decade[merSIA$Year >= 1990] <- "1990-2006"


merSIA$decade<-as.factor(merSIA$decade)

#split some species label by region
merSIA$labSPgeo<-NA
merSIA$labSPgeo[merSIA$Location.x == 'Hippa Island' & merSIA$Species == 'LSPE'] <- 'LSPE.N'
merSIA$labSPgeo[merSIA$Location.x == 'Storm Island' & merSIA$Species == 'LSPE'] <- 'LSPE.M'
merSIA$labSPgeo[merSIA$Location.x == 'Thorton Island' & merSIA$Species == 'LSPE'] <- 'LSPE.S'
merSIA$labSPgeo[merSIA$Location.x == 'Thomas Island' & merSIA$Species == 'LSPE'] <- 'LSPE.S'
merSIA$labSPgeo[merSIA$Location.x == 'Cleland Island' & merSIA$Species == 'LSPE'] <- 'LSPE.S'

merSIA$labSPgeo[merSIA$Location.x == 'Cleland Island' & merSIA$Species == 'RHAU'] <- 'RHAU.S'
merSIA$labSPgeo[merSIA$Location.x == 'Lucy Island' & merSIA$Species == 'RHAU'] <- 'RHAU.N'
merSIA$labSPgeo[merSIA$Location.x == 'Pine Island' & merSIA$Species == 'RHAU'] <- 'RHAU.M'


merSIA$labSPgeo[merSIA$Species == 'DCCO'] <- 'DCCO'
merSIA$labSPgeo[merSIA$Species == 'PECO'] <- 'PECO'
merSIA$labSPgeo[merSIA$Species == 'GBHE'] <- 'GBHE'
merSIA$labSPgeo[merSIA$Species == 'ANMU'] <- 'ANMU'

merSIA$labSPgeo<-as.factor(merSIA$labSPgeo)
summary(merSIA$labSPgeo)
summary(merSIA$decade)
####
####Backing up
merSIAxx<-merSIA
#recover from backup
#merSIA<-merSIAxx

##subset deleting the  collumns I dont need
##names(merSIA)
#c("Project..","SIA..Lab.d13C...d15N.Analysis..","NWRC..Lab.ID.","Lab.Name","Tray.Name","Well.Id","d15N.Comment","d13C...d15N.Dry.Wt...mg.","SIA.Lab.dS.Analysis..","Obs")
merSIA[,c("Project..","SIA..Lab.d13C...d15N.Analysis..","Specimen...or...in.pool","dS.Dry.Weight..mg.",
          "NWRC..Lab.ID.","Lab.Name","Tray.Name","Well.Id","d15N.Comment",
          "d13C...d15N.Dry.Wt...mg.","SIA.Lab.dS.Analysis..","SIA.Avg..of.Replicate.","Obs", "Location.y", "SP", "Year.y")] <- list(NULL)

#SAVE FILE
#write.csv(merSIA, file="merSIALAST.csv")


##read final file LOAD LOAD LOAD
#merSIA<-read.csv("merSIALAST.csv")#if load this file, check for all FACTORS

##produce dataframe with NO NAs
#delete NAs for d13C
merSIACN <-merSIA[complete.cases(merSIA[c("d13Ccor","Delta.15N")]),]
#delete NAs for dS
merSIACNS <-merSIACN[complete.cases(merSIACN[c("Delta.S")]),]
#delete GBHE
merSIACNSnGB<-merSIACNS[merSIACNS$Species %in% c("ANMU","LSPE", 
                                                 "PECO", "RHAU", "DCCO"),]

#delete NAs for AASIA after CNS clean
merSIABAA <-merSIACNS[complete.cases(merSIACNS[c("dC_Ala")]),]

#delete NAs for AASIA disregarding CSN
merSIABAAxxx<-merSIA[complete.cases(merSIA[c("dC_Ala")]),]

merSIA <-droplevels(merSIA)
merSIACN <-droplevels(merSIACN)
merSIACNS <-droplevels(merSIACNS)
merSIABAA <-droplevels(merSIABAA)
merSIABAAxxx<-droplevels(merSIABAAxxx)

#database for S disregarding CN 
merSIAS <-merSIA[complete.cases(merSIA[c("Delta.S")]),]


#database for AA andS  disregarding CN
merSIAamS <-merSIAS[complete.cases(merSIAS[c("dC_Ala")]),]
table(merSIAamS$Species, merSIAamS$decade)
table(merSIAamS$Species)
merSIAamS<-droplevels(merSIAamS)

#test
