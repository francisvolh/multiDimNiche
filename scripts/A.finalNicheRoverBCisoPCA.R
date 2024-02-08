#FINAL VERSION OF CODE TO RUN NICHE ROVER for 1D 2D, 3D, 5D
# NOTE: run time for 5D would take about 30 hours

# Currently set up to  select desired run, with the ndim object,
# must customize for desired run type and time


library(factoextra)
library(vegan)
library(nicheROVER)
library(ggplot2)
library(ggpubr)
library(hdrcde)

#########################################################################
#Call raw subset data 

MM<-read.csv("raw_data/AASIA_raw_data.csv") # read object with original filtered dataframe to run PCA on

##PCA of all AASIA Cs and Ns ###################

#CHOOSE SET OF AASIA: all carbons
res.pca <- prcomp(MM[, c("dC_Ala","dC_Val","dC_Gly",
                         "dC_Ile","dC_Leu","dC_Pro",
                         "dC_Asp","dC_Phe","dC_Glu",
                         "dC_Lys")], scale = TRUE)

#Append PC1 score to data ALL Css
MM$PC1Css <- as.numeric(scores(res.pca, choices=c(1)))
MM$PC2Css <- as.numeric(scores(res.pca, choices=c(2)))
MM$PC3Css <- as.numeric(scores(res.pca, choices=c(3)))

#run and append PCs for ALL Nss
res.pca <- prcomp(MM[, c("dN_Ala","dN_Asp","dN_Glu","dN_Gly",
                         "dN_Ile","dN_Leu","dN_Lys","dN_Phe",
                         "dN_Pro","dN_Val" )], scale = TRUE)

MM$PC1Nss <- as.numeric(scores(res.pca, choices=c(1)))
MM$PC2Nss <- as.numeric(scores(res.pca, choices=c(2)))
MM$PC3Nss <- as.numeric(scores(res.pca, choices=c(3)))

########Standardizing variables##################################################

MM$stdDelta.S<-NULL
MM$stdDelta.S<-(MM$Delta.S-(mean(MM$Delta.S)))/sd(MM$Delta.S)

#write.csv(MM, "raw_data/AASIA_data.csv")
#####
set.seed(4321)
nsamples<-100000

ndim<-1 ## Choose dimensions to run model which will be the same as nisotope, but I want to keep them separate

# nicheRover model will run in line 112

#ndim.vector<-c("1d", "2d", "3d", "5d")
##can code as to loop in 1 go over all the dimensions selected - TO DO 

##
##

#Libraries
library(abind)
library(hdrcde)
library(nicheROVER)

print(paste("Started at", Sys.time()))
#######################################################
#functions for Rover metrics

Mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

niche.size <- function(Sigma, alpha = .95) {
  n <- nrow(Sigma)
  sz <- as.numeric(determinant(Sigma, log = TRUE)$modulus)
  sz <- .5 * (sz + n * (log(pi) + log(qchisq(alpha, df = n)))) - lgamma(.5*n+1)
  exp(sz)
}


post.range <- function(post.obj, Dmult) {
  Nobs <- dim(post.obj$Sigma)[3]
  Range.B <- matrix(NA, nrow=Nobs, ncol=2)
  for(j in 1:Nobs) {
    vary <- post.obj$Sigma [1,1,j]
    varx <- post.obj$Sigma [2,2,j]
    meany <- post.obj$mu[j,1]
    meanx <- post.obj$mu[j,2]
    sdy <- sqrt(vary)
    sdx <- sqrt(varx)
    
    bottom.y=meany-(sdy*Dmult)
    top.y=meany+(sdy*Dmult)
    left.x=meanx-(sdx*Dmult)
    right.x=meanx+(sdx*Dmult)
    
    Yrange=abs(top.y-bottom.y)
    Xrange=abs(left.x-right.x)
    
    Range.B[j,] <- cbind(Yrange, Xrange)
  }
  colnames(Range.B) <-c("Y Range", "X Range")
  return(Range.B) 
}


stat.par <- function(post.est, stat=median, Dmult=2,
                     alpha=0.95, cred=c(0.05, 0.95)){ 
  
  modes <- matrix(NA, nrow = length(post.est), ncol=16)
  #rown <- list()
  yname <- colnames(post.est[[2]]$mu)[1]
  xname <- colnames(post.est[[2]]$mu)[2]
  nn=Filter(Negate(is.null), post.est)   #extract all non-null elements of the list
  
  for(i in 1:length(nn)){
    
    #First calculate the Modes, then credible intervals
    tmpn <- apply(nn[[i]]$Sigma, 3, niche.size, alpha=alpha)
    ns <- hdr(tmpn, prob=95, all.modes=F)$mode #calculate mode of the 95% highest density region
    ns.cred <- quantile(tmpn, prob = cred)
    lower.ns <- as.numeric(ns.cred[1])
    upper.ns <- as.numeric(ns.cred[2])
    
    mu.hdr <- apply(nn[[i]]$mu, 2, hdr, prob=95, all.modes=F) #mode of means for x and y
    y.mode <- as.numeric(mu.hdr[[1]][2])
    x.mode <- as.numeric(mu.hdr[[2]][2])
    mu.cred <- apply(nn[[i]]$mu, 2, quantile, prob = cred, na.rm = TRUE)
    lower.mu.y <- mu.cred[1,1]
    upper.mu.y <- mu.cred[2,1]
    lower.mu.x <- mu.cred[1,2]
    upper.mu.x <- mu.cred[2,2]
    
    rng <- post.range(nn[[i]], Dmult)
    rg.hdr <- apply(rng, 2, hdr, prob=95, all.modes=F)
    rgy.mode <- as.numeric(rg.hdr[[1]][2])
    rgx.mode <- as.numeric(rg.hdr[[2]][2])
    mu.cred <- apply(rng, 2, quantile, prob = cred, na.rm = TRUE)
    lower.rgy.y <- mu.cred[1,1]
    upper.rgy.y <- mu.cred[2,1]
    lower.rgy.x <- mu.cred[1,2]
    upper.rgy.x <- mu.cred[2,2]
    
    group <-names(nn[i])
    
    mr <- c(group,
            ns, lower.ns, upper.ns,
            y.mode, lower.mu.y, upper.mu.y,
            x.mode, lower.mu.x, upper.mu.x,
            rgy.mode, lower.rgy.y, upper.rgy.y,
            rgx.mode, lower.rgy.x, upper.rgy.x)
    
    #rown[i] <-names(nn[i])
    modes[i,] <- mr
  }
  colnames(modes) <-c("Group","NS", "NS.Lower", "NS.Upper",
                      paste("Mean",yname, sep="."), paste("Lower.Mean",yname, sep="."), paste("Upper.Mean",yname, sep="."),
                      paste("Mean",xname, sep="."), paste("Lower.Mean",xname, sep="."), paste("Upper.Mean",xname, sep="."),
                      paste("Range",yname, sep="."), paste("Lower.Rng",yname, sep="."), paste("Upper.Rng",yname, sep="."),
                      paste("Range",xname, sep="."), paste("Lower.Rng",xname, sep="."), paste("Upper.Rng",xname, sep="."))
  #rownames(modes) <- rown
  return(modes)
}

#######################################################
#RUN nicherover for corresponding dimension from MM dataframe
print(paste("Processing nicheROver Bayesian model", nsamples, "iterations"))
system.time({
  if (ndim == 1) {
    dataCNS<-MM[,c("Species","PC1Nss")]
    BCiso.par <- tapply(1:nrow(dataCNS), dataCNS$Species,
                        function(ii) niw.post(nsamples = nsamples,
                                              X = dataCNS[ii, 2:ncol(dataCNS)])) 
  }else{
    
    
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
  }
  
})



# Variables to have names on columns ---------------------------------------------------------------
spp <- sort(unique(dataCNS$Species)) # species names

nspp <- length(unique(dataCNS$Species)) # number of species

isotopes <- colnames(dataCNS)[2:ncol(dataCNS)] # isotope label names

nisotope <- length(isotopes) # numner of isotopes/dimensions used

ncol <- 1:nspp # seq of species numbers
nrows <- 1:nsamples # seq of sample numbers

print(paste("model run for", isotopes, "isotopes for", nsamples, "iterations"))
#######################################################
# Then extract the niche metrics for each species
cr.p = c(0.05, 0.95) 
if (ndim != 1) {
  system.time({
    rover.metrics = stat.par(BCiso.par, stat=Mode, Dmult=2,
                             alpha=0.95, cred=cr.p)
    
    print(rover.metrics)
    
  })
}else{"Not possible to run default metrics for 1D"}


#Overlap metrics
system.time({
  over.stat <- overlap(BCiso.par, nreps = nsamples, nprob = 1000, alpha = c(0.95, 
                                                                            0.99))
  
  over.mean <- apply(over.stat, c(1:2, 4), mean) * 100
  round(over.mean, 2)
  
})
print(round(over.mean, 2))

# The mean overlap metrics calculated across iteratations for both niche
# region sizes (alpha = .95 and alpha = .99) can be calculated and displayed
# in an array.

system.time({
  over.cred <- apply(over.stat * 100, c(1:2, 4), quantile, prob = c(0.025, 0.975), 
                     na.rm = TRUE)
  
})
round(over.cred[, , , 1])  # display alpha = .95 niche region
##ACM FVO CODIGUE
#Extract object with Niche Sizes from nicheRover object to Calculate percentiles and graph
nichRobj<-BCiso.par
alpha<-0.95
niches<-NULL
nn=Filter(Negate(is.null), nichRobj)

#loop comes from Homemade function code for niche calculations Martin and Heidi
system.time({
  for (i in 1:length(unique(dataCNS$Species))) {
    tmpn <- apply(nn[[i]]$Sigma, 3, niche.size, alpha=alpha)
    niches<-cbind(niches, tmpn)
  }
})

colnames(niches)<-spp

#for 1D only
if (ndim == 1) {
  oneDimNiche<-NULL
  #for (j in 1:nisotope) {
  # isotopo<-matrix.data[,,j]
  for (i in 1:nspp) {
    mode1 <- mean(niches[,i])
    quantile1 <- quantile(niches[,i],probs=c(0.025,0.5,0.975)) 
    #maximo <- max(isotopo[,i])
    #minimo <- min(isotopo[,i])
    #rango <- max(isotopo[,i]) - min(isotopo[,i])
    modaquantile <- c(mode1,quantile1 
                      #maximo, minimo, rango
    )
    oneDimNiche<-rbind(oneDimNiche, modaquantile)
    
    
  } 
  #}
  oneDimNiche<- as.data.frame(oneDimNiche)
  colnames(oneDimNiche) <- c("NicheSize","2.5%","50%","97.5%")
  
  row.names(oneDimNiche) <- spp 
  
  print(oneDimNiche) 
}

#niches object #its all data of niches

#Calculate Bhattacharrya Coefficient for Niche sizes distrib
Bhatta.NicheFIN<-NULL
system.time({
  for (i in 1:nspp) {
    Bhatta.nicheini<-NULL
    for (j in 2:ncol(niches)) {
      bhatcoefnic<-dispRity::bhatt.coeff(niches[,i], niches[,j])
      
      Bhatta.nicheini<-cbind(Bhatta.nicheini,bhatcoefnic) 
    }
    dimnames(Bhatta.nicheini)<-list(spp[i],spp[2:5])
    Bhatta.NicheFIN<-rbind(Bhatta.NicheFIN, Bhatta.nicheini)
  }
  Bhatta.NicheFIN
})


save(niches, file = paste0("niches_", ndim, "_100k.RData"))


#######################################################
####Prepare object for Layman metrics code
output <- NULL
for (i in 1:nspp) {
  data <- BCiso.par[[i]]
  mu <- as.data.frame(data[1])
  output <- abind(output, mu, along = 3)
  
}

output2 <- NULL
for (i in 1:nisotope) {
  output1.1 <- output[,i,]
  output2 <- abind(output2, output1.1, along = 3)
}
dim(output2)

# mean mu y sigma for plots 3D MAY FAIL IN OTHER DIMENSIONS
output3 <- NULL

for (i in 1:nisotope) {
  promedio <- colMeans(output2[,,i])
  output3 <- rbind(output3, promedio)
}
meanmu <- t(output3)

output4 <- NULL
for (i in 1:nspp) {
  data <- BCiso.par[[i]]
  meansigma <- apply(data$Sigma, c(1,2), mean)
  output4 <- abind(output4, meansigma, along = 3)
}

#######################################################
#load from NicheRover output from Loader code 
matrix.data<-output2 #all iterated values

# Resultados ---------------------------------------------------------------

#Centroid locations modes and quantiles

## notice that all dimensions independent calculations like iso locations, 
#should be the same at the max number of dimensions compared to less number of dimensions
system.time({
  iso.location<-array(NA, dim = c(ncol(matrix.data),7,dim(matrix.data)[3]))#for this number of columns =4
  for (j in 1:nisotope) {
    isotopo<-matrix.data[,,j]
    for (i in 1:nspp) {
      mode1 <- mean(isotopo[,i])
      quantile1 <- quantile(isotopo[,i],probs=c(0.025,0.5,0.975)) 
      maximo <- max(isotopo[,i])
      minimo <- min(isotopo[,i])
      rango <- max(isotopo[,i]) - min(isotopo[,i])
      modaquantile <- c(mode1,quantile1, maximo, minimo, rango)
      iso.location[i,,j]<-modaquantile
    } 
  }
  dimnames(iso.location)<-list(spp,
                               c("media","2.5%","50%","97.5%", "max", "min", "rango"), isotopes)
  print(iso.location)
})

##########################################################################################
##########################################################################################

#Calculate Bhattacharrya Coefficient for Isotopes locations (mu)
Bhatta.ResFIN<-NULL
for (z in 1:nisotope) {
  Bhatta.Res<-NULL
  for (i in 1:nspp) {
    Bhatta.ini<-NULL
    for (j in 2:nspp) {
      bhatcoef<-dispRity::bhatt.coeff(matrix.data[,i,z], matrix.data[,j,z])
      Bhatta.ini<-cbind(Bhatta.ini,bhatcoef) 
    }
    dimnames(Bhatta.ini)<-list(spp[i],spp[2:5])
    Bhatta.Res<-rbind(Bhatta.Res, Bhatta.ini)
  }
  Bhatta.ResFIN<-abind(Bhatta.ResFIN,Bhatta.Res,along = 3)
}
#Bhatta.ResFIN

################################################################################
################################################################################
################################################################################
#Create function to calculate distance probability test of centroids
if (nisotope == 1) {
  distance<-function(a,b){
    l<-abs(a[1]-b[1])
    return(l)
  }
}else{
  if (nisotope == 2) {
    distance<-function(a,b){
      l<-sqrt((a[,1]-b[,1])^2+(a[,2]-b[,2])^2)
      return(l)
    }
  }else{
    if (nisotope == 3){
      distance<-function(a,b){
        l<-sqrt((a[,1]-b[,1])^2+(a[,2]-b[,2])^2+(a[,3]-b[,3])^2)
        return(l)
      }
    }else{
      if (nisotope == 5) {
        distance<-function(a,b){
          l<-sqrt((a[,1]-b[,1])^2+(a[,2]-b[,2])^2+(a[,3]-b[,3])^2+(a[,4]-b[,4])^2+(a[,5]-b[,5])^2)
          return(l)
        }
      }
    }
  } 
}


################################################################################
#from rossman

# Distance between centroids 
# to test splits the posterior distributions output
# into a null and a test arrays

system.time({
  if (ndim == 1) { # creating a different for an array with only 1 dimension (couldnt fix the code to run without the if else)
    null_mu<-array(output2[seq(1,nrow(output2),by=2),,], #the data from post dists
                   dim=c(nrow(output2)/2, #dimension for the rows of data
                         nspp#, #dimension for the colums of spp
                         #nisotope
                   )) # dimension for number of isotopes (1 or more sets of 2 entry (species cols x samples rows) matrices)
    test_mu<-array(output2[seq(2,nrow(output2),by=2),,],
                   dim=c(nrow(output2)/2,
                         nspp#,
                         #nisotope
                   ))
    
    dist.diff1<-array(NA, 
                      dim=c(nspp,
                            nspp,
                            nrow(output2)/2))
    
    for (j in 1:nspp){
      for (m in 1:nspp){
        dist.diff1[j,m,]<-distance(test_mu[,j],test_mu[,m]) - distance(test_mu[,j],null_mu[,j]) - distance(test_mu[,m],null_mu[,m])
      }}
    
    dist.diff2<-array(NA, dim=c(nspp,nspp,nrow(output2)/2))
    
    for (j in 1:nspp){
      for (m in 1:nspp){
        dist.diff2[j,m,]<-distance(test_mu[,j],test_mu[,m])
      }}
    
    dist.sum<-rep(list(matrix(NA,ncol=nspp,nrow=nspp)),3)
    
    for (i in 1:3){
      for (j in 1:nspp){
        for (m in 1:nspp){
          dist.sum[[i]][j,m]<-quantile(dist.diff2[j,m,],probs=c(0.025,0.5,0.975)[i])
        }}}
    
    names(dist.sum)<-c("q2.5","q50","q97.5")
    for(i in 1:3){
      colnames(dist.sum[[i]])<-spp
      rownames(dist.sum[[i]])<-spp
    }
    pdist<-matrix(NA, nrow=nspp,ncol=nspp)
    for(j in 1:nspp){
      for(m in 1:nspp){
        pdist[j,m]<-sum(dist.diff1[j,m,]>0)/(nrow(output2)/2)
      }
    }
    colnames(pdist)<-row.names(pdist)<-spp
  }else{
    null_mu<-array(output2[seq(1,nrow(output2),by=2),,], #the data from post dists
                   dim=c(nrow(output2)/2, #dimension for the rows of data
                         nspp, #dimension for the colums of spp
                         nisotope)) # dimension for number of isotopes (1 or more sets of 2 entry (species cols x samples rows) matrices)
    test_mu<-array(output2[seq(2,nrow(output2),by=2),,],
                   dim=c(nrow(output2)/2,
                         nspp,
                         nisotope))
    
    dist.diff1<-array(NA, 
                      dim=c(nspp,
                            nspp,
                            nrow(output2)/2))
    
    for (j in 1:nspp){
      for (m in 1:nspp){
        dist.diff1[j,m,]<-distance(test_mu[,j,],test_mu[,m,]) - distance(test_mu[,j,],null_mu[,j,]) - distance(test_mu[,m,],null_mu[,m,])
      }}
    
    dist.diff2<-array(NA, dim=c(nspp,nspp,nrow(output2)/2))
    
    for (j in 1:nspp){
      for (m in 1:nspp){
        dist.diff2[j,m,]<-distance(test_mu[,j,],test_mu[,m,])
      }}
    
    dist.sum<-rep(list(matrix(NA,ncol=nspp,nrow=nspp)),3)
    
    for (i in 1:3){
      for (j in 1:nspp){
        for (m in 1:nspp){
          dist.sum[[i]][j,m]<-quantile(dist.diff2[j,m,],probs=c(0.025,0.5,0.975)[i])
        }}}
    
    names(dist.sum)<-c("q2.5","q50","q97.5")
    for(i in 1:3){
      colnames(dist.sum[[i]])<-spp
      rownames(dist.sum[[i]])<-spp
    }
    pdist<-matrix(NA, nrow=nspp,ncol=nspp)
    for(j in 1:nspp){
      for(m in 1:nspp){
        pdist[j,m]<-sum(dist.diff1[j,m,]>0)/(nrow(output2)/2)
      }
    }
    colnames(pdist)<-row.names(pdist)<-spp
  }
  
  
})

#print("probability of having different centroids")
#print(pdist, digits=3)

#print("Distance between centroids")
#print(dist.sum, digits=3)

################################################################################
################################################################################
################################################################################

#1 Calculate universe centroid
#objetivo 1

#isotope <- colnames(dataCNS[,2:ncol(dataCNS)])
isotope<-isotopes

outputCD <- NULL

for (i in 1:nisotope) {
  meanballs <- mean(matrix.data[,,i])
  cuantiles <- quantile(matrix.data[,,i], probs=c(0.025,0.5,0.975))
  outputCD <- rbind(outputCD, c(meanballs, cuantiles))
}

dimnames(outputCD) <- list(isotopes, c("mean", "2.5%", "50%", "97.5%"))

#print("Universe Centroid")
#print(outputCD)


#2.1 Calculate distance of each point to Universe centroid
#objetivo 2
system.time({
  
  output5 <- NULL
  
  
  for (j in ncol) {
    for (z in nrows) {
      if (nisotope == 1) {
        centroid <- abs(matrix.data[z,j,1] - outputCD[1])
      }else{
        if (nisotope == 2) {
          centroid <- sqrt((matrix.data[z,j,1] - outputCD[1])^2 + (matrix.data[z,j,2] - outputCD[2])^2)
        }else{
          if (nisotope == 3){
            centroid <- sqrt((matrix.data[z,j,1] - outputCD[1])^2 + (matrix.data[z,j,2] - outputCD[2])^2 + 
                               (matrix.data[z,j,3] - outputCD[3])^2)
          }else{
            if (nisotope == 5) {
              centroid <- sqrt((matrix.data[z,j,1] - outputCD[1])^2 + (matrix.data[z,j,2] - outputCD[2])^2 + 
                                 (matrix.data[z,j,3] - outputCD[3])^2 + (matrix.data[z,j,4] - outputCD[4])^2 +
                                 (matrix.data[z,j,5] - outputCD[5])^2)
            }
          }
        }
      }
      
      
      output5 <- c(output5, centroid)
    }
  }
  
  #output5# are all distances of each point/iteration to the universe Centroid
  
  matrixfinal <- matrix(output5, nrow = nrow(matrix.data[,,1]), ncol = length(ncol))
  quantmatrix <- c(mean(matrixfinal), quantile(matrixfinal, probs = c(0.025, 0.5, 0.975)))
  names(quantmatrix) <- c("mean", "2.5%", "50%", "97.5%")
})
##save for each DIMENSION and ITERATION type
save(matrixfinal, file = paste0("output5_", nisotope, "_100k.RData")) #this was carrying an error name saving "matrix" which didnt exist

#print("Mean overall distance to Centroid to Universe")
#print(quantmatrix) #this is just a matrix of distances


#2.2 Obtain quantile results for Distance to centroid
# objetivo 2.2

output6 <- NULL

for (i in ncol) {
  percentil <- quantile(matrixfinal[,i], probs = c(0.025, 0.5, 0.975))
  promedio <- mean(matrixfinal[,i])
  
  output6 <- rbind(output6, c(promedio, percentil))
}

dimnames(output6) <- list(spp, c("mean", "2.5%", "50%", "97.5%"))

#print("Distance to Universe centroid per group")
#print(output6) 



##NND for all community new
print(paste("Started NND calculations at", Sys.time()))
system.time({
  
  output12 <- NULL
  
  for (z in ncol) {
    output11 <- NULL
    for (i in nrows) {
      if (nisotope == 1) {
        xspp <- abs(matrix.data[i,z,1] - matrix.data[,,1])
      }else{
        if (nisotope == 2) {
          xspp <- sqrt((matrix.data[i,z,1] - matrix.data[,,1])^2 + (matrix.data[i,z,2] - matrix.data[,,2])^2)
        }else{
          if (nisotope == 3){
            xspp <- sqrt((matrix.data[i,z,1] - matrix.data[,,1])^2 + (matrix.data[i,z,2] - matrix.data[,,2])^2 + 
                           (matrix.data[i,z,3] - matrix.data[,,3])^2)
          }else{
            if (nisotope == 5) {
              xspp <- sqrt((matrix.data[i,z,1] - matrix.data[,,1])^2 + (matrix.data[i,z,2] - matrix.data[,,2])^2 + 
                             (matrix.data[i,z,3] - matrix.data[,,3])^2 + (matrix.data[i,z,4] - matrix.data[,,4])^2
                           + (matrix.data[i,z,5] - matrix.data[,,5])^2)
            }
          }
        }
      }
      xspp[xspp == 0] <- NA
      
      distmin <- min(xspp, na.rm = TRUE)
      output11 <- rbind(output11, distmin)
      
      
    }
    
    output12 <- cbind(output12, output11)
    
  }
  
  meanoutput12 <- mean(output12)
  sdoutput12 <- sd(output12)
  quantileoutput12 <- quantile(output12,  probs = c(0.025, 0.5, 0.975))
  resumenoutput12 <- c(meanoutput12, sdoutput12, quantileoutput12)
  
  names(resumenoutput12) <- c("mean", "sd", "2.5%", "50%", "97.5")
})

##save for each DIMENSION and ITERATION type
save(output12, file = paste0("output12_", ndim, "_100k.RData"))

#old version if save
####save(output12, file = "moutput12_2D_100k.RData") ###BEWARE OF OVERWRITING FILES FOR DIFFERENTE DIMENSIONS

#print("Mean NND to universe")
#print(resumenoutput12)

## NEW NEW calculate SDNND
samplessd <- NULL
for (i in 1:nsamples) {
  samples <- sd(sample(output12, size = 10))
  
  samplessd <- c(samplessd, samples)
}

resumesample <- c(mean(samplessd), quantile(samplessd,  probs = c(0.025, 0.5, 0.975)))
names(resumesample) <- c("mean", "2.5%", "50%", "97.5")

save(samplessd, file = paste0("samplessd_", ndim, "_100k.RData"))



#NEW NEW isotope range 
matrixdiffsample <- NULL
matrixrange<-NULL
for (j in 1:ndim) {
  matrixslice <- matrix.data[,,j]
  diffsample <- NULL
  for (i in 1:nsamples) {
    sampleslice <- sample(matrixslice, 10)
    diff <- max(sampleslice) - min(sampleslice)
    
    diffsample <- rbind(diffsample, diff)
  }
  matrixrange<-abind(matrixrange, diffsample, along=3) 
  
  meandiffsample <- mean(diffsample)
  quantilesample <- quantile(diffsample,  probs = c(0.025, 0.5, 0.975))
  
  resultsample <- c(meandiffsample, quantilesample)
  names(resultsample) <- c("mean", "2.5%", "50%", "97.5")
  
  matrixdiffsample <- rbind(matrixdiffsample, resultsample)
}

rownames(matrixdiffsample) <- isotopes

save(matrixrange, file = paste0("matrixrange_", ndim, "_100k.RData"))

########################################
########################################
########################################
#all outputs at once
print("Centroid locations modes and quantiles per group")
print(iso.location)

print("Probability Test for differences between distributions of Centroid locations")
print(Bhatta.ResFIN)

print("probability of having different centroids")
print(pdist, digits=3)

print("Distance between centroids")
print(dist.sum, digits=3)

print("Universe Centroid")
print(outputCD)

print("Mean overall distance to Centroid to Universe")
print(quantmatrix)

print("Distance to Universe centroid per group")
print(output6) 

#print("Mean NND within species")
#print(output4) 

#print("Mean NND of per species to other species")
#print(output8) #each species, separately, vs other species + quantiles

print("Mean NND to universe")
print(resumenoutput12)

print("SDNND to universe")
print(resumesample)

print("Dimension Range")
print(matrixdiffsample)

print(paste("Finished at", Sys.time()))
