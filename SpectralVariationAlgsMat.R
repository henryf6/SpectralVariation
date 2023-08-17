## ******************************************************************** ##
## SpectralVariationAlgs.R 
##
## Author: Henry Frye
## Date Created: Sep 2022
##
##
## Purpose:
## Create functions for spectral variation analyses using a more 
## generalized matrix approach
## ******************************************************************** ##

## ******************************************************************** ##
####Set Up####
## ******************************************************************** ##

library(tidyverse)
library(foreach)
#note that the inter and intra spectral angle variability can be a bit slow to run
# since they rely on double loops and I couldn't figure out how to quite get them to run
# on apply functions.


#convert radians to degrees
rad2deg <- function(rad) {(rad * 180) / (pi)}

#Calculate AVW with variable wavelength window
AVWCalc <-  function(spectra, wavelengths){#spectra dimensions and wavelengths must match
  
  AVW <- vector(length= dim(spectra)[1])
  Numer <- vector(length = dim(spectra)[2])
  Denom <- vector(length = dim(spectra)[2])
  
  for(j in 1:dim(spectra)[1]) {
    for(i in 1:length(wavelengths)) {
      Numer[i] <- spectra[j,i]
      Denom[i] <- spectra[j,i] / wavelengths[i]
    }
    
    AVW[j] <-  sum(sapply(Numer, sum, na.rm = TRUE)) /  sum(sapply(Denom, sum, na.rm = TRUE))
  }
  return(AVW)
}


#Spectral Angle calculation for matrices, cross-checked against hsdar's function
SpectralAngle <- function(test, reference, margin = 1) {
  if(identical(test, reference) == TRUE) { 
    return(0)
  } else {
    top = apply(test*reference, margin, sum, na.rm = TRUE)
    b1 = sqrt(apply(test^2, margin, sum, na.rm = TRUE))
    b2 = sqrt(apply(reference^2, margin, sum, na.rm = TRUE))
    return(acos(top / (b1*b2)))
  }
}


#Intra-specific variability that requires only a spectral matrix, rows are samples and
# columns are wavelenghts, note that function only returns the unique pairwise distances
# hence the upper triangle at the end.
####in this version I tried to merge the hsdar and matrix options, but ultimately left them to different functions


IntraSpecVarTest <- function(spectra) {
  #if(class(spectra) == "Speclib") {
  # empty <- matrix(nrow = length(spectra@ID), ncol = length(spectra@ID))
  
  #  for(i in 1:length(spectra@ID)){
  #   for(j in 1:length(spectra@ID)){
  #    empty[i,j] <- sam(spectra[i,],spectra[j,])
  #  }
  #}
  
  #} else {
  empty <- matrix(nrow = nrow(spectra), ncol = nrow(spectra))
  
  for(i in 1:nrow(spectra)){
    for(j in 1:nrow(spectra)){
      empty[i,j] <- SpectralAngle(spectra[i,],spectra[j,])
      #myfun_v <- Vectorize(myfun) 
    }
  }
  #}
  values <- empty[lower.tri(empty)]
  ValueDeg <- rad2deg(values)
  return(ValueDeg)
  
}


#Inter variability between groups of spectral angle. This calculates both inter and intra
# while removing the duplicate values in the intra-group case.

InterSpecVarTest <- function(spectra1, spectra2) {
  
  if(identical(spectra1,spectra2) == TRUE){
    
    empty <- matrix(nrow = nrow(spectra1), ncol = nrow(spectra2))
    
    for(i in 1:nrow(spectra1)){
      for(j in 1:nrow(spectra2)){
        empty[i,j] <- SpectralAngle(spectra1[i,],spectra2[j,]) }}
    
    values <- empty[lower.tri(empty)]
    ValueDeg <- rad2deg(values)
    return(ValueDeg)
    
  } else {
    empty <- matrix(nrow = nrow(spectra1), ncol = nrow(spectra2))
    
    for(i in 1:nrow(spectra1)){
      for(j in 1:nrow(spectra2)){
        empty[i,j] <- SpectralAngle(spectra1[i,],spectra2[j,]) }}
    
  }
  ValueDeg <- rad2deg(as.vector(empty))
  return(ValueDeg)
}



#This is the same version of the function above, but now with a parallel computing option
#Inter variability between groups of spectral angle. This calculates both inter and intra
# while removing the duplicate values in the intra-group case.

library(doParallel)
cores=detectCores()
cl <- makeCluster(cores[1]-2) #not to overload your computer
registerDoParallel(cl)


InterSpecVarParallel <- function(spectra1, spectra2) {
  
  if(identical(spectra1,spectra2) == TRUE){
    
    # empty <- matrix(nrow = nrow(spectra1), ncol = nrow(spectra2))
    # 
    # foreach(i = nrow(spectra1), .combine = 'cbind') %:% 
    #   foreach(j = nrow(spectra2), .combine = 'c') %dopar% {
    #     SpectralAngle(i,j) }
    
    emptyeach <- foreach(i = 1:nrow(spectra1), .combine = "rbind",  .export = "SpectralAngle") %:% 
      foreach(j = 1:nrow(spectra2), .combine = "c",  .export = "SpectralAngle") %dopar% {
        SpectralAngle(spectra1[i,],spectra2[j,])
      }
    
    
    
    values <- emptyeach[lower.tri(emptyeach)]
    ValueDeg <- rad2deg(values)
    return(ValueDeg)
    
  } else {
    # empty <- matrix(nrow = nrow(spectra1), ncol = nrow(spectra2))
    # 
    # foreach(i = nrow(spectra1), .combine = 'cbind') %:% 
    #   foreach(j = nrow(spectra2), .combine = 'c') %dopar% {
    #     SpectralAngle(i,j) }
    # 
    emptyeach <- foreach(i = 1:nrow(spectra1), .combine = "rbind", .export = "SpectralAngle") %:% 
      foreach(j = 1:nrow(spectra2), .combine = "c",  .export = "SpectralAngle") %dopar% {
        SpectralAngle(spectra1[i,],spectra2[j,])
      }
    
  }
  ValueDeg <- rad2deg(as.vector(emptyeach))
  return(ValueDeg)
}
