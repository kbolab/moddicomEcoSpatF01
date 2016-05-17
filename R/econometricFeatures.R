#' cdscdscd
#' 
#' @description  dfdsfdfdsfdsfds
#' @export
#' @useDynLib moddicomEcoSpatF01 
#' @import fields
#' @examples \dontrun{
#' 
#' # create an mmButo object and load some DICOM series 
#' 
#' obj<-mmButo()
#' obj$openDICOMFolder(obj = obj, pathToOpen='./DICOMSeries/study4Radiomics' );
#' 
#' # get a ROI Voxel list
#' GTV<-obj$getROIVoxel(ROIName = "GTV")
#' 
#' # and calculate the wMatrix for each of them
#' wMatrix<-F01.wMatrix( inputGTV = GTV, 2,2,2)
#' 
#' } 
F01.wMatrix <- function(inputGTV,dStar, erosionMarginX, erosionMarginY){
  objS<-services();
  Wlist <- list()
  inputGTV.eroded <-objS$applyErosion(inputGTV,erosionMarginX,erosionMarginY, margin.z=0)
  #create "eroded" mmButoStructureVoxelList object
  #loop on patients
  
  for (k in seq(1,length(inputGTV))){
    cat(c("\n Patient: ",k ) )
    patID<-names(inputGTV)[[k]] 
    Wlist[[ patID ]] <-list()
    vcNonEroso <- inputGTV[[k]]$masked.images$voxelCube
    vcEroso <- inputGTV.eroded[[k]]$masked.images$voxelCube
    listaVoxelErosi<-which(!is.na(vcEroso),arr.ind = TRUE)
    
    #ora seleziono le coord (x,y) per ogni slice z e calcolo la matrice delle distanze euclidee tra le coppie di punti (funzione: rdist) e la matrice W definita nel ciclo for(i)for(j)
    for(z in seq(1,max(listaVoxelErosi[,3]))){
      cat(c("."))
      slice <- listaVoxelErosi[which(listaVoxelErosi[,3]==z),]
      if (class(slice) == "integer" || length(slice) == 0) next
      slice2D <- slice[,1:2]
      D <- matrix()
      D <- rdist(slice2D)
      W <- matrix(0,nrow=nrow(D), ncol=ncol(D))
      Wstar <- matrix(0,nrow=nrow(W), ncol=ncol(W))
      # invoca la funzione in C forzando esplicitamente il tipo
      # dei parametri (importante)
      lista.risultato<-.C("internalLoop",
                          as.matrix(D),
                          as.integer(nrow(D)),
                          as.double(dStar),
                          as.matrix(W));
      #normalizzo i coefficienti di W
      for (i in seq(1,nrow(lista.risultato[[4]]))){
        Wnorm <- 0
        for (j in seq(1,ncol(lista.risultato[[4]]))){ 
          Wnorm <- Wnorm + lista.risultato[[4]][i,j]
        }
        Wstar[i,] <-lista.risultato[[4]][i,]/Wnorm
      }
      #costruisco le k*z(k) (nPazienti*nslicePaziente) matrici W
      Wstar <- round(Wstar,2)
      Wlist[[ patID ]][[z]] <- list("wStar"= Wstar,"coords"= slice, "vcEroso"=vcEroso[,,z])
    }
    #chiuso il ciclo per le slice del paziente kesimo
  }
  #chiuso il ciclo sui pazienti
  return(Wlist)
}

#' cdscdg45g54scd
#' 
#' @description  dfdsg54g54g54g3fdfdsfdsfds
#' @export
F01.moranGrayMean <- function(Wlist){
  moranGray <- list()
  #ciclo sui pazienti
  for(k in seq(1:length(Wlist))){
    patID<-names(Wlist)[k]
    moranGray[[ patID ]] <- list()
    #ciclo sulle slice
    for(z in seq(1:length(Wlist[[k]]))){
      wStar <- Wlist[[k]][[z]]$wStar
      coords <- Wlist[[k]][[z]]$coords
      vcEroso <- Wlist[[k]][[z]]$vcEroso
      x <- coords[,1]
      y <- coords[,2]
      grayLevels <- numeric()
      #il minimo mi serve per portare a zero il valore minimo se negativo
      minumGray <- min(vcEroso,na.rm=T)
      for(i in seq(1:length(x))){
        if (minumGray < 0 ) {grayLevels[i] <- vcEroso[x[i],y[i]] - minumGray}
        else {grayLevels[i] <- vcEroso[x[i],y[i]]}
      }
      #hist(grayLevels)
      grayMean <- mean(grayLevels)
      # definisco la grandezza da "laggare" come differenza rispetto al valor medio
      grayVar <- grayLevels - grayMean
      if (length(wStar) == 0) next
      #calcolo lo scalare di correlazione spaziale I
      laggedGray <- wStar %*% grayVar
      I <- grayVar %*% laggedGray
      norm <- grayVar %*% grayVar
      Inorm <- as.numeric(I/norm)
      moranGray[[ patID ]][[z]] <- list("I_coeff"= Inorm)
    }
  }
  return(moranGray)
}

#' cdscdgg5445g54scd
#' 
#' @description  dfdsg5445g5g54g54g3fdfdsfdsfds
#' @export
F01.statsI <- function(moranGray){
  Imean <- numeric()
  for(k in seq(1:length(moranGray))){
    Icoeff <- numeric()
    for(z in seq(1:length(moranGray[[k]]))){
      if (length(moranGray[[k]][[z]]) == 0) {Icoeff[z] <- NA}
      else {Icoeff[z] <- moranGray[[k]][[z]]}
    }
    Icoeff <- sapply(Icoeff,unlist)
    Icoeff <- Icoeff[which(complete.cases(Icoeff))]
    slice <- seq(1:length(Icoeff))
    plot(slice,Icoeff)
    Imean[k] <- mean(Icoeff)
  }
  return(Imean)
}
