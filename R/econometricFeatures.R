#' Moran 2D spatial analysison DICOM images (extension of moddicomV2 package)
#'
#' @description   Package of functions to extract and test significance level of slice-wise spatial features (Moran analysis) from DICOM images.
#'                \itemize{
#'                \item \code{whatDoYouWantBack=c("wMatrix")} it gives back: 1) wStar: the normalized W matrix for the current slice; 2) Coords: x and y coordinates of the pixels which contributed to W* (if dimension of W* is NxN, then length of Coords is N); 3) VcEroso: the image (voxelCube) after the erosion.
#'                \item \code{whatDoYouWantBack=c("moranGrayMean") } it gives back: 1) I_coeff: the I value computed on each slice; 2) I_nullHp: the expected I coefficient for null hp (no spatial autocorr.); 3) varI_nullHp: the expected variance of I coefficient for null hp.
#'                \item \code{whatDoYouWantBack=c("extractedI") } it gives back: 1) nCoords: number of contributing pixels; 2) Icoeff: value of I coefficient for the slice; 3) pValue: statistical significance of the I coefficient.
#'                \item \code{whatDoYouWantBack=c("features")  } it gives back: 1) Icoeff: vector of per-slice I coefficients; 2) meanI: mean of per-slice I values: 3) medianI: median of per-slice I values; 4) maxI: maximum of per-slice I values; 5) minI: minimum of per-slice I values; 6) rangeI:  maxI - minI; 7) maxISq: squared maximum of per-slice I values; 8) minISq: squared minimum of per-slice I values; 9) weightedMeanI: weighted mean of per-slice I values; weights are number of contributing pixels.
#'                \item \code{whatDoYouWantBack=c("plotFeatures")} plots the absolute frequency distribution for each feature
#'                \item \code{whatDoYouWantBack=c("dataframeI")} it gives back a dataframe with the extracted features as columns (a record for each patient)
#'
#'                }
#' @param Parameters for calculateCluster methoda are:
#'   \itemize{
#'    \item \code{d* } distance threshold to biuld the W matrix
#'    \item \code{erosionMarginX } number of pixel to "discard" at the border along the x direction
#'    \item \code{erosionMarginY } number of pixel to "discard" at the border along the y direction
#'    \item \code{pixelSpacingX } pixel spacing in the x direction (default =1)
#'    \item \code{pixelSpacingY } pixel spacing in the y direction (default =1)
#'    \item \code{whatDoYouWantBack } a character vector including one or more of the following values:"wMatrix","moranGrayMean","extractedI","features","plotFeatures","dataframeI".
#'   }
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
#' # calculate the wMatrix for each of them (distance threshold = 2, erosion along x and y =2, pixel spacing along x and y =1)
#' wMatrix <-F01.pipeline.01( inputGTV = GTV, 2,2,2,whatDoYouWantBack=c("wMatrix"))
#'
#' }
#'
#' # calculate slice-wise I coefficients and null Hypotesis values for each of them (distance threshold = 2, erosion along x and y =2, pixel spacing along x and y =1)
#' m.g.mean <-F01.pipeline.01( inputGTV = GTV, 2,2,2,whatDoYouWantBack=c("moranGrayMean"))
#'
#' }
#'
#'
#' # calculate slice-wise I coefficients and Moran test values for each of them (distance threshold = 2, erosion along x and y =2, pixel spacing along x and y =1)
#' stats.I <-F01.pipeline.01( inputGTV = GTV, 2,2,2,whatDoYouWantBack=c("extractedI"))
#'
#' }
#'
#'
#' # calculate features from significant I coefficients (distance threshold = 2, erosion along x and y =2, pixel spacing along x and y =1)
#' features.I <-F01.pipeline.01( inputGTV = GTV, 2,2,2,whatDoYouWantBack=c("features"))
#'
#' }
#'
#'
#' # calculate features from significant I coefficients and put them in a dataframe
#' dataframe.I <-F01.pipeline.01( inputGTV = GTV, 2,2,2,whatDoYouWantBack=c("dataframeI"))
#'
#' }


F01.pipeline.01<-function(inputGTV,dStar, erosionMarginX, erosionMarginY, pixelSpacingX=1, pixelSpacingY=1, whatDoYouWantBack=c("statsI")) {

  # operazioni preliminari
  objS<-services();
  stats.I<-list()
  final.result<-list()
  inputGTV.eroded <-objS$applyErosion(inputGTV,erosionMarginX,erosionMarginY, margin.z=0)

  for(possibiliResults in whatDoYouWantBack) final.result[[possibiliResults]]<-list();

  # cicla su tutti i pazienti
  for (k in seq(1,length(inputGTV))){
    cat(c("\n Patient: ",k ) )
    patID<-names(inputGTV)[[k]]
    vcNonEroso <- inputGTV[[k]]$masked.images$voxelCube
    vcEroso <- inputGTV.eroded[[k]]$masked.images$voxelCube
    listaVoxelErosi<-which(!is.na(vcEroso),arr.ind = TRUE)
    # per prima cosa calcola la matrice W
    if (length(listaVoxelErosi) != 0){
      wMatrix<-calcola.matrice.W(voxelCube,dStar, erosionMarginX, erosionMarginY, vcNonEroso, vcEroso, listaVoxelErosi, pixelSpacingX, pixelSpacingY)
      if("wMatrix" %in% whatDoYouWantBack) final.result$wMatrix[[patID]]<-wMatrix

     # ora calcola il moranGrayMean
      mGM<-calcola.moranGrayMean(wMatrix)
      if("moranGrayMean" %in% whatDoYouWantBack) final.result$moranGrayMean[[patID]]<-mGM

     # ora calcola I
      stats.I<- calcola.statsT(mGM)
      if("statsI" %in% whatDoYouWantBack) final.result$statsI[[patID]]<-stats.I

      # extract denoised (more than 20 contributing pixels) amd significant (pValue less than 0.05) moran Is
      extractedI <- extract.I(wMatrix,mGM,stats.I)
      if("extractedI" %in% whatDoYouWantBack) final.result$extractedI[[patID]]<-extractedI

      #
      features <- calcola.features(extractedI)
      if("features" %in% whatDoYouWantBack) final.result$features[[patID]] <- features
    }
  }

  if ("plotFeatures" %in% whatDoYouWantBack) plot.features(final.result$features)
  dataFrameI <- crea.dataframeI(final.result$features)
  if ("dataFrameI" %in% whatDoYouWantBack) final.result$dataFrameI <- dataFrameI

  return(final.result)
}


crea.dataframeI <- function (features){
  rowsId <-sub(".*/", "", names(features))
  nomi <- names(features[[1]])[2:length(features[[1]])]
  tmp <- list()
    for (pat in 1:length(rowsId)){
      tmp[[pat]] <- rbind(features[[pat]][2:length(nomi)+1])
    }
  data <- do.call(rbind.data.frame, tmp)
  data <- cbind(rowsId, data)
  data[2:length(data)] <- sapply(data[2:length(data)], as.numeric)
  return(data)
}


plot.features <- function(features){
  meanIall <- numeric()
  medianIall <- numeric()
  maxIall <- numeric()
  minIall <- numeric()
  rangeIall <- numeric()
  maxISqall <- numeric()
  minISqall <- numeric()
  weightedMeanIall <- numeric()

    for (z in 1:length(features)){
    meanIall[z] <- features[[z]]$meanI
    medianIall[z] <- features[[z]]$medianI
    maxIall[z] <- features[[z]]$maxI
    minIall[z] <- features[[z]]$minI
    rangeIall[z] <- features[[z]]$rangeI
    maxISqall[z] <- features[[z]]$maxISq
    minISqall[z] <- features[[z]]$minISq
    weightedMeanIall[z] <- features[[z]]$weightedMeanI
    }


  png("meanIall.png")
  hist(meanIall)
  dev.off()
  png("medianIall.png")
  hist(medianIall)
  dev.off()
  png("maxIall.png")
  hist(maxIall)
  dev.off()
  png("minIall.png")
  hist(minIall)
  dev.off()
  png("rangeIall.png")
  hist(rangeIall)
  dev.off()
  png("maxISqall.png")
  hist(maxISqall)
  dev.off()
  png("minISqall.png")
  hist(minISqall)
  dev.off()
  png("weightedMeanIall.png")
  hist(weightedMeanIall)
  dev.off()
  return()
}


calcola.features <- function(extractedI){
  Icoeff <- numeric()
  nPixel <- numeric()
  for (z in 1:length(extractedI)){
    Icoeff[z] <- extractedI[[z]]$Icoeff
    nPixel[z] <- extractedI[[z]]$nCoords
  }
  meanI <- mean(Icoeff,na.rm=T)
  medianI <- median(Icoeff,na.rm=T)
  maxI <- max(Icoeff,na.rm=T)
  minI <- min(Icoeff,na.rm=T)
  rangeI <- maxI - minI
  maxISq <- maxI^2
  minISq <- minI^2
  weightedMeanI <- sum(Icoeff * nPixel,na.rm=T)/sum(nPixel,na.rm=T)

  features <- list("Icoeff" = Icoeff, "meanI" = meanI, "medianI" = medianI, "maxI" = maxI, "minI" = minI, "rangeI" = rangeI, "maxISq" = maxISq, "minISq" = minISq, "weightedMeanI" = weightedMeanI)

  return(features)
}


extract.I <- function(Wlist,moranGray,testResults){
  moranGray <-  moranGray[!sapply(moranGray, is.null)]
  nCoords <- numeric()
  extractedI <- list()
    for (slice in 1:length(moranGray)){
      Npixel <- numeric()
      Npixel[slice] <- dim(Wlist[[slice]]$coords)[1]
       if (Npixel[slice] > 20 & testResults$pValues[slice] <= 0.05){
           nCoords[slice] <- Npixel[slice]
           extractedI[[slice]] <- list("sliceId"=slice, "nCoords" = nCoords[slice], "Icoeff"= moranGray[[slice]]$I_coeff, "pValue"= testResults$pValues[slice])
       }
    }
  return(extractedI)
}

calcola.statsT<-function(moranGray){
  Imean <- numeric()
  Icoeff <- numeric()
  Inull <- numeric()
  varNull <- numeric()
  zScores <- numeric()
  pValues <- numeric()
    for(z in seq(1:length(moranGray))){
      if (length(moranGray[[z]]) == 0) {Icoeff[z] <- NA}
      else {Icoeff[z] <- moranGray[[z]]$I_coeff
      Inull[z] <- moranGray[[z]]$I_nullHp
      varNull[z] <- moranGray[[z]]$varI_nullHp}
    }
  Icoeff <- sapply(Icoeff,unlist)
  Icoeff <- Icoeff[which(complete.cases(Icoeff))]
  Inull <- sapply(Inull,unlist)
  Inull <- Inull[which(complete.cases(Inull))]
  varNull <- sapply(varNull,unlist)
  varNull <- varNull[which(complete.cases(varNull))]
  slice <- seq(1:length(Icoeff))
  #plot(slice,Icoeff)
    for(z in seq(1:length(Icoeff))){
      if (length(Icoeff) == length(Inull)){
        zScores[z] <- (Icoeff[z] - Inull[z])/sqrt(varNull[z])
        pValues[z] <- 2*pnorm(-abs(zScores[z]))
      }
      else skip
    }
  testResults <- list("zScores"= zScores, "pValues"= pValues)
  #Imean <- mean(Icoeff)
  return(testResults)
}


calcola.moranGrayMean<-function(Wlist) {
  #ciclo sulle slice
  moranGray<-list()
    for(z in seq(1:length(Wlist))){
     if(length(Wlist[[z]]$wStar) < 400) next
      Wstar <- Wlist[[z]]$wStar
      coords <- Wlist[[z]]$coords
      vcEroso <- Wlist[[z]]$vcEroso
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
      #calcolo lo scalare di correlazione spaziale I
      laggedGray <- Wstar %*% grayVar
      laggedGray <- laggedGray[which(!is.na(laggedGray))]
      grayVar <- grayVar[which(!is.na(laggedGray))]
      I <- grayVar %*% laggedGray
      norm <- grayVar %*% grayVar
      Inorm <- as.numeric(I/norm)

      #### CALCOLO ANCHE IL VALOR MEDIO E LA VARIANZA DI I SOTTO L'IPOTESI NULLA
      if (is.nan(I) == T) next
      t_Wstar <- t(Wstar) # 'a trassposta
      N <- nrow(Wstar)
      expectI = 1 / ( - N + 1 );
      S1 <- 0
      S2 <- 0
      S4 <- 0
      S5 <- 0
      sumSquaresW <- 0
      varI <- 0
      lista.varMoranNull<-.C("varMoranNull",
                             as.double(Wstar),
                             as.double(t_Wstar),
                             as.integer(N),
                             as.double(S1),
                             as.double(S2));

      S1 <- lista.varMoranNull[[4]]
      S2 <- lista.varMoranNull[[5]]
      S3 <-  nrow(Wstar) * sum((grayVar^4)) * (1/sum(grayVar^2)^2)
      sumSquaresW <- sum(Wstar)^2
      S4 <- (N^2 - 3 * N + 3) * S1 - N * S2 + 3 * sumSquaresW
      S5 <- (N^2 - N) * S1 - 2 * N * S2 + 6 * sumSquaresW
      varI <- ((N * S4 - S3 * S5 ) / ((N - 1) * (N - 2) * (N - 3) * sumSquaresW)) - expectI^2

      moranGray[[z]] <- list("I_coeff"= Inorm, "I_nullHp"= expectI, "varI_nullHp"= varI)
    }
  return(moranGray)
}

calcola.matrice.W<-function(voxelCube,dStar, erosionMarginX, erosionMarginY, vcNonEroso, vcEroso, listaVoxelErosi, pixelSpacingX, pixelSpacingY) {
  Wlist<-list()
  #ora seleziono le coord (x,y) per ogni slice z e calcolo la matrice delle distanze euclidee tra le coppie di punti (funzione: rdist) e la matrice W definita nel ciclo for(i)for(j)
    for(z in seq(1,max(listaVoxelErosi[,3]))){
      cat(c("."))
      slice <- listaVoxelErosi[which(listaVoxelErosi[,3]==z),]
      if (class(slice) == "integer" || length(slice) == 0) next
      slice2D <- slice[,1:2]
      slice2D[,1] <- slice2D[,1]*pixelSpacingX
      slice2D[,2] <- slice2D[,2]*pixelSpacingY
      D <- matrix()
      D <- rdist(slice2D)
      W <- matrix(0,nrow=nrow(D), ncol=ncol(D))
      Wstar <- matrix(0,nrow=nrow(W), ncol=ncol(W))
      # invoca la funzione in C forzando esplicitamente il tipo
      # dei parametri (importante)
      lista.risultato<-.C("internalLoop",
                          as.double(D),
                          as.integer(nrow(D)),
                          as.double(dStar),
                          as.matrix(W));
      rm(D)
      rm(W)
      aaa <- lista.risultato[[4]]
      rm(lista.risultato)
      #normalizzo i coefficienti di W
      Wnorm <- numeric()
      Wnorm <- rowSums(aaa)
        for (i in seq(1,nrow(aaa))){
          if (Wnorm[i]!=0){
          Wstar[i,] <-aaa[i,]/Wnorm[i]
          }
        }
      rm(aaa)
      #costruisco le k*z(k) (nPazienti*nslicePaziente) matrici W
      Wstar <- round(Wstar,2)
      Wlist[[z]] <- list("wStar"= Wstar,"coords"= slice, "vcEroso"= vcEroso[,,z])
      rm(Wstar)
      rm(slice)
      gc()
    }
  return(Wlist)
}

