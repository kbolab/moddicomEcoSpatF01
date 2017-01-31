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


F01.pipeline.01_yAxis_xFixed<-function(inputGTV, dStar, erosionMarginX=0, erosionMarginY=0, whatDoYouWantBack=c("extractedI","dataFrameI","features")) {

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

    # se x è fisso (piano x-y) --> scambio x con z e z con x
    pixelSpacingY <- inputGTV[[k]]$geometricalInformationOfImages$pixelSpacing[2]
    pixelSpacingX <- as.numeric(inputGTV[[k]]$geometricalInformationOfImages$SliceThickness)
    # per prima cosa calcola la matrice W
    if (length(listaVoxelErosi) != 0){
      wMatrix<-calcola.matrice.Wy_xFixed(voxelCube,dStar, erosionMarginX, erosionMarginY, vcNonEroso, vcEroso, listaVoxelErosi, pixelSpacingX, pixelSpacingY)
      if("wMatrix" %in% whatDoYouWantBack) final.result$wMatrix[[patID]]<-wMatrix

      # ora calcola il moranGrayMean
      mGM<-calcola.moranGrayMeany_xFixed(wMatrix)
      if("moranGrayMean" %in% whatDoYouWantBack) final.result$moranGrayMean[[patID]]<-mGM

      #ora calcola I
      stats.I<- calcola.statsTy_xFixed(mGM)
      if("statsI" %in% whatDoYouWantBack) final.result$statsI[[patID]]<-stats.I

      #extract denoised (more than 20 contributing pixels) amd significant (pValue less than 0.05) moran Is
      extractedI <- extract.Iy_xFixed(wMatrix,mGM,stats.I)
      if("extractedI" %in% whatDoYouWantBack) final.result$extractedI[[patID]]<-extractedI

      #
      features <- calcola.featuresy_xFixed(extractedI)
      if("features" %in% whatDoYouWantBack) final.result$features[[patID]] <- features
    }
  }

  # if ("plotFeatures" %in% whatDoYouWantBack) plot.features(final.result$features)

  dataFrameI <- crea.dataframeIy_xFixed(final.result$features)
  if ("dataFrameI" %in% whatDoYouWantBack) final.result$dataFrameI <- dataFrameI

  return(final.result)
}



crea.dataframeIy_xFixed <- function (features){
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


calcola.featuresy_xFixed <- function(extractedI){
  Icoeff <- numeric()
  nPixel <- numeric()
  for (z in 1:length(extractedI)){
    Icoeff[z] <- extractedI[[z]]
  }
  meanI <- mean(Icoeff,na.rm=T)
  medianI <- median(Icoeff,na.rm=T)
  maxI <- max(Icoeff,na.rm=T)
  minI <- min(Icoeff,na.rm=T)
  rangeI <- maxI - minI
  maxISq <- maxI^2
  minISq <- minI^2

  features <- list("Icoeff" = Icoeff, "meanI" = meanI, "medianI" = medianI, "maxI" = maxI, "minI" = minI, "rangeI" = rangeI,
                   "maxISq" = maxISq, "minISq" = minISq)

  return(features)
}


extract.Iy_xFixed <- function(Wlist,moranGray,testResults){
  # nCoords <- numeric()
  # extractedI <- list()
  # Npixel <- list()
  Islice <- numeric()
  for(z in 1:length(moranGray)){
    # for (y in 1:length(moranGray[[z]]$I_coeffList)){
    #   Npixel[[z]] <- list()
    #   extractedI[[z]] <- list()
    #   if (!is.null(dim(Wlist[[z]]$coords[[y]]))){
    #     Npixel[[z]][[y]] <- dim(Wlist[[z]]$coords[[y]])[1]
    #   }
    #   else {Npixel[[z]][[y]] <- NA}
    #   if(!is.na(Npixel[[z]][[y]])){
    #     if (!is.null(testResults$pValuesList[[z]][[y]]) && testResults$pValuesList[[z]][[y]] <= 0.05){
    #       extractedI[[z]][[y]] <- list("sliceId"=z, "IcoeffList"= moranGray[[z]]$I_coeffList[[y]], "pValueList"= testResults$pValuesList[[z]][[y]])
    #     }
    #     else {extractedI[[z]][[y]] <- NA}
    #   }
    #   else {extractedI[[z]][[y]] <- NA}
    # }
    Islice[z] <- mean(unlist(moranGray[[z]]$I_coeffList))

  }
  return(Islice)
}


calcola.statsTy_xFixed<-function(moranGray){

  Imean <- numeric()
  Icoeff <- list()
  Inull <- list()
  varNull <- list()
  zScoresList <- list()
  pValuesList <- list()

  for(z in seq(1:length(moranGray))){

    Inull[[z]]<-list()
    varNull[[z]]<-list()
    Icoeff[[z]]<-list()
    zScoresList[[z]] <- list()
    pValuesList[[z]] <- list()

    for(x in seq(1:length(moranGray[[z]]$I_coeffList))){
      if (length(moranGray[[z]]$I_coeffList) == 0 || is.null(moranGray[[z]]$I_coeffList[[x]])) {Icoeff[[z]][[x]] <- NA}
      else {Icoeff[[z]][[x]] <- moranGray[[z]]$I_coeffList[[x]]
      Inull[[z]][[x]] <- moranGray[[z]]$IList_nullHp[[x]]
      varNull[[z]][[x]] <- moranGray[[z]]$varIList_nullHp[[x]]
      zScoresList[[z]][[x]] <- (Icoeff[[z]][[x]] - Inull[[z]][[x]])/sqrt(varNull[[z]][[x]])
      pValuesList[[z]][[x]] <- 2*pnorm(-abs(zScoresList[[z]][[x]]))
      }
    }
  }

  testResults <- list("zScoresList"= zScoresList, "pValuesList"= pValuesList)
  #Imean <- mean(Icoeff)
  return(testResults)
}

####################

## A helper function that tests whether an object is either NULL _or_
## a list of NULLs
is.NullOb <- function(x) is.null(x) | all(sapply(x, is.null))

## Recursively step down into list, removing all such objects
rmNullObs <- function(x) {
  x <- Filter(Negate(is.NullOb), x)
  lapply(x, function(x) if (is.list(x)) rmNullObs(x) else x)
}

###################


calcola.moranGrayMeany_xFixed<-function(Wlist) {
  #ciclo sulle slice
  moranGray <-list()
  varIList <- list()
  for(z in seq(1:length(Wlist))){
    Wstar <- Wlist[[z]]$wStarList
    Wstar <- rmNullObs(Wstar)
    coords <- Wlist[[z]]$coordsList
    coords <- rmNullObs(coords)
    vcEroso <- Wlist[[z]]$vcErosoList
    vcEroso <- rmNullObs(vcEroso)
    InormList <- list()
    t_Wstar <-list()
    expectIList <- list()
    for(x in seq(1:length(coords))){
      if(length(coords)==0) next
      y <- coords[[x]]
      grayLevels <- numeric()
      #il minimo mi serve per portare a zero il valore minimo se negativo
      vcEroso[[x]] <- vcEroso[[x]][!is.na(vcEroso[[x]])]
      minumGray <- min(vcEroso[[x]])
      if(length(dim(y)[1])==0 || dim(y)[1] <= 5) next
      for(i in seq(1:dim(y)[1])){
        if (minumGray < 0 ) {grayLevels[i] <- vcEroso[[x]][i] - minumGray}
        else {grayLevels[i] <- vcEroso[[x]][i]}
      }
      grayMean <- mean(grayLevels,na.rm=T)
      # definisco la grandezza da "laggare" come differenza rispetto al valor medio
      grayVar <- grayLevels - grayMean
      #calcolo lo scalare di correlazione spaziale I
      laggedGray <- Wstar[[x]] %*% grayVar
      laggedGray <- laggedGray[which(!is.na(laggedGray))]
      grayVar <- grayVar[which(!is.na(laggedGray))]
      I <- grayVar %*% laggedGray
      norm <- grayVar %*% grayVar
      Inorm <- as.numeric(I/norm)
      InormList[[x]] <- Inorm
      #### CALCOLO ANCHE IL VALOR MEDIO E LA VARIANZA DI I SOTTO L'IPOTESI NULLA
      if (is.null(InormList[[x]]) == T) next
      t_Wstar[[x]] <- t(Wstar[[x]]) # 'a trassposta
      N <- nrow(Wstar[[x]])
      expectI = 1 / ( - N + 1 );
      S1 <- 0
      S2 <- 0
      S4 <- 0
      S5 <- 0
      sumSquaresW <- 0
      varI <- 0
      lista.varMoranNull<-.C("varMoranNull",
                             as.double(Wstar[[x]]),
                             as.double(t_Wstar[[x]]),
                             as.integer(N),
                             as.double(S1),
                             as.double(S2));

      S1 <- lista.varMoranNull[[4]]
      S2 <- lista.varMoranNull[[5]]
      S3 <-  nrow(Wstar[[x]]) * sum((grayVar^4)) * (1/sum(grayVar^2)^2)
      sumSquaresW <- sum(Wstar[[x]])^2
      S4 <- (N^2 - 3 * N + 3) * S1 - N * S2 + 3 * sumSquaresW
      S5 <- (N^2 - N) * S1 - 2 * N * S2 + 6 * sumSquaresW
      varI <- ((N * S4 - S3 * S5 ) / ((N - 1) * (N - 2) * (N - 3) * sumSquaresW)) - expectI^2
      expectIList[[x]] <- expectI
      varIList[[x]] <- varI
    }
    moranGray[[z]] <- list("I_coeffList"= InormList ,"IList_nullHp"= expectIList, "varIList_nullHp"= varIList)
  }
  return(moranGray)
}




calcola.matrice.Wy_xFixed<-function(voxelCube,dStar, erosionMarginX, erosionMarginY, vcNonEroso, vcEroso, listaVoxelErosi, pixelSpacingX, pixelSpacingY) {
  Wlist<-list()
  #ora seleziono le coord (x,y) per ogni slice z e calcolo la matrice delle distanze euclidee tra le coppie di punti (funzione: rdist) e la matrice W definita nel ciclo for(i)for(j)
  for(x in seq(1,max(listaVoxelErosi[,1]))){
    cat(c("."))
    slice <- listaVoxelErosi[which(listaVoxelErosi[,1]==x),]
    WstarList <- list()
    coordsList <- list()
    vcErosoList <- list()
    for(z in seq(1,max(listaVoxelErosi[,3]))){
      if(is.null(dim(slice))){slicey <- slice[which(slice[3]==z)]}
      else {slicey <- slice[which(slice[,3]==z),]}
      if (class(slicey) == "integer" || length(slicey) == 0) next
      slicey <- slicey*pixelSpacingY
      D <- matrix()
      D <- rdist(slicey)
      W <- matrix(0,nrow=nrow(D), ncol=ncol(D))
      Wstar <- matrix(0,nrow=nrow(W), ncol=ncol(W))
      # invoca la funzione in C forzando esplicitamente il tipo
      # dei parametri (importante)
      lista.risultato<-.C("internalLoop",
                          as.double(D),
                          as.integer(nrow(D)),
                          as.double(dStar),
                          as.matrix(W));
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
      #costruisco le k*z(k) (nPazienti*nslicePaziente) matrici W
      Wstar <- round(Wstar,2)
      WstarList[[z]] <- Wstar
      coordsList[[z]] <- slicey
      vcErosoList[[z]] <- vcEroso[x,,z]
    }

    Wlist[[x]] <- list("wStarList"= WstarList,"coordsList"= coordsList, "vcErosoList"= vcErosoList)

  }
  return(Wlist)
}
