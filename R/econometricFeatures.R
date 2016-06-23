#' cdscdscd
#'
#' @description  dfdsfdfdsfdsfds yuhyuhyu
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
#'
F01.pipeline.01<-function(inputGTV,dStar, erosionMarginX, erosionMarginY, whatDoYouWantBack=c("statsI")) {

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
      wMatrix<-calcola.matrice.W(voxelCube,dStar, erosionMarginX, erosionMarginY, vcNonEroso, vcEroso, listaVoxelErosi)
      if("wMatrix" %in% whatDoYouWantBack) final.result$wMatrix[[patID]]<-wMatrix

     # ora calcola il moranGrayMean
      mGM<-calcola.moranGrayMean(wMatrix)
      if("moranGrayMean" %in% whatDoYouWantBack) final.result$moranGrayMean[[patID]]<-mGM

     # ora calcola I
      stats.I<- calcola.statsT(mGM)
      if("statsI" %in% whatDoYouWantBack) final.result$statsI[[patID]]<-stats.I
    }
  }
  return(final.result)
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
  plot(slice,Icoeff)
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
    if (length(Wstar) == 0) next
    #calcolo lo scalare di correlazione spaziale I
    laggedGray <- Wstar %*% grayVar
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

calcola.matrice.W<-function(voxelCube,dStar, erosionMarginX, erosionMarginY, vcNonEroso, vcEroso, listaVoxelErosi) {
  Wlist<-list()
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

    Wlist[[z]] <- list("wStar"= Wstar,"coords"= slice, "vcEroso"=vcEroso[,,z])
  }
  return(Wlist)
}

