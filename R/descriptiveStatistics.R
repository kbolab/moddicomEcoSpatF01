#' Moran 2D spatial analysison DICOM images (extension of moddicomV2 package)
#' @useDynLib moddicomEcoSpatF01
#' @export slice.wise.moran.distribution
#' @export non.parametric.anisotrpy.test
#' @export single.directional.moran.scatterplot.xy
#' @export single.directional.moran.scatterplot.xz
#' @export single.directional.moran.scatterplot.yz
#' @export directional.moran.scatterplot
#' @export significant.moran.map



slice.wise.moran.distribution <- function(aaa_xy,aaa_xz,aaa_yx,boxplot=T){


  m_xy <- list()
  for (i in 1:length(aaa_xy$extractedI)){
    m_xy[[i]] <-  aaa_xy$extractedI[[i]]$Icoeff
  }

  allPat_xy <-  unlist(m_xy)
  allPat_xy_mean <- mean(allPat_xy,na.rm = T)
  allPat_xy_median <- median(allPat_xy,na.rm = T)
  allPat_xy_mode <- Mode(allPat_xy)

  m_xz <- list()
  for (i in 1:length(aaa_xz$extractedI)){
    m_xz[[i]] <-  aaa_xz$extractedI[[i]]
  }

  allPat_xz <-  unlist(m_xz)
  allPat_xz_mean <- mean(allPat_xz,na.rm = T)
  allPat_xz_median <- median(allPat_xz,na.rm = T)
  allPat_xz_mode <- Mode(allPat_xz)


  m_yz <- list()
  for (i in 1:length(aaa_yz$extractedI)){
    m_yz[[i]] <-  aaa_yz$extractedI[[i]]
  }

  allPat_yz <-  unlist(m_yz)
  allPat_yz_mean <- mean(allPat_yz,na.rm = T)
  allPat_yz_median <- median(allPat_yz,na.rm = T)
  allPat_yz_mode <- Mode(allPat_yz)


  res <- list("xy_mean"=allPat_xy_mean,"xy_median"=allPat_xy_median,"xy_mode"=allPat_xy_mode,
              "xz_mean"=allPat_xy_mean,"xz_median"=allPat_xy_median,"xz_mode"=allPat_xy_mode,
              "xz_mean"=allPat_xy_mean,"xz_median"=allPat_xy_median,"xz_mode"=allPat_xy_mode)

  if(boxplot==T){
    boxplot(allPat_xy,allPat_xz,allPat_yz,xaxt="n",xlab="",main="boxplot of Moran coefficients (all-slices)")
    axis(1, at=1:3, labels=c("xy","xz","yz"))
  }

  return(res)
}

non.parametric.anisotrpy.test <- function(aaa_xy,aaa_xz,aaa_yz){

  m_xy <- list()
  for (i in 1:length(aaa_xy$extractedI)){
    m_xy[[i]] <-  aaa_xy$extractedI[[i]]$Icoeff
  }
  allPat_xy <-  unlist(m_xy)

  m_xz <- list()
  for (i in 1:length(aaa_xz$extractedI)){
    m_xz[[i]] <-  aaa_xz$extractedI[[i]]
  }
  allPat_xz <-  unlist(m_xz)

  m_yz <- list()
  for (i in 1:length(aaa_yz$extractedI)){
    m_yz[[i]] <-  aaa_yz$extractedI[[i]]
  }
  allPat_yz <-  unlist(m_yz)

  d_xy <- density(allPat_xy,na.rm = T)
  d_xz <- density(allPat_xz,na.rm = T)
  d_yz <- density(allPat_yz,na.rm = T)

  plot(d_xz,col="red",xlim=c(-1,1),ylim=c(0,7),main="A simple non-parametric test of anisotropy")
  lines(d_yz,col="blue")
  lines(d_xy,col="green")
  legend("topleft", legend = c("xz plane","yz plane","xy plane"), col = c("red","blue","green"), lwd = 1,cex = 1)

  return()

}

single.directional.moran.scatterplot.xy <- function(aaa){

  rit_xy<-list()
  vox_xy<-list()
  for(pat in 1:length(aaa$wMatrix)){
    rit_xy[[pat]] <-list()
    vox_xy[[pat]] <-list()
    for(au in 1:length(aaa$wMatrix[[pat]])){
      if(length(aaa$wMatrix[[pat]][[au]]$wStar) != 0){
        rit_xy[[pat]][[au]] <- aaa$wMatrix[[pat]][[au]]$wStar %*% aaa$wMatrix[[pat]][[au]]$vcEroso[which(!is.na(aaa$wMatrix[[pat]][[au]]$vcEroso))]
        vox_xy[[pat]][[au]] <- aaa$wMatrix[[pat]][[au]]$vcEroso[which(!is.na(aaa$wMatrix[[pat]][[au]]$vcEroso))]
      }
    }
  }

  dati_xy <- as.data.frame(matrix(0,nrow = length(unlist(rit_xy)),ncol = 2))
  dati_xy[,1] <- unlist(vox_xy)
  dati_xy[,2] <- unlist(rit_xy)
  plot(dati_xy[,1],dati_xy[,2],xlab='original voxel intensity', ylab='spatially lagged voxel intensity',main="Directional Moran scatterplot xy")
  lm_xy <- lm(dati_xy[,2]~dati_xy[,1])
  abline(lm_xy$coefficients[1], lm_xy$coefficients[2],col="green")

  return(lm_xy)

}


single.directional.moran.scatterplot.xz <- function(aaa){

  rit_xz<-list()
  vox_xz<-list()
  for(pat in 1:length(aaa$wMatrix)){
    rit_xz[[pat]] <-list()
    vox_xz[[pat]] <-list()
    for(au in 1:length(aaa$wMatrix[[pat]])){
      rit_xz[[pat]][[au]] <- list()
      vox_xz[[pat]][[au]] <- list()
      for(cazz in 1:length(aaa$wMatrix[[pat]][[au]]$wStarList))
        if(length(aaa$wMatrix[[pat]][[au]]$wStarList[[cazz]]) != 0){
          rit_xz[[pat]][[au]][[cazz]] <- aaa$wMatrix[[pat]][[au]]$wStarList[[cazz]] %*% aaa$wMatrix[[pat]][[au]]$vcErosoList[[cazz]][which(!is.na(aaa$wMatrix[[pat]][[au]]$vcErosoList[[cazz]]))]
          vox_xz[[pat]][[au]][[cazz]] <- aaa$wMatrix[[pat]][[au]]$vcErosoList[[cazz]][which(!is.na(aaa$wMatrix[[pat]][[au]]$vcErosoList[[cazz]]))]
        }
    }
  }

  dati_xz <- as.data.frame(matrix(0,nrow = length(unlist(rit_xz)),ncol = 2))
  dati_xz[,1] <- unlist(vox_xz)
  dati_xz[,2] <- unlist(rit_xz)
  plot(dati_xz[,2],dati_xz[,1],xlab='original voxel intensity', ylab='spatially lagged voxel intensity',main="Directional Moran scatterplot xz")
  lm_xz <- lm(dati_xz[,1]~dati_xz[,2])
  abline(lm_xz$coefficients[1], lm_xz$coefficients[2],col="red")

  return(lm_xz)
  }


single.directional.moran.scatterplot.yz <- function(aaa){

  rit_yz<-list()
  vox_yz<-list()
  for(pat in 1:length(aaa$wMatrix)){
    rit_yz[[pat]] <-list()
    vox_yz[[pat]] <-list()
    for(au in 1:length(aaa$wMatrix[[pat]])){
      rit_yz[[pat]][[au]] <- list()
      vox_yz[[pat]][[au]] <- list()
      for(cazz in 1:length(aaa$wMatrix[[pat]][[au]]$wStarList))
        if(length(aaa$wMatrix[[pat]][[au]]$wStarList[[cazz]]) != 0){
          rit_yz[[pat]][[au]][[cazz]] <- aaa$wMatrix[[pat]][[au]]$wStarList[[cazz]] %*% aaa$wMatrix[[pat]][[au]]$vcErosoList[[cazz]][which(!is.na(aaa$wMatrix[[pat]][[au]]$vcErosoList[[cazz]]))]
          vox_yz[[pat]][[au]][[cazz]] <- aaa$wMatrix[[pat]][[au]]$vcErosoList[[cazz]][which(!is.na(aaa$wMatrix[[pat]][[au]]$vcErosoList[[cazz]]))]
        }
    }
  }

  dati_yz <- as.data.frame(matrix(0,nrow = length(unlist(rit_yz)),ncol = 2))
  dati_yz[,1] <- unlist(vox_yz)
  dati_yz[,2] <- unlist(rit_yz)
  plot(dati_yz[,2],dati_yz[,1],xlab='original voxel intensity', ylab='spatially lagged voxel intensity',main="Directional Moran scatterplot yz")
  lm_yz <- lm(dati_yz[,1]~dati_yz[,2])
  abline(lm_yz$coefficients[1], lm_yz$coefficients[2],col="blue")

  return(lm_yz)

}

directional.moran.scatterplot <- function(lm_xy,lm_xz,lm_yz){

  plot(1,type='n',xlim=c(1,600),ylim=c(0,600),xlab='original voxel intensity', ylab='spatially lagged voxel intensity',main="Directional Moran scatteplot")
  abline(lm_xz$coefficients[1], lm_xz$coefficients[2],col="red")
  abline(lm_yz$coefficients[1], lm_yz$coefficients[2],col="blue")
  abline(lm_xy$coefficients[1], lm_xy$coefficients[2],col="green")
  legend("bottomright", legend = c("xz plane","yz plane","xy plane"), col = c("red","blue","green"), lwd = 1,cex = 1)

}

significant.moran.map <- function(aaa_xy,aaa_xz,aaa_yz,image_number=1,threshold=0.05){

  p_xz <- list()
  pvalue_xz <- list()
  for(k in 1:length(aaa_xz$statsI)){
    p_xz[[k]] <- list()
    for(j in 1:length(aaa_xz$statsI[[k]]$pValuesList)){
      p_xz[[k]][[j]] <- median(unlist(aaa_xz$statsI[[k]]$pValuesList[[j]]),na.rm = T)
    }
    pvalue_xz[[k]] <- unlist(p_xz[[k]])
  }

  p_yz <- list()
  pvalue_yz <- list()
  for(k in 1:length(aaa_yz$statsI)){
    p_yz[[k]] <- list()
    for(j in 1:length(aaa_xz$statsI[[k]]$pValuesList)){
      p_yz[[k]][[j]] <- median(unlist(aaa_yz$statsI[[k]]$pValuesList[[j]]),na.rm = T)
    }
    pvalue_yz[[k]] <- unlist(p_yz[[k]])
  }

  pvalue_xy <- list()
  for(k in 1:length(aaa_xy$statsI)){
    pvalue_xy[[k]] <- aaa_xy$statsI[[k]]$pValues
  }

  #DRAW GRAPH (for a single image only)
  pValuesM <- cbind(pvalue_xy[[image_number]],pvalue_yz[[image_number]],pvalue_xz[[image_number]])
  pValuesM.copy <- matrix(FALSE,ncol = ncol(pValuesM),nrow = nrow(pValuesM))
  pValuesM.copy[which(pValuesM <= threshold)]<- TRUE
  z <- as.logical(levels(as.factor(pValuesM.copy)))
  image(x=c(0:2),y=c(1:dim(pValuesM.copy)[1]),z=t(pValuesM.copy),col=z,main="Significance of slice-by-slice moran (black if p<=0.05)",ylab="slice number", xaxt="n",xlab="")
  axis(1, at=0:2, labels=c("xy","yz","xz"))
  grid(nx = dim(pValuesM.copy)[2], ny = dim(pValuesM.copy)[1], col = "lightgray", lty = "solid",lwd = par("lwd"), equilogs = FALSE)

}

#########################################################
Mode <- function(x) {
  ux <- unique(x)
  ux <- ux[which(!is.na(ux))]
  ux[which.max(tabulate(match(x, ux)))]
}
#########################################################

