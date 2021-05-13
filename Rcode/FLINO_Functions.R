
library("plyr")

## https://www.rdocumentation.org/packages/NOISeq/versions/2.16.0/topics/Normalization
# BiocManager::install("NOISeq")
# BiocManager::install("fCI")
library(fCI)  
library(NOISeq)
library(qsmooth)




medianRemoveNAs = function(x) { return (median(x,na.rm = TRUE))}

meanRemoveNAs = function(x) { return (mean(x,na.rm = TRUE))}

sdRemoveNAs = function(x) { return (sd(x,na.rm = TRUE))}

countNonNAs = function(x) { return (length(x[!is.na(x)]))}

rmseRemoveNAs = function(x) { 
  mx = mean(x,na.rm = TRUE)
  return (sqrt(mean((x - mx)^2,na.rm = TRUE)))
}

cvRemoveNAs = function(x) { 
  x_sd = sd(x,na.rm = TRUE)
  x_mean = mean(x,na.rm = TRUE)
  return (x_sd / x_mean)
}


log10Tck <- function(side, type){
  lim <- switch(side, 
                x = par('usr')[1:2],
                y = par('usr')[3:4],
                stop("side argument must be 'x' or 'y'"))
  at <- floor(lim[1]) : ceiling (lim[2])
  return(switch(type, 
                minor = outer(1:9, 10^(min(at):max(at))),
                major = 10^at,
                stop("type argument must be 'major' or 'minor'")
  ))
}


plotslideView = function(slideDat, channel) {
  
  par(mfrow=c(3,2))
  
  #slideDat = merge(slideDat,slide_position_patient_id,by=c("slide","position"),all.x=TRUE)
  slideDat = join(slideDat,slide_position_patient_id,by=c("slide","position"))
  
  for(slide in c("S18030274","S18030275","S18030276")) {
    
    
    db = slideDat[slideDat[,"slide"] == slide, ]
    
    gw = 0.6
    
    db[,"scaled"] =  0.9 * (db[,"marker"] -min(slideDat[,"marker"], na.rm = TRUE)) / (max(slideDat[,"marker"], na.rm = TRUE) -min(slideDat[,"marker"], na.rm = TRUE))
    db[ ,"Missing"] = FALSE
    db[is.na(db[,"scaled"]),"Missing"] = TRUE
    ii = which(!is.na(db[,"scaled"]))
    db[ii,"color"] =  rgb(db[ii ,"scaled"],db[ii ,"scaled"],db[ii ,"scaled"],1)
    ii = which(is.na(db[,"scaled"]))
    db[ii,"color"] =  "yellow"
    if (length(ii) > 0) {
      print(paste0(marker, " - " , slide))
      print(ii)
    }
    
    slideScaled =  0.9 * (median(db[,"marker"], na.rm = TRUE) -min(slideDat[,"marker"], na.rm = TRUE)) / (max(slideDat[,"marker"], na.rm = TRUE) -min(slideDat[,"marker"], na.rm = TRUE))
    slideColor =  rgb(slideScaled,slideScaled,slideScaled,1)
    
    
    
    plot(db[,"x"],-db[,"y"],pch=NA,main=paste0(marker, " - " , slide), xlim=c(min(slideDat[,"x"], na.rm = TRUE)-2*gw,max(slideDat[,"x"], na.rm = TRUE)+2*gw),ylim=c(-max(slideDat[,"y"], na.rm = TRUE)-2*gw,-min(slideDat[,"y"], na.rm = TRUE)+2*gw), xlab="",ylab="")
    rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = slideColor)
    for(i in 1:nrow(db)) {
      positionCol = db[,"color"]
      rect(db[,"x"]-gw, -db[,"y"]-gw, db[,"x"]+gw, -db[,"y"]+gw, col = positionCol,border = positionCol)
    }
    text(db[,"x"],-db[,"y"],db[,"position"],col='white')
    
  }
  
  par(mfrow=c(1,1))
  
}


setDefaultPlotParm = function() {
  acex <<- 0.9
  tcex <<-  1
  lcex <<-  0.75
  cexLegend <<- 1
  SHOW_LEGEND <<-  FALSE
  t1cex <<- 1
  t1bcex <<- 0.9
  t2cex <<- 0.9
  t2bcex <<- 0.8
  t3cex <<- 0.8
  symSize <<- 0.7
  sympch <<- 20
  symcol <<- "lightsteelblue4"
  symJitter <<- 0.15
  ADD_HORIZONTAL_MIN_BAR_VALUE <<- TRUE
  DISPLAY_DATA_POINTS <<- TRUE
  DISPLAY_DATA_LINE_SEGMENTS <<- FALSE
  DISPLAY_ERROR_BARS <<- FALSE
  DISPLAY_QUANTILE_BARS <<- TRUE
  DISPLY_GRID <<- TRUE
  DRAW_BOX_FRAME <<- TRUE
  USE_STANDARD_DEVIATION_FOR_ERROR_BARS <<- TRUE
  
  ADD_HORIZONTAL_BAR_AT_GIVEN_VALUE  <<- FALSE
  HORIZONTAL_BAR_VALUE <<- 0.1
  MEDIAN_BAR_VERSION <<- TRUE
}



plotBMdata = function(dat, marker, pd_ylim = NA, cellPositions = NA) {
  
  
  otherPositions = unique(dat[,"position"])
  otherPositions = otherPositions[!otherPositions %in% cellPositions]
  
  # myXlab= paste0(marker, " Round ", round , " ", channel ," Channel")
  # myYlab=paste0(marker, " (FOV Median (Grid Mean pixel intensity))")
  if (length(pd_ylim) < 2) {
    boxplot(marker~slide, data=dat, frame=FALSE, names=getSlideNames(slides),yaxt='n', outpch = NA,xlab="",ylab="")
    axis(2,las=2)
  } else {
    boxplot(marker~slide, data=dat, frame=FALSE, names=getSlideNames(slides),yaxt='n', outpch = NA, ylim=pd_ylim,xlab="",ylab="")
    axis(2,las=2)
  }
  
  
  jitterMethod = "jitter"
  for(si in 1:length(slides)) {
    stripchart(dat[dat[,"slide"] == slides[si] & dat[,"position"] %in% otherPositions,"marker"],at =si, vertical = TRUE, method = jitterMethod,
               pch = 1, col = "grey20",
               add = TRUE,cex=1)
  }
  
  for(si in 2:length(slides)) {
    for (position in cellPositions) {
      x0 = si-1
      x1 = si
      y0 = dat[dat[,"slide"] == slides[si-1] & dat[,"position"] == position,"marker"]
      y1 = dat[dat[,"slide"] == slides[si] & dat[,"position"] == position,"marker"]
      if (length(y0) > 0 & length(y1) > 0) {
        segments(x0, y0, x1, y1, col = "blue", lty = par("lty"), xpd = FALSE)
      }
      
      
    }
  }
  
  
  jitterMethod = "overplot" # "stack" 
  for(si in 1:length(slides)) {
    stripchart(dat[dat[,"slide"] == slides[si] & dat[,"position"] == cellPositions[1],"marker"],at =si, vertical = TRUE, method = jitterMethod,
               pch = 21, bg='white', col = "blue",
               add = TRUE,cex=1.8)
    
    stripchart(dat[dat[,"slide"] == slides[si] & dat[,"position"] == cellPositions[2],"marker"],at =si, vertical = TRUE, method = jitterMethod,
               pch = 23, bg='white', col = "blue",
               add = TRUE,cex=1.8)
    
    stripchart(dat[dat[,"slide"] == slides[si] & dat[,"position"] == cellPositions[3],"marker"],at =si, vertical = TRUE, method = jitterMethod,
               pch = 22, bg='white', col = "blue",
               add = TRUE,cex=1.8)
    
    stripchart(dat[dat[,"slide"] == slides[si] & dat[,"position"] == cellPositions[4],"marker"],at =si, vertical = TRUE, method = jitterMethod,
               pch = 24, bg='white', col = "blue",
               add = TRUE,cex=1.8)
  }
  
  
  
}


normalizeNonzeroByStrata<-function(markerdata,INDEX,FUN){
  ## Coerce markerdata to data frame if not a data frame or matrix
  if(!is.data.frame(markerdata) & !is.matrix(markerdata)) markerdata<-as.data.frame(markerdata)
  nr<-nrow(markerdata)
  nc<-ncol(markerdata)
  ## Define INDEX if missing
  if(is.null(INDEX)) INDEX=rep(1,nr)
  FUN<-match.fun(FUN)
  
  ## Split row indices by INDEX
  id<-1:nr
  idl<-split(id,INDEX)
  ## Normalize by strata on nonzero values by columns
  res<-do.call("cbind",lapply(1:nc,function(i){
    m<-markerdata[,i]
    ## Compute overall statistic
    f0<-FUN(m[which(m!=0)])
    ## Normalize by strata
    ret<-do.call("rbind",lapply(idl,function(j){
      m1<-m[j]
      id0<-which(m1!=0)
      m1[id0]<-m1[id0]-FUN(m1[id0])+f0
      m2<-cbind(j,m1)
      m2
    }))
    rownames(ret)<-NULL
    ret<-ret[order(ret[,1]),2]
    ret
  }))
  rownames(res)<-rownames(markerdata)
  res<-data.frame(res)
  names(res)<-colnames(markerdata)
  return(res)
}


get_markers = function(dat.marker_SEG) {
  cns = colnames(dat.marker_SEG)
  cns = cns[startsWith(cns,"Mean.AF000.")]
  markers = str_split_fixed(cns,"\\.",4)[,4]
  markers = markers[order(markers)]
  return(markers)
}

get_bm = function(prefix = "Mean.AF000.", markers) {
  bm = paste0(prefix ,markers,".",markers)
  return (bm)
}


get_dat0 = function(dat.marker_SEG, dat.Corr_SEG, cy3rounds, cy5rounds, bm) {
  
  dat0 = dat.marker_SEG
  
  cm = colnames(dat.Corr_SEG)
  cm = cm[endsWith(cm,".Corr.Corr")]
  
  if (length(cm) > 0) {
    #dat0 = merge(dat0,dat.Corr_SEG[,c("slide","position","Cell.ID",cm)],by=c("slide","position","Cell.ID" ),all.x=TRUE)
    dat0 = join(dat0,dat.Corr_SEG[,c("slide","position","Cell.ID",cm)],by=c("slide","position","Cell.ID" ))
    
    cm2 = cm
    cm2 = gsub("Mean.RG","",cm2)
    cm2 = gsub(".Corr.Corr","",cm2)
    cm2 = as.numeric(cm2)
    markerRoundDapiCorrLookup = rbind(cy3rounds,cy5rounds)
    markerRoundDapiCorrLookup = unique(markerRoundDapiCorrLookup)
    #markerRoundDapiCorrLookup = merge(markerRoundDapiCorrLookup,data.frame(CorrColumnName=cm,Round=cm2),by="Round",all.x=TRUE)
    markerRoundDapiCorrLookup = join(markerRoundDapiCorrLookup,data.frame(CorrColumnName=cm,Round=cm2),by="Round")
    
    rm(cm2)
  } else {
    markerRoundDapiCorrLookup = data.frame()
  }
  
  return(dat0)
}



quantileNormalization = function(dat, quantileProb = 0.75, FILTER_ZEROS = TRUE) {
  if (FILTER_ZEROS) {
    q3 <- apply(dat, 2, function(x){quantile(x[x>0], quantileProb)})
  } else {
    q3 <- apply(dat, 2, function(x){quantile(x, quantileProb)})
  }
  dat.norm <- t(t(dat) / q3)

  #uqua method
  # d <- q3 * sum(rowSums(dat))/sum(q3)
  # dat.norm <- t(t(dat)/d) * 10^6

  return (dat.norm)
}

sincerosModified = function (datos, k) 
{
  datos = as.matrix(datos)
  datos0 <- as.matrix(datos)
  if (is.null(k)) {
    mini0 <- min(datos[noceros(datos, num = FALSE, k = 0)])
    kc <- mini0/2
    datos0[datos0 == 0] <- kc
  }
  else {
    datos0[datos0 == 0] <- k
  }
  datos0
}

uquaModified = function (datos, quantileProb = 0.75, long = 1000, lc = 0, k = 0) 
{
  if (!is.null(ncol(long))) {
    mynames = long[, 1]
    long = long[, 2]
    names(long) = mynames
  }
  L <- (long/1000)^lc
  datos = datos/L
  datos0 <- sincerosModified(datos, k)
  if (ncol(as.matrix(datos)) > 1) {
    sumatot <- rowSums(datos)
    supertot <- sum(sumatot)
    counts0 <- which(sumatot == 0)
    if (length(counts0) > 0) {
      datitos <- datos[-counts0, ]
    }
    else {
      datitos <- datos
    }
    q3 <- apply(datitos, 2, quantile, probs = quantileProb)
    d <- q3 * supertot/sum(q3)
    datos.norm <- t(t(datos0)/d) * 10^6
  }
  else {
    datos.norm <- datos0/L
  }
  na.omit(datos.norm)
}

               




quantile75 = function(x) { return (quantile(x)[4])}

get_datNorm = function(dat0, BMdat, method = "TMM", SCALE_NORM_TO_RAW_MEDIAN=TRUE, UNEVEN_MATRIX_LENGTH_METHOD = "SAMPLING_TO_MIN_LENGTH") {
  
  if (method == "MEDIAN" | 
      method == "Q2NORM" | 
      method == "Q3NORM" | 
      method == "Q375NORM" | 
      method == "Q625NORM" | 
      method == "Q875NORM" |
      method == "Q90NORM"  |
      method == "Q95NORM"  |
      method == "Q99NORM"  |
      method == "MAX") {
    # These methods can work on data in which the number
    # of normalization objects is not balanced across slides
    if (method == "Q3NORM") {
      BMdat.norm=normalizeNonzeroByStrata(markerdata=BMdat,INDEX=dat0$slide,FUN=quantile75)
    } else if (method == "Q2NORM" | method == "MEDIAN") {
      BMdat.norm=normalizeNonzeroByStrata(markerdata=BMdat,INDEX=dat0$slide,FUN=median)
    } else if (method == "Q375NORM") {
      myFUNC = function(x) { return (quantile(x,probs = c(0,0.375,1))[2])}
      BMdat.norm=normalizeNonzeroByStrata(markerdata=BMdat,INDEX=dat0$slide,FUN=myFUNC)
    } else if (method == "Q625NORM") {
      myFUNC = function(x) { return (quantile(x,probs = c(0,0.625,1))[2])}
      BMdat.norm=normalizeNonzeroByStrata(markerdata=BMdat,INDEX=dat0$slide,FUN=myFUNC)
    } else if (method == "Q875NORM") {
      myFUNC = function(x) { return (quantile(x,probs = c(0,0.875,1))[2])}
      BMdat.norm=normalizeNonzeroByStrata(markerdata=BMdat,INDEX=dat0$slide,FUN=myFUNC)
    } else if (method == "Q90NORM") {
      myFUNC = function(x) { return (quantile(x,probs = c(0,0.90,1))[2])}
      BMdat.norm=normalizeNonzeroByStrata(markerdata=BMdat,INDEX=dat0$slide,FUN=myFUNC)
    } else if (method == "Q95NORM") {
      myFUNC = function(x) { return (quantile(x,probs = c(0,0.95,1))[2])}
      BMdat.norm=normalizeNonzeroByStrata(markerdata=BMdat,INDEX=dat0$slide,FUN=myFUNC)
    } else if (method == "Q99NORM") {
      myFUNC = function(x) { return (quantile(x,probs = c(0,0.99,1))[2])}
      BMdat.norm=normalizeNonzeroByStrata(markerdata=BMdat,INDEX=dat0$slide,FUN=myFUNC)
    } else if (method == "MAX") {
      BMdat.norm=normalizeNonzeroByStrata(markerdata=BMdat,INDEX=dat0$slide,FUN=max)
    } else {
      BMdat.norm = BMdat
      cat("-------------------------------", "\n")
      cat("UNKNOWN NORM METHOD ", elist[["PARM_NORM_METHOD"]], "\n")
    }

  } else {

    bm = colnames(BMdat)
    mBMdat = cbind(dat0[,c("slide","position","Cell.ID")],BMdat)
    mBMdat[,"order"] = c(1:nrow(mBMdat))
    oBMdat = mBMdat[,c("slide","position","Cell.ID","order")]
    for(nm in 1:ncol(BMdat)) {
      
      # determine the number of object IDs to use for each slide
      slideSegIDs = aggregate(mBMdat[,"Cell.ID"] ,by=list(mBMdat[,"slide"]),FUN=length)
      minSegIDs = min(slideSegIDs[,"x"])
      maxSegIDs = max(slideSegIDs[,"x"])
      vibase = c()
      
      slides = unique(mBMdat[,"slide"])
      d1 = mBMdat[mBMdat[,"slide"] == slides[1], c("position", "Cell.ID", bm[nm]) ]
      colnames(d1)[1] = paste0(slides[1],"_",colnames(d1)[1])
      colnames(d1)[2] = paste0(slides[1],"_",colnames(d1)[2])
      colnames(d1)[3] = slides[1]
      d1[,3] = as.numeric(d1[,3])
      if (UNEVEN_MATRIX_LENGTH_METHOD == "SAMPLING_TO_MIN_LENGTH") {
        d2 = d1[sample(c(1:nrow(d1)), size=minSegIDs, replace =FALSE), ]
      } else if (UNEVEN_MATRIX_LENGTH_METHOD == "SAMPLING_TO_MAX_LENGTH") {
        addNum = maxSegIDs - nrow(d1)
        vibase = c(vibase,c(1:nrow(d1)))
        d2 = rbind(d1,d1[sample(c(1:nrow(d1)), size=addNum, replace =TRUE), ])
      } else if (UNEVEN_MATRIX_LENGTH_METHOD == "FIRST_ORDERING_MIN_LENGTH") {
        d2 = d1[c(1:minSegIDs), ]
      }
      for (slide in slides[-1]) {
        
        d1 = mBMdat[mBMdat[,"slide"] == slide, c("position", "Cell.ID",bm[nm]) ]
        d1[,3] = as.numeric(d1[,3])
        colnames(d1)[1] = paste0(slide,"_",colnames(d1)[1])
        colnames(d1)[2] = paste0(slide,"_",colnames(d1)[2])
        colnames(d1)[3] = slide

        if (UNEVEN_MATRIX_LENGTH_METHOD == "SAMPLING_TO_MIN_LENGTH") {
          d1 = d1[sample(c(1:nrow(d1)), size=minSegIDs, replace =FALSE), ]
        } else if (UNEVEN_MATRIX_LENGTH_METHOD == "SAMPLING_TO_MAX_LENGTH") {
          addNum = maxSegIDs - nrow(d1)
          vibase = c(vibase,c((nrow(d2)+1):(nrow(d2) + nrow(d1))))
          d1 = rbind(d1,d1[sample(c(1:nrow(d1)), size=addNum, replace =TRUE), ])
        } else if (UNEVEN_MATRIX_LENGTH_METHOD == "FIRST_ORDERING_MIN_LENGTH") {
          d1 = d1[c(1:minSegIDs), ]
        }
        
        d2 = cbind(d2,d1)
      }
      
      
      if (method == "MRN") {
        # Median Ratio Normalization (MRN)
        # library(fCI)  
        # deseq.median.ratio.normalization
        d1 <- deseq.median.ratio.normalization(d2[,slides])
      } else if (method == "TMM") {
        # Trimmed Mean of M (TMM) normalization (Robinson and Oshlack, 2010)
        # library(NOISeq)
        d1 <- tmm(d2[,slides])
      } else if (method == "UPPER_QUARTILE") {
        # Upper Quartile  normalization(Bullard 2010)
        # library(NOISeq)
        d1 <- uqua(d2[,slides])
      } else if (method == "UQUA_MODIFIED") {
        # Upper Quartile  normalization(Bullard 2010)
        # library(NOISeq)
        d1 <- uquaModified(d2[,slides])
      } else if (method == "Q75N") {
        # quantile I-norm = alpha I-raw  normalization
        d1 <- quantileNormalization(d2[,slides], quantileProb = 0.75, FILTER_ZEROS = FALSE)
      } else if (method == "Q25NZ") {
        # quantile I-norm = alpha I-raw  normalization
        d1 <- quantileNormalization(d2[,slides], quantileProb = 0.25, FILTER_ZEROS = TRUE)
      } else if (method == "Q375NZ") {
        # quantile I-norm = alpha I-raw  normalization
        d1 <- quantileNormalization(d2[,slides], quantileProb = 0.375, FILTER_ZEROS = TRUE)
      } else if (method == "Q50NZ") {
        # quantile I-norm = alpha I-raw  normalization
        d1 <- quantileNormalization(d2[,slides], quantileProb = 0.5, FILTER_ZEROS = TRUE)
      } else if (method == "Q625NZ") {
        # quantile I-norm = alpha I-raw  normalization
        d1 <- quantileNormalization(d2[,slides], quantileProb = 0.625, FILTER_ZEROS = TRUE)
      } else if (method == "Q75NZ") {
        # quantile I-norm = alpha I-raw  normalization
        d1 <- quantileNormalization(d2[,slides], quantileProb = 0.75, FILTER_ZEROS = TRUE)
      } else if (method == "Q875NZ") {
        # quantile I-norm = alpha I-raw  normalization
        d1 <- quantileNormalization(d2[,slides], quantileProb = 0.875, FILTER_ZEROS = TRUE)
      } else if (method == "Q90NZ") {
        # quantile I-norm = alpha I-raw  normalization
        d1 <- quantileNormalization(d2[,slides], quantileProb = 0.90, FILTER_ZEROS = TRUE)
      } else if (method == "Q95NZ") {
        # quantile I-norm = alpha I-raw  normalization
        d1 <- quantileNormalization(d2[,slides], quantileProb = 0.95, FILTER_ZEROS = TRUE)
      } else if (method == "Q99NZ") {
        # quantile I-norm = alpha I-raw  normalization
        d1 <- quantileNormalization(d2[,slides], quantileProb = 0.99, FILTER_ZEROS = TRUE)
      } else if (method == "QMAXNZ") {
        # quantile I-norm = alpha I-raw  normalization
        d1 <- quantileNormalization(d2[,slides], quantileProb = 1, FILTER_ZEROS = TRUE)
      } else if (method == "RPKM") {
        # RPKM normalization(Mortazavi, 2008)
        # library(NOISeq)
        d1 <- rpkm(d2[,slides])
      } else if (method == "QUANTILE") {
        # Smooth quantile normalization (Hicks 2018)
        # Hicks 2018 
        if (length(slides) == 2) {
          dd = d2[,slides]
          dd[,"dupSlide1"] = dd[,1]
          dd[,"dupSlide2"] = dd[,2]
          qs <- qsmooth(object = dd, group_factor = c(1:ncol(dd)))
          d1 = data.frame(qsmoothData(qs)) 
          d1 = d1[,slides]
        } else {
          qs <- qsmooth(object = d2[,slides], group_factor = c(1:ncol(d2[,slides])))
          d1 = data.frame(qsmoothData(qs)) 
        }
        
      } else {
        d1 = d2[,slides]
        cat("-------------------------------", "\n")
        cat("UNKNOWN NORM METHOD ", elist[["PARM_NORM_METHOD"]], "\n")
      }

      colnames(d1) = paste0(colnames(d1),"_","norm")
      d2 = cbind(d2,d1)

      d1 = d2[,c(paste0(slides[1],"_","position"), paste0(slides[1],"_","Cell.ID"),paste0(slides[1],"_","norm"))]
      colnames(d1)[1] = "position"
      colnames(d1)[2] = "Cell.ID"
      colnames(d1)[3] = bm[nm]
      d1[,"slide"] = slides[1]
      d3 = d1
      for (slide in slides[-1]) {
        d1 = d2[,c(paste0(slide,"_","position"), paste0(slide,"_","Cell.ID"),paste0(slide,"_","norm"))]
        colnames(d1)[1] = "position"
        colnames(d1)[2] = "Cell.ID"
        colnames(d1)[3] = bm[nm]
        d1[,"slide"] = slide
        d3 = rbind(d3,d1)
      }
      
      if (UNEVEN_MATRIX_LENGTH_METHOD == "SAMPLING_TO_MAX_LENGTH") {
        d3 = d3[vibase,]
      }
      
      #oBMdat = merge(oBMdat,d3,by=c("slide","position","Cell.ID"),all.x=TRUE)
      oBMdat = join(oBMdat,d3,by=c("slide","position","Cell.ID"))
      
    }
    
    oBMdat = oBMdat[order(oBMdat[,"order"]), ]
    row.names(oBMdat) = NULL
    BMdat.norm = oBMdat[,bm]
  }

  if (SCALE_NORM_TO_RAW_MEDIAN) {
    for(nm in 1:ncol(BMdat)) {
      nf = median(BMdat[,nm],na.rm=TRUE) / median(BMdat.norm[,nm],na.rm=TRUE)
      BMdat.norm[,nm] = BMdat.norm[,nm] * nf
    }
  }
  
  BMdat.norm = extrapolateNormalization(dat0, BMdat, BMdat.norm, DISPLAY_FITS = FALSE)
  
  return (BMdat.norm)
}


calcIntensityCorrectionFit = function(pd1, rawCi, normCi, fxlim, fylim, slide, DISPLAY_FITS, polynomial_degree = 2) {
  
  #Note correction fit is likely to occur in log space, but could be absolute space
  # Min pixel intensity is 1 so log of 1 is zero and for absolute scale its pretty close to zero
  
  # Before adding the origin point compute the range coverage factor
  # rangeCoverageFactor = (maxValue - minValue) /(maxValue - zero)
  rangeCoverageFactor = 1 - min(pd1[,"x"]) / max(pd1[,"x"])

  # add point at origin
  i = nrow(pd1) + 1
  pd1[i,"x"] = 0
  pd1[i,"y"] = 0
  
  if (polynomial_degree == 1 || rangeCoverageFactor < 0.1 || nrow(unique(pd1)) < 3) {
    fit <- lm(y ~ x, data = pd1)
  } else if (polynomial_degree == 2) {
    fit <- lm(y ~ x + I(x^2), data = pd1)
  } else if (polynomial_degree == 3) {
    fit <- lm(y ~ x + I(x^2) + I(x^3), data = pd1)
  }
 
  
  #print(summary(fit))

  if (DISPLAY_FITS) {
    cat("--------------------------------","\n")
    cat(slide,"\n")
    cat("--------------------------------","\n")
    print(summary(fit))
    newdat = data.frame(x = seq(min(pd1$x), max(pd1$x), length.out = 100))
    newdat$pred = predict(fit, newdata = newdat)
    plot(y ~ x,data=pd1,xlab=rawCi,ylab=normCi, main=slide, xlim=fxlim,ylim=fylim)
    grid()
    with(newdat, lines(x = x, y = pred, col='red',lwd=2))

  }
  
  return(fit)
  
}


# This function looks will use pair values found in both BMdat and BMdat.norm to build empirical extrapolation functions for each slide or slide-position
# These functions will then be used to infer values for missing pairs (e.g. BMdat has a value but it is missing in BMdat.norm)
getNormalizationFunctions = function(dat0, BMdat, BMdat.norm, cbm, DISPLAY_FITS = FALSE, polynomial_degree = 2) {
  

  dbr = cbind(dat0[,c("slide","position", "Cell.ID")],BMdat)
  dbn = cbind(dat0[,c("slide","position", "Cell.ID")],BMdat.norm)

  dbdif = dbn[,c("slide","position", "Cell.ID")]


  rawCi = paste0("Raw.",gsub("Mean.AF000.","",cbm))
  normCi = paste0("Norm.",gsub("Mean.AF000.","",cbm))
  
  dbdif[,paste0("NormMinusRaw.",gsub("Mean.AF000.","",cbm))] = dbn[,cbm] - dbr[,cbm]
  dbdif[,rawCi] = dbr[,cbm]
  dbdif[,normCi] = dbn[,cbm]
  
  SlideNormFuncList = list()
  pd1 = na.omit(dbdif[,c(rawCi,normCi)])
  fxlim=c(0,max(pd1[,1]))
  fylim=c(0,max(pd1[,2]))

  dslides = unique(dbdif[,"slide"])

  for (slide in dslides) {
    d = dbdif
    d = d[d[,"slide"] == slide, ]
    pd1 = na.omit(d[,c(rawCi,normCi)])
    # plot(pd1, main=slide, xlim=fxlim,ylim=fylim)
    colnames(pd1) = c("x","y")
    
    fit <- calcIntensityCorrectionFit(pd1, rawCi, normCi, fxlim, fylim, slide, DISPLAY_FITS, polynomial_degree)

    SlideNormFuncList[[slide]] = fit
  }

  
  return (SlideNormFuncList)
}


# This function looks will use pair values found in both BMdat and BMdat.norm to build empirical extrapolation functions for each slide or slide-position
# These functions will then be used to infer values for missing pairs (e.g. BMdat has a value but it is missing in BMdat.norm)
applyNormalizationFunctions = function(dat0, BMdat, BMdat.norm, cbm, SlideNormFuncList) {
  
  
  dslides = names(SlideNormFuncList)
  
  for (slide in dslides) {
    fit =SlideNormFuncList[[slide]]
    
    # for each entry that is !NA in BMdat use prediction
    mi = which(!is.na(BMdat[,cbm]) & dat0[,"slide"] == slide) 
    if (length(mi) > 0) {
      
      pvs = predict(fit, newdata = data.frame(x = BMdat[mi,cbm]))
      pvs[!is.na(pvs) & pvs < 0] = 0
      BMdat.norm[mi,cbm] = pvs
    }
  }
  
  return (BMdat.norm)
  
}
  

# This function looks will use pair values found in both BMdat and BMdat.norm to build empirical extrapolation functions for each slide or slide-position
# These functions will then be used to infer values for missing pairs (e.g. BMdat has a value but it is missing in BMdat.norm)
extrapolateNormalization = function(dat0, BMdat, BMdat.norm, DISPLAY_FITS = FALSE, polynomial_degree = 2) {
  
  dbr = cbind(dat0[,c("slide","position", "Cell.ID")],BMdat)
  dbn = cbind(dat0[,c("slide","position", "Cell.ID")],BMdat.norm)
  
  # aggregate(x = dbr[,bm], by = list(dbr[,"slide"]), FUN = "medianRemoveNAs")
  # aggregate(x = dbn[,bm], by = list(dbn[,"slide"]), FUN = "medianRemoveNAs")

  dbdif = dbn[,c("slide","position", "Cell.ID")]
  for(cbm in bm) {
    
    rawCi = paste0("Raw.",gsub("Mean.AF000.","",cbm))
    normCi = paste0("Norm.",gsub("Mean.AF000.","",cbm))
      
    dbdif[,paste0("NormMinusRaw.",gsub("Mean.AF000.","",cbm))] = dbn[,cbm] - dbr[,cbm]
    dbdif[,rawCi] = dbr[,cbm]
    dbdif[,normCi] = dbn[,cbm]
    
    
    pd1 = na.omit(dbdif[,c(rawCi,normCi)])
    fxlim=c(0,max(pd1[,1]))
    fylim=c(0,max(pd1[,2]))
    dslides = unique(dbdif[,"slide"])
    for (slide in dslides) {
      d = dbdif
      d = d[d[,"slide"] == slide, ]
      pd1 = na.omit(d[,c(rawCi,normCi)])
      # plot(pd1, main=slide, xlim=fxlim,ylim=fylim)
      colnames(pd1) = c("x","y")
      
      fit <- calcIntensityCorrectionFit(pd1, rawCi, normCi, fxlim, fylim, slide, DISPLAY_FITS, polynomial_degree)

      # for each entry that is NA in BMdat.norm and !NA in BMdat for a given slide fill in the missing value using prediction
      mi = which(is.na(BMdat.norm[,cbm]) & !is.na(BMdat[,cbm]) & dat0[,"slide"] == slide) 
      if (length(mi) > 0) {
        pvs = predict(fit, newdata = data.frame(x = BMdat[mi,cbm]))
        pvs[!is.na(pvs) & pvs < 0] = 0
        BMdat.norm[mi,cbm] = pvs
      }
    }
  }

  
  
  return (BMdat.norm)
}

#######################################################################################
#######################################################################################

getEvaluationList = function(elistNum=1) {
  
  #######################################################################################
  ##### fill out parameters with default values
  #######################################################################################
  
  elist = list()
  
  # use to toggle specific evaluation runs on and off for pipeline development purposes
  elist[["PARM_RUN"]] = TRUE
  
  elist[["PLOT_TO_JPEG"]] = FALSE
  elist[["PLOT_EVAL_DATA"]] = FALSE
  elist[["SAVE_EVAL_DATA"]] = FALSE
  elist[["SAVE_EVAL_DATA_FILENAME"]] = "EVALUTION_DATA"
  
  elist[["GEN_QQ_PLOT"]] = FALSE
  elist[["GEN_DOT_NUM_SLIDES_PLOT"]] = FALSE
  elist[["GEN_HIST_PLOT"]] = TRUE
  elist[["GEN_SEG_OBJ_MAP"]] = FALSE
  elist[["GEN_SEG_OBJ_INT_PLOT"]] = TRUE
  
  
  
  # "NORMALIZE_ALL_FOVS_EVALUATE_ALL_FOVS", "NORMALIZE_1_FOVS_EVALUATE_ALL_FOVS", "NORMALIZE_2_FOVS_EVALUATE_ALL_FOVS", "NORMALIZE_5_FOVS_EVALUATE_ALL_FOVS"
  #  "NORMALIZE_100_OBJS_EVALUATE_ALL_FOVS", 
  elist[["PARM_NORMALIZE_EVALUATE_PIPELINE"]] = "NORMALIZE_ALL_FOVS_EVALUATE_ALL_FOVS"
  
  # PARM_NORM_METHODS = c("MEDIAN","QUANTILE")
  elist[["PARM_NORM_METHOD"]] = "MEDIAN"
  elist[["PARM_SCALE_NORM_TO_RAW_MEDIAN"]] = TRUE
  elist[["PARM_UNEVEN_MATRIX_LENGTH_METHOD"]] = "SAMPLING_TO_MIN_LENGTH"
  
  
  # elist[["PARM_SEG_OBJ_NAMES"]] = c("Grid512","Grid256","Grid128","Grid64","NucleiSCA")
  elist[["PARM_SEG_OBJ_NAME"]] = "Grid256"    
  
  #elist[["PARM_EVALUATE_SLIDES"]] = c("S18030274","S18030275", "S18030276")
  elist[["PARM_EVALUATE_SLIDE"]] = "S18030274"
  elist[["PARM_EVALUATE_POSITION"]] = NULL #45
  elist[["PARM_FILTER_MIN_DAPI_RNDtoRND_CORR"]] = 75
  elist[["PARM_FILTER_MIN_DAPI_INTENSITY"]] = 1
  elist[["PARM_FILTER_MIN_OBJ_AREA"]] = 0
  elist[["PARM_FILTER_MAX_OBJ_AREA"]] = 1e38
  
  elist[["PARM_NORM_VE_MIN_CORR"]] = 0
  
  elist[["PARM_EVAL_SEG_OBJ"]] = FALSE
  elist[["PARM_EVAL_SEG_OBJ_NAME"]] = "NucleiSCA"
  
  elist[["PARM_EVAL_MIN_DAPI_RNDtoRND_CORR"]] = 90
  elist[["PARM_EVAL_MIN_DAPI_INTENSITY"]] = 1
  elist[["PARM_EVAL_MIN_OBJ_AREA"]] = 85
  elist[["PARM_EVAL_MAX_OBJ_AREA"]] = 471
  
  # following contains data with three exposure times of 20, 50, and 100
  elist[["PARM_EVAL_EXPOSURE_TIME"]] = NA
  elist[["PARM_EVAL_MIN_NUM_DAPI_RNDS"]] = 14
  elist[["PARM_EVAL_DAPI_RNDS"]] = c("RG003", "RG005", "RG007", "RG008", "RG010", "RG012", "RG014", "RG015", "RG017", "RG018", "RG020", "RG021", "RG022", "RG024")
  

  elist[["PARM_EVAL_QUANTILE_PROBS"]] = c(0, 0.25,0.5,0.75, 1)
  elist[["PARM_PLOT_GOOD_GRIDS"]] = FALSE
  elist[["PARM_PLOT_POSITION"]] = "45"
  
  #positions = c(NA,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,32,33,34,35,36,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,67,69,70,71,73,74,75,76,77,78,79,81,82,83,84)
  
  #######################################################################################
  ##### customize parameters for each specfic test scenario
  #######################################################################################
  
  cns = colnames(eRuns)
  for(cn in cns) {
    elist[[cn]] =   eRuns[elistNum,cn]
  }
  
  if (!is.na(elist[["PARM_EVALUATE_POSITION"]])) {
    if (!is.na(elist[["PARM_EVALUATE_POSITION"]]) & (length(elist[["PARM_EVALUATE_POSITION"]]) == 1  & length(grep("_",elist[["PARM_EVALUATE_POSITION"]])) > 0)) {
      elist[["PARM_EVALUATE_POSITION"]] = unlist(strsplit(elist[["PARM_EVALUATE_POSITION"]], "_"))
    }
    elist[["PARM_PLOT_POSITION"]] = elist[["PARM_EVALUATE_POSITION"]][1]
  }
  
  
  if (!is.na(elist[["PARM_EVAL_EXPOSURE_TIME"]])) {
    if (elist[["PARM_EVAL_EXPOSURE_TIME"]] == 20) {
      elist[["PARM_EVAL_MIN_NUM_DAPI_RNDS"]] = 2
      elist[["PARM_EVAL_DAPI_RNDS"]] = c("RG003", "RG005")
    } else if (elist[["PARM_EVAL_EXPOSURE_TIME"]] == 50) {
      elist[["PARM_EVAL_MIN_NUM_DAPI_RNDS"]] = 2
      elist[["PARM_EVAL_DAPI_RNDS"]] = c("RG007", "RG024")
    } else if (elist[["PARM_EVAL_EXPOSURE_TIME"]] == 100) {
      elist[["PARM_EVAL_MIN_NUM_DAPI_RNDS"]] = 10
      elist[["PARM_EVAL_DAPI_RNDS"]] = c("RG008", "RG010", "RG012", "RG014", "RG015", "RG017", "RG018", "RG020", "RG021", "RG022")
    } else if (elist[["PARM_EVAL_EXPOSURE_TIME"]] == "100_2A") {
      elist[["PARM_EVAL_MIN_NUM_DAPI_RNDS"]] = 2
      elist[["PARM_EVAL_DAPI_RNDS"]] = c("RG008", "RG010")
    } else if (elist[["PARM_EVAL_EXPOSURE_TIME"]] == "100_2B") {
      elist[["PARM_EVAL_MIN_NUM_DAPI_RNDS"]] = 2
      elist[["PARM_EVAL_DAPI_RNDS"]] = c("RG012", "RG014")
    } else if (elist[["PARM_EVAL_EXPOSURE_TIME"]] == "100_2C") {
      elist[["PARM_EVAL_MIN_NUM_DAPI_RNDS"]] = 2
      elist[["PARM_EVAL_DAPI_RNDS"]] = c("RG015", "RG017")
    } else if (elist[["PARM_EVAL_EXPOSURE_TIME"]] == "100_2D") {
      elist[["PARM_EVAL_MIN_NUM_DAPI_RNDS"]] = 2
      elist[["PARM_EVAL_DAPI_RNDS"]] = c("RG018", "RG020")
    } else if (elist[["PARM_EVAL_EXPOSURE_TIME"]] == "100_2E") {
      elist[["PARM_EVAL_MIN_NUM_DAPI_RNDS"]] = 2
      elist[["PARM_EVAL_DAPI_RNDS"]] = c("RG021", "RG022")
    }
  }
  
  if (length(elist[["PARM_EVAL_DAPI_RNDS"]]) == 1  & length(grep(":",elist[["PARM_EVAL_DAPI_RNDS"]])) > 0) {
    elist[["PARM_EVAL_DAPI_RNDS"]] = unlist(strsplit(elist[["PARM_EVAL_DAPI_RNDS"]], ":"))
  }

  
  if (!"PARM_SEG_OBJ_FILE_NAME" %in% names(elist)) {
    elist[["PARM_SEG_OBJ_FILE_NAME"]] = paste0("FOLFOX_DAPI_", elist[["PARM_SEG_OBJ_NAME"]] , ".RData")
  }
  
  if (!"PARM_EVAL_SEG_OBJ_FILE_NAME" %in% names(elist)) {
    elist[["PARM_EVAL_SEG_OBJ_FILE_NAME"]] = paste0("FOLFOX_DAPI_", elist[["PARM_EVAL_SEG_OBJ_NAME"]] , ".RData")
  }
  
  return (elist)
}


getNumEvaluationLists = function() {
  return (nrow(eRuns))
}

##############################################################
##############################################################

# This function looks will use pair values found in both BMdat and BMdat.norm to build empirical extrapolation functions for each slide or slide-position
# These functions will then be used to infer values for missing pairs (e.g. BMdat has a value but it is missing in BMdat.norm)
genImageNormFunc = function(dat0, BMdat, BMdat.norm, cbm, DISPLAY_FITS = FALSE, polynomial_degree = 2) {
  
  dbr = cbind(dat0[,c("slide","position", "Cell.ID")],BMdat)
  dbn = cbind(dat0[,c("slide","position", "Cell.ID")],BMdat.norm)
  
  dbdif = dbn[,c("slide","position", "Cell.ID")]
  
  
  rawCi = paste0("Raw.",gsub("Mean.AF000.","",cbm))
  normCi = paste0("Norm.",gsub("Mean.AF000.","",cbm))
  
  dbdif[,paste0("NormMinusRaw.",gsub("Mean.AF000.","",cbm))] = dbn[,cbm] - dbr[,cbm]
  dbdif[,rawCi] = dbr[,cbm]
  dbdif[,normCi] = dbn[,cbm]
  
  SlideNormFuncList = list()
  pd1 = na.omit(dbdif[,c(rawCi,normCi)])
  fxlim=c(0,max(pd1[,1]))
  fylim=c(0,max(pd1[,2]))
  
  dslides = unique(dbdif[,"slide"])
  
  for (slide in dslides) {
    d = dbdif
    d = d[d[,"slide"] == slide, ]
    pd1 = na.omit(d[,c(rawCi,normCi)])
    # plot(pd1, main=slide, xlim=fxlim,ylim=fylim)
    colnames(pd1) = c("x","y")
    
    fit <- calcIntensityCorrectionFit(pd1, rawCi, normCi, fxlim, fylim, slide, DISPLAY_FITS, polynomial_degree)
    
    SlideNormFuncList[[slide]] = fit
  }
  
  
  return (SlideNormFuncList)
}

##############################################################
##############################################################

# This function looks will use pair values found in both BMdat and BMdat.norm to build empirical extrapolation functions for each slide or slide-position
# These functions will then be used to infer values for missing pairs (e.g. BMdat has a value but it is missing in BMdat.norm)
applyImageNorm = function(BMdat.norm, cbm, SlideNormFuncList) {
  
  dslides = names(SlideNormFuncList)
  
  for (slide in dslides) {
    fit =SlideNormFuncList[[slide]]
    mi = which(!is.na(BMdat.norm[,cbm]) & BMdat.norm[,"slide"] == slide) 
    if (length(mi) > 0) {
      pvs = predict(fit, newdata = data.frame(x = BMdat.norm[mi,cbm]))
      pvs[!is.na(pvs) & pvs < 0] = 0
      BMdat.norm[mi,cbm] = pvs
    } else {
      cat("Warning mi is zero length\n")
    }
  }
  
  return (BMdat.norm)
  
}

##############################################################
##############################################################

applyGridBased_FLINO = function(cellObjectDataFile = file.path(pathToData, "FOLFOX_MarkerData_NucleiSCA.RData"),
                       gridObjectDataFile = file.path(pathToData,  "FOLFOX_MarkerData_Grid32.RData"),
                       slides = c("S18030274", "S18030275", "S18030276")) {
  
  cat("Performing Grid-Based Object Normalization","\n")
  
  # *********************************************************************************************
  # Load Grid32  Object Marker Data
  cat("Loading...","\n")
  cat(gridObjectDataFile,"\n")
  load(gridObjectDataFile, verbose = FALSE)
  
  # No filtering used
  dat.Corr_SEG = data.frame()
  
  bms = colnames(dat.marker_SEG)
  bms = bms[startsWith(bms,"Mean.")]
  dat0_AllSlides = get_dat0(dat.marker_SEG, dat.Corr_SEG, cy3rounds, cy5rounds, bms)
  dat0_AllSlides = dat0_AllSlides[dat0_AllSlides[,"slide"] %in% slides, ]
  cat("Using all TMA positions to perform Grid Object Normalization","\n")
  #positions = c(81,82,83,84)
  #dat0_AllSlides = dat0_AllSlides[dat0_AllSlides[,"position"] %in% positions, ]
  # *********************************************************************************************
  
  BMdat = dat0_AllSlides
  # Transform  to log space
  for(cbm in bms) {
    BMdat[,cbm] = log2(1 + BMdat[,cbm])
  }
  
  
  quantileProb = 0.75 # NORMALIZING GRIDS OBJECTS - USE 75% Quantile
  slideQuantile = aggregate(BMdat[,bms],by=list(BMdat[,"slide"]),FUN= function(x){quantile(x[x>0], quantileProb, na.rm = TRUE)})
  colnames(slideQuantile)[1] = "slide"
  BMdat.norm = BMdat
  for(cbm in bms) {
    for(slide in unique(BMdat.norm[,"slide"])) {
      ii = which(BMdat.norm[,"slide"] == slide)
      BMdat.norm[ii, cbm] <- BMdat.norm[ii, cbm] / slideQuantile[slideQuantile[,"slide"]==slide,cbm] 
    }
  }
  rm(slideQuantile)

  
  
  for(cbm in bms) {
    nf = median(BMdat[,cbm],na.rm = TRUE) / median(BMdat.norm[,cbm],na.rm = TRUE)
    BMdat.norm[,cbm] = BMdat.norm[,cbm] * nf
  }
  
  
  grid.fit.bms = bms
  grid.SlideNormFuncList_LogSpace = list()
  FIT_POLYNOMIAL_DEGREE = 2
  showDisplayFits = FALSE
  for(cbm in grid.fit.bms) {
    grid.SlideNormFuncList_LogSpace[[cbm]] = genImageNormFunc(BMdat[,c("slide","position", "Cell.ID")], BMdat[,-c(1:3)], BMdat.norm[,-c(1:3)], cbm, DISPLAY_FITS = showDisplayFits, polynomial_degree = FIT_POLYNOMIAL_DEGREE)
  }
  
  # Transform back to absolute space
  for(cbm in bms) {
    BMdat[,cbm] = (2^BMdat[,cbm] - 1)
    BMdat.norm[,cbm] = (2^BMdat.norm[,cbm] - 1)
  }
  
  aggDat = BMdat
  aggDat = aggregate(aggDat[,bms], by = list(aggDat[,"slide"], aggDat[,"position"]), FUN = "medianRemoveNAs")
  colnames(aggDat)[1] = "slide"
  colnames(aggDat)[2] = "position"
  BMdat.FOV = aggDat
  
  aggDat = BMdat.norm
  aggDat = aggregate(aggDat[,bms], by = list(aggDat[,"slide"], aggDat[,"position"]), FUN = "medianRemoveNAs")
  colnames(aggDat)[1] = "slide"
  colnames(aggDat)[2] = "position"
  BMdat.norm.FOV = aggDat
  
  grid.BMdat.FOV = BMdat.FOV
  grid.BMdat.norm.FOV = BMdat.norm.FOV
  grid.bms = bms
  grid.BMdat = BMdat
  grid.BMdat.norm = BMdat.norm
  
  
  cat("Applying Grid Object Normalization to Segmented Cell Objects","\n")
  
  # *********************************************************************************************
  # Load Segmented Cell Object Marker Data
  cat("Loading...","\n")
  cat(cellObjectDataFile,"\n")
  load(cellObjectDataFile, verbose = FALSE)

  # No filtering used
  dat.Corr_SEG = data.frame()
  
  bms = colnames(dat.marker_SEG)
  bms = bms[startsWith(bms,"Mean.")]
  dat0_AllSlides = get_dat0(dat.marker_SEG, dat.Corr_SEG, cy3rounds, cy5rounds, bms)
  dat0_AllSlides = dat0_AllSlides[dat0_AllSlides[,"slide"] %in% slides, ]
  # positions = c(81,82,83,84)
  # dat0_AllSlides = dat0_AllSlides[dat0_AllSlides[,"position"] %in% positions, ]
  # *********************************************************************************************
  
  BMdat = dat0_AllSlides
  # Transform  to log space
  for(cbm in bms) {
    BMdat[,cbm] = log2(1 + BMdat[,cbm])
  }
  
  nuclei.BMdat.gridNorm = BMdat
  for(cbm in grid.fit.bms) {
    nuclei.BMdat.gridNorm = applyImageNorm(nuclei.BMdat.gridNorm, cbm, grid.SlideNormFuncList_LogSpace[[cbm]])
  }
  
  for(cbm in grid.fit.bms) {
    nf = median(BMdat[,cbm],na.rm = TRUE) / median(nuclei.BMdat.gridNorm[,cbm],na.rm = TRUE)
    nuclei.BMdat.gridNorm[,cbm] = nuclei.BMdat.gridNorm[,cbm] * nf
  }
  
  
  # Transform back to absolute space
  for(cbm in bms) {
    BMdat[,cbm] = (2^BMdat[,cbm] - 1)
    BMdat.norm[,cbm] = (2^BMdat.norm[,cbm] - 1)
    nuclei.BMdat.gridNorm[,cbm] = (2^nuclei.BMdat.gridNorm[,cbm] - 1)
  }
  aggDat = BMdat
  aggDat = aggregate(aggDat[,bms], by = list(aggDat[,"slide"], aggDat[,"position"]), FUN = "medianRemoveNAs")
  colnames(aggDat)[1] = "slide"
  colnames(aggDat)[2] = "position"
  BMdat.FOV = aggDat
  
  aggDat = BMdat.norm
  aggDat = aggregate(aggDat[,bms], by = list(aggDat[,"slide"], aggDat[,"position"]), FUN = "medianRemoveNAs")
  colnames(aggDat)[1] = "slide"
  colnames(aggDat)[2] = "position"
  BMdat.norm.FOV = aggDat
  
  aggDat = nuclei.BMdat.gridNorm
  aggDat = aggregate(aggDat[,bms], by = list(aggDat[,"slide"], aggDat[,"position"]), FUN = "medianRemoveNAs")
  colnames(aggDat)[1] = "slide"
  colnames(aggDat)[2] = "position"
  nuclei.BMdat.gridNorm.FOV = aggDat
  
  nuclei.BMdat.FOV = BMdat.FOV
  nuclei.BMdat = BMdat
  
  
  nuclei.BMdat <<- nuclei.BMdat
  nuclei.BMdat.gridNorm <<- nuclei.BMdat.gridNorm
  nuclei.BMdat.FOV <<- nuclei.BMdat.FOV
  nuclei.BMdat.gridNorm.FOV <<- nuclei.BMdat.gridNorm.FOV

}

##############################################################
##############################################################


