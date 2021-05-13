while(dev.next()>1) dev.off()
options(stringsAsFactors=FALSE)
rm(list=objects())

library("stringr")
library("png")

pathToData = "Data"
pathToRcode= "Rcode" 
pathToFigures = "figures"
pathToImages = "images"

source(file.path(pathToRcode,"FLINO_Functions.R"))


SCREEN = "SCREEN"
TIFF = "TIFF"
EPS = "EPS"

PLOT_TO = TIFF

GENERATE_PANEL_1 = FALSE  # NOT IMPLEMENTED : Overview of the MxIF data processing workflows for normalizing images and intensities of biological features across virtual slides
GENERATE_PANEL_2 = TRUE   # Performance of normalization methods for DAPI segmented nuclei objects
GENERATE_PANEL_3 = TRUE   # The effect of control sample number on error correction of TMA slide images
GENERATE_PANEL_4 = TRUE  # NOT IMPLEMENTED : Application of grid-object normalization to BAX staining images of four cell lines across three physical slides. 

GENERATE_PANEL_S1 = TRUE
GENERATE_PANEL_S2 = TRUE 
GENERATE_PANEL_S3 = TRUE   # Supporting figure for main text figure 2
GENERATE_PANEL_S4 = TRUE   # Impact of filtering of objects for DAPI round-to-round correlation on performance of normalization
GENERATE_PANEL_S5 = TRUE   # Performance of normalization methods utilizing a grid-based approach
GENERATE_PANEL_S6a = TRUE
GENERATE_PANEL_S6b = TRUE
GENERATE_PANEL_S6c = TRUE
GENERATE_PANEL_S6d = TRUE
GENERATE_PANEL_S6e = TRUE


CV_METRIC_TYPE = "CV_OF_EVALUTATION_OBJECTS"
#CV_METRIC_TYPE = "CV_OF_NORMALIZATION_OBJECTS"

if (CV_METRIC_TYPE == "CV_OF_NORMALIZATION_OBJECTS") {
  CV_METRIC_RawAbsSpace = "SegObjErrCV_RawAbsSpace"
  CV_METRIC_NormAbsSpace = "SegObjErrCV_NormAbsSpace"
  CV_METRIC_RawLogSpace = "SegObjErrCV_RawLogSpace"
  CV_METRIC_NormLogSpace = "SegObjErrCV_NormLogSpace"
} else if (CV_METRIC_TYPE == "CV_OF_EVALUTATION_OBJECTS") {
  CV_METRIC_RawAbsSpace = "EvalSegObj_SegObjErrCV_RawAbsSpace"
  CV_METRIC_NormAbsSpace = "EvalSegObj_SegObjErrCV_NormAbsSpace"
  CV_METRIC_RawLogSpace = "EvalSegObj_SegObjErrCV_RawLogSpace"
  CV_METRIC_NormLogSpace = "EvalSegObj_SegObjErrCV_NormLogSpace"
}

get_display_PARM_METRIC = function(metric) {
  if (metric == CV_METRIC_NormLogSpace)
    return ("log")
  else if (metric == CV_METRIC_NormAbsSpace)
    return ("absolute")
  else
    return ("?")
  
}



#### Can build data file used to generate all of the plots but it will take 
#### a very very long time (1 month) without running in a HPC environment.
#### There are over 26,870 evaluation runs. 
if (!file.exists(file.path(pathToData,"results_eRuns_FLINO.txt"))) {
  analysisRun = "eRuns_FLINO_Manuscript.txt"
  source(file.path(pathToRcode,"Run_FLINO_Evaluator.R"))
  datFLINOresults = results
} else {
  #### Load results file from 26,870 evaluation runs. 
  datFLINOresults = read.delim(file=file.path(pathToData, "results_eRuns_FLINO.txt"), header = TRUE, sep = "\t", blank.lines.skip = TRUE, na.strings = c("NA","NaN"))
}



#### Build file used to build figure 2 if data does not already exist. 
#### This could take 10 minutes of CPU time
if (!file.exists(file.path(pathToData,"NucleiSCA_Q50NZ_14VS.RData"))) {
  
  analysisRun = "eRuns_NucleiSCA_Q50NZ_14VS.txt"
  source(file.path(pathToRcode,"Run_FLINO_Evaluator.R"))
  
  #pd = cbind(dat0[,c("slide","position")],BMdat[,"Mean.AF000.dapi.dapi"],BMdat.norm[,"Mean.AF000.dapi.dapi_NormInLogSpace"])
  pd = cbind(dat0[,c("slide","position")],EvalSegObj_BMdat[,"Mean.AF000.dapi.dapi"],EvalSegObj_BMdat.norm[,"Mean.AF000.dapi.dapi_NormInLogSpace"])
  colnames(pd) = c("slide","position","Mean.AF000.dapi.dapi","Mean.AF000.dapi.dapi_NormInLogSpace")
  pd = aggregate(pd[,c("Mean.AF000.dapi.dapi","Mean.AF000.dapi.dapi_NormInLogSpace")],by=list(pd[,"slide"], pd[,"position"]), medianRemoveNAs)
  colnames(pd) = c("slide","position","RawIntensity","NormIntensity")
  medianRawIntensities = medianRemoveNAs(EvalSegObj_BMdat[,"Mean.AF000.dapi.dapi"])
  medianNormIntensities = medianRemoveNAs(EvalSegObj_BMdat.norm[,"Mean.AF000.dapi.dapi_NormInLogSpace"])
  slides = unique(pd[,"slide"])
  
  save(results, pd, slides, medianRawIntensities, medianNormIntensities, file=file.path(pathToData,"NucleiSCA_Q50NZ_14VS.RData"))
  
} else {
  load(file.path(pathToData,"NucleiSCA_Q50NZ_14VS.RData"))
}

getSlideNames = function(slides) {
  slides = gsub("S18030274","Slide A1",slides)
  slides = gsub("S18030275","Slide A2",slides)
  slides = gsub("S18030276","Slide A3",slides)
  return (slides)
}


plotVirtualSlideData = function(pdResults, pdTrials, pdPlotColumns, plotBottom = 0, plotTop = 1.2, XLAB = "Slide-to-Slide Staining Intensity Corretion Method") {
  # plotBottom = 0
  # plotTop = 0.8
  
  gList = list()
  tabbedMeans = data.frame()
  tabbedSD = data.frame()
  tabbedCount = data.frame()
  tabbedMedians = data.frame()
  tabbedMin = data.frame()
  tabbedMax = data.frame()
  tabbed25 = data.frame()
  tabbed75 = data.frame()
  for(gi in 1:nrow(pdPlotColumns)) {
    
    pd = pdResults
    pd = pd[pd[,"PARM_EVALUATE_SLIDE"] == "S18030274", ]
    pd = pd[pd[,"PARM_NORM_METHOD"] == pdPlotColumns[gi, "PARM_NORM_METHOD"], ]
    
    if (!is.na(pdPlotColumns[gi,"PARM_EVAL_SEG_OBJ"])) {
      pd =  pd[pd[,"PARM_EVAL_SEG_OBJ"] == pdPlotColumns[gi,"PARM_EVAL_SEG_OBJ"], ] 
    }
    pd =  pd[pd[,"PARM_NORMALIZE_EVALUATE_PIPELINE"] == pdPlotColumns[gi,"PARM_NORMALIZE_EVALUATE_PIPELINE"], ]
    if (!is.na(pdPlotColumns[gi,"PARM_FILTER_MIN_DAPI_RNDtoRND_CORR"])) {
      pd =  pd[pd[,"PARM_FILTER_MIN_DAPI_RNDtoRND_CORR"] == pdPlotColumns[gi,"PARM_FILTER_MIN_DAPI_RNDtoRND_CORR"], ]
    }
    pd =  pd[pd[,"PARM_SEG_OBJ_NAME"] == pdPlotColumns[gi,"PARM_SEG_OBJ_NAME"], ]  
    
    if (!is.na(pdPlotColumns[gi,"PARM_EVAL_SEG_OBJ_NAME"])) {
      pd =  pd[pd[,"PARM_EVAL_SEG_OBJ_NAME"] == pdPlotColumns[gi,"PARM_EVAL_SEG_OBJ_NAME"], ]   
    }
    
    gList[[gi]] =  pd[, pdPlotColumns[gi,"PARM_METRIC"]]
    
    if (gi == 1) {
      tabbedMeans = data.frame(mean(gList[[gi]]))
      tabbedSD = data.frame(sd(gList[[gi]]))
      tabbedCount = data.frame(length(gList[[gi]]))
      tabbedMedians = data.frame(median(gList[[gi]]))
      
      tabbedMin = data.frame(as.numeric(quantile(gList[[gi]],prob=0)[1]))
      tabbedMax = data.frame(as.numeric(quantile(gList[[gi]],prob=1)[1]))
      tabbed25 = data.frame(as.numeric(quantile(gList[[gi]],prob=0.25)[1]))
      tabbed75 = data.frame(as.numeric(quantile(gList[[gi]],prob=0.75)[1]))
    } else {
      tabbedMeans = cbind(tabbedMeans, data.frame(mean(gList[[gi]])))
      tabbedSD = cbind(tabbedSD, data.frame(sd(gList[[gi]])))
      tabbedCount = cbind(tabbedCount, data.frame(length(gList[[gi]])))
      tabbedMedians = cbind(tabbedMedians, data.frame(median(gList[[gi]])))
      
      tabbedMin = cbind(tabbedMin, data.frame(as.numeric(quantile(gList[[gi]],prob=0)[1])))
      tabbedMax = cbind(tabbedMax, data.frame(as.numeric(quantile(gList[[gi]],prob=1)[1])))
      tabbed25 = cbind(tabbed25, data.frame(as.numeric(quantile(gList[[gi]],prob=0.25)[1])))
      tabbed75 = cbind(tabbed75, data.frame(as.numeric(quantile(gList[[gi]],prob=0.75)[1])))
    }
  }
  
  tabbedMeans= as.matrix(tabbedMeans)
  tabbedSD= as.matrix(tabbedSD)
  tabbedCount= as.matrix(tabbedCount)
  tabbedSE = tabbedSD /  sqrt(tabbedCount)
  tabbedMedians = as.matrix(tabbedMedians)
  
  tabbedMin = as.matrix(tabbedMin)
  tabbedMax = as.matrix(tabbedMax)
  tabbed25 = as.matrix(tabbed25)
  tabbed75 = as.matrix(tabbed75)
  
  if (MEDIAN_BAR_VERSION) {
    barValues = tabbedMedians
    quantileBar = tabbedMeans
  } else {
    barValues = tabbedMeans
    quantileBar = tabbedMedians
  }

  barCenters <<- barplot(height = barValues,
                         names.arg = rep("",times=ncol(barValues)),
                         beside = TRUE, las = 1,
                         ylim = c(plotBottom, plotTop),
                         cex.names = acex,
                         cex.axis = acex,
                         cex.lab = acex,
                         main = "",
                         ylab = "",
                         xlab = "",
                         col="grey",
                         border = "black", axes = TRUE,
                         legend.text = SHOW_LEGEND,
                         args.legend = list(title = "Conditions", 
                                            x = "topright",
                                            cex = cexLegend))
 
  if (DISPLY_GRID) {
    grid(nx=NA, ny=NULL)
  }
  
  if(ADD_HORIZONTAL_BAR_AT_GIVEN_VALUE) {
    abline(h=HORIZONTAL_BAR_VALUE,lty=2,lwd=1,col="red")
  } else if(ADD_HORIZONTAL_MIN_BAR_VALUE) {
    abline(h=min(barValues,na.rm = TRUE),lty=2,lwd=1,col="red")
  }
  

  if (DRAW_BOX_FRAME) {
    box()
  }
  

  
  barCenters <<- barplot(height = barValues,add=TRUE,
                         names.arg = rep("",times=ncol(barValues)),
                         beside = TRUE, las = 1,
                         ylim = c(plotBottom, plotTop),
                         cex.names = acex,
                         cex.axis = acex,
                         cex.lab = acex,
                         main = "",
                         ylab = "",
                         xlab = "",
                         col="grey",
                         border = "black", axes = TRUE,
                         legend.text = SHOW_LEGEND,
                         args.legend = list(title = "Conditions", 
                                            x = "topright",
                                            cex = cexLegend))
  
  
  #(1=bottom, 2=left, 3=top, 4=right).
  
  if (DISPLAY_DATA_POINTS) {
    
    for(gi in 1:length(gList)) {
      
      if (length(gList[[gi]]) > 0) {
        stripchart(gList[[gi]],at = barCenters[gi], vertical = TRUE, method = "jitter",
                   jitter = symJitter,
                   pch = sympch, col = symcol,
                   add = TRUE,cex=symSize)
      }
    }
    
  }
  
  
  if (DISPLAY_ERROR_BARS) {
    
    if (USE_STANDARD_DEVIATION_FOR_ERROR_BARS) {
      tabbedError = tabbedSD
    } else {
      tabbedError = 2 * tabbedSE
    }
    
    
    mylwd = 2
    barCol = "black"
    
    segments(barCenters, tabbedMeans - tabbedError, barCenters,
             tabbedMeans + tabbedError, lwd = mylwd, col = barCol)
    
    arrows(barCenters, tabbedMeans - tabbedError, barCenters,
           tabbedMeans + tabbedError, lwd = mylwd, angle = 90,
           code = 3, length = 0.05 * 1, col = barCol)
    
  }
  
  if (DISPLAY_QUANTILE_BARS) {
    
    barCol = "black"
    mylwd = 1
    widthMod = 0.5 * (barCenters[2] - barCenters[1]) / 2
    for(gi in 1:length(gList)) {
      segments(barCenters[gi]-widthMod/2,tabbedMin[gi],barCenters[gi]+widthMod/2,tabbedMin[gi],lwd = mylwd, col = barCol)
      segments(barCenters[gi]-widthMod/2,tabbedMax[gi],barCenters[gi]+widthMod/2,tabbedMax[gi],lwd = mylwd, col = barCol)
    }
    segments(barCenters, tabbedMin, barCenters, tabbedMax, lwd = mylwd, col = barCol)

    mylwd = 2
    widthMod = 0.5 * (barCenters[2] - barCenters[1]) / 2
    for(gi in 1:length(gList)) {
      segments(barCenters[gi]-widthMod/2,quantileBar[gi],barCenters[gi]+widthMod/2,quantileBar[gi],lwd = mylwd, col = barCol)
    }

  }
  
  
  for(gi in 1:length(gList)) {
    mtext( pdPlotColumns[gi, "Display_X_Line_0"], line=0, side =1,at=barCenters[gi], cex=t1bcex)
    mtext( pdPlotColumns[gi, "Display_X_Line_1"], line=1, side =1,at=barCenters[gi], cex=t1bcex)  
    mtext( pdPlotColumns[gi, "Display_X_Line_2"], line=2, side =1,at=barCenters[gi], cex=t1bcex)  
  }
  mtext(XLAB, line=3.2, side =1,at=mean(barCenters), cex=t1cex)
  
  
  if (length(pdTrials) == 1) {
    mtext(pdTrials, line=1, side =3, cex=t3cex)
  }
  
  
  
  pdTable = rbind(tabbedMeans,tabbedSD,tabbedCount,tabbedSE,tabbedMedians,tabbedMin,tabbedMax,tabbed25, tabbed75)
  colnames(pdTable) = paste0("v",c(1:ncol(pdTable)))
  row.names(pdTable) = c("mean","sd","n","se","median","min","max","25_percent","75_percent")
  return(pdTable)
}




###############################################################
###############################################################

# NOT IMPLEMENTED

if (GENERATE_PANEL_1) {
  panelName = "Fig_1"
  
  if (PLOT_TO == TIFF) {
    #jpeg(file.path(pathToFigures,paste0("Dapi",".slide",".jpeg")), quality = 100, width = 720, height = 480)
    myFile = file.path(pathToFigures,paste0(panelName,".tif"))
    tiff(myFile, height = 3.75, width = 7.5, units = 'in', compression="lzw",  type = "windows", res = 300)
    #par(cex=1.3)
  } else if (PLOT_TO == EPS) {
    myFile = file.path(pathToFigures,paste0(panelName,".eps"))
    setEPS()
    postscript(myFile,  height = 3.75, width = 7.5)
  }
  
  
  if (PLOT_TO != SCREEN) {
    dev.off()
  }
}


###############################################################
###############################################################


if (GENERATE_PANEL_2) {
  panelName = "Fig_2"
  
  if (PLOT_TO == TIFF) {
    #jpeg(file.path(pathToFigures,paste0("Dapi",".slide",".jpeg")), quality = 100, width = 720, height = 480)
    myFile = file.path(pathToFigures,paste0(panelName,".tif"))
    tiff(myFile, height = 5.4, width = 7.5, units = 'in', compression="lzw",  type = "windows", res = 300)
    #par(cex=1.3)
  } else if (PLOT_TO == EPS) {
    myFile = file.path(pathToFigures,paste0(panelName,".eps"))
    setEPS()
    postscript(myFile,  height = 3.75, width = 7.5)
  }
  
  filterCorr = 0
  value_PARM_NORMALIZE_EVALUATE_PIPELINE  = "NORMALIZE_ALL_FOVS_EVALUATE_ALL_FOVS"
  value_PARM_METRIC = CV_METRIC_NormLogSpace
  NUM_GROUPS = 7
  pdPlotColumns = data.frame(
    PARM_NORM_METHOD=c("MEDIAN",
                       "MEDIAN",
                       "MEDIAN",
                       "Q50NZ",
                       "Q50NZ",
                       "MRN",
                       "TMM"),
    Display_X_Line_0=c("No",
                       "Median",
                       "Median",
                       "Q50",
                       "Q50",
                       "MRN",
                       "TMM"),
    PARM_EVAL_SEG_OBJ=c(rep(TRUE, times=NUM_GROUPS)),
    PARM_SEG_OBJ_NAME=c(rep("NucleiSCA", times=NUM_GROUPS)),
    PARM_EVAL_SEG_OBJ_NAME=c(rep("NucleiSCA", times=NUM_GROUPS)),
    PARM_FILTER_MIN_DAPI_RNDtoRND_CORR= c(rep(filterCorr, times=NUM_GROUPS)),
    PARM_NORMALIZE_EVALUATE_PIPELINE=c(rep(value_PARM_NORMALIZE_EVALUATE_PIPELINE,times=NUM_GROUPS)),
    PARM_METRIC=c(CV_METRIC_RawLogSpace,
                  CV_METRIC_NormAbsSpace,
                  CV_METRIC_NormLogSpace,
                  CV_METRIC_NormAbsSpace,
                  CV_METRIC_NormLogSpace,
                  CV_METRIC_NormLogSpace,
                  CV_METRIC_NormLogSpace),
    
    Display_X_Line_1=c("Correction",
                       "Absolute",
                       "Log",
                       "Absolute",
                       "Log",
                       "Log",
                       "Log"),
    
    Display_X_Line_2= c(rep("", times=NUM_GROUPS))) 
  
  maxVirtualSlides = c(2,3,10,14)
  pdTrials = unique(datFLINOresults[,"PARM_EVAL_DAPI_RNDS"])
  pdTrials = pdTrials[(str_count(pdTrials, pattern = ":") + 1) %in% maxVirtualSlides]
  pdResults = datFLINOresults
  pdResults = pdResults[is.na(pdResults[,"PARM_EVALUATE_POSITION"]), ]
  pdResults = pdResults[pdResults[,"PARM_EVAL_DAPI_RNDS"] %in% pdTrials,]
  
  

  #par(fig=c(0,0.333+xOverlay,0.5+2*yOverlay,1), new=FALSE)
  #c(x1, x2, y1, y2)
  # make labels and margins smaller
  #par(mai=c(0.1,0.1,0.2,0.1))
  par(mar = c(3, 5, 1, 1)) #  bottom, left, top and right  c(5.1, 4.1, 4.1, 2.1).
  
  setDefaultPlotParm()
  DISPLY_GRID = FALSE
  DRAW_BOX_FRAME = FALSE
  acex <<- 1
  tcex <<-  1
  lcex <<- 1
  t2cex <<- 1
  ploty = 1.3
  pys = c(1.25,0.3, 0.25, 0.25, rep(0.2,times=NUM_GROUPS-4))
  ADD_HORIZONTAL_BAR_AT_GIVEN_VALUE  = TRUE
  HORIZONTAL_BAR_VALUE = 0.0738
  pdTable = plotVirtualSlideData(pdResults, pdTrials, pdPlotColumns[1:NUM_GROUPS,], 0 , ploty, "")
  abline(h=0)
  mtext("Error across virtual slides (Coefficient of Variation)", line=3.5, side =2, cex=t2cex)
  mtext("for dapi segmented nuclear object intensities", line=2.5, side =2, cex=t2cex)
  for(gi in c(1:7)) {
    text(format(pdTable[5,gi],digits = 3),  x = barCenters[gi],y =pys[gi], cex = tcex, xpd=NA)
  }
  
  text("A",x=-1,y=1.3,cex=1.5, font=2, xpd=NA)
  text("B",x=(barCenters[1]+barCenters[2])/2,y=1.3,cex=1.5, font=2, xpd=NA)
  text("C",x=(barCenters[1]+barCenters[2])/2,y=0.9,cex=1.5, font=2, xpd=NA)
  
  par(new=TRUE)
  

  par(mar = c(20, 15, 1, 1)) #  bottom, left, top and right
  
  
  load(file.path(pathToData,"NucleiSCA_Q50NZ_14VS.RData"))
  
  plot(pd[,"position"],pd[,"RawIntensity"],yaxt='n',xaxt='n',ylab="",xlab="",pch=NA,xlim=c(0,85),ylim=c(0,10000))
  axis(2,las=2)
  axis(1,at=c(0,10,20,30,40,50,60,70,80),labels = c("","","","","","","","",""))
  for(slide in slides) {
    ii = which(pd[,"slide"] == slide)
    points(pd[ii,"position"],pd[ii,"RawIntensity"],type="l" )
  }
  abline(h=medianRawIntensities,col="blue", lty=5)
  mtext("Median Intensity",side=2,line=3.5)
  
  cv1 = format(results[,"EvalSegObj_SegObjErrCV_RawLogSpace"],digits=2)
  text(x= 80, y= 8600,paste0("Uncorrected (CV = ", cv1, ")\n "),adj=1,cex=1)
  
  par(new=TRUE)
  
  par(mar = c(13, 15, 8, 1)) #  bottom, left, top and right
  
  
  plot(pd[,"position"],pd[,"NormIntensity"],yaxt='n',xaxt='n',ylab="",xlab="",pch=NA,xlim=c(0,85),ylim=c(0,10000))
  axis(2,las=2)
  axis(1,at=c(0,10,20,30,40,50,60,70,80))
  for(slide in slides) {
    ii = which(pd[,"slide"] == slide)
    points(pd[ii,"position"],pd[ii,"NormIntensity"],type="l" )
  }
  abline(h=medianNormIntensities,col="red", lty=5)
  mtext("Median Intensity",side=2,line=3.5)
  mtext("TMA Sample Position",side=1,line=2.3)
  cv2 = format(results[,"EvalSegObj_SegObjErrCV_NormLogSpace"],digits=1)
  text(x= 80, y= 8600,paste0("Q50 Log Space Normalized\n14 Virtual Slides (CV = ", cv2, ")"),adj=1,cex=1)
  
  if (PLOT_TO != SCREEN) {
    dev.off()
  }
}


###############################################################
###############################################################


if (GENERATE_PANEL_S3) {
  panelName = "Fig_S3"
  
  if (PLOT_TO == TIFF) {
    #jpeg(file.path(pathToFigures,paste0("Dapi",".slide",".jpeg")), quality = 100, width = 720, height = 480)
    myFile = file.path(pathToFigures,paste0(panelName,".tif"))
    tiff(myFile, height = 5.4, width = 7.5, units = 'in', compression="lzw",  type = "windows", res = 300)
    #par(cex=1.3)
  } else if (PLOT_TO == EPS) {
    myFile = file.path(pathToFigures,paste0(panelName,".eps"))
    setEPS()
    postscript(myFile,  height = 3.75, width = 7.5)
  }
  
  filterCorr = 0
  value_PARM_NORMALIZE_EVALUATE_PIPELINE  = "NORMALIZE_ALL_FOVS_EVALUATE_ALL_FOVS"
  value_PARM_METRIC = CV_METRIC_NormLogSpace
  NUM_GROUPS = 11
  pdPlotColumns = data.frame(
    PARM_NORM_METHOD=c("MEDIAN",
                       "QUANTILE",
                       "QUANTILE",
                       "MEDIAN",
                       "MEDIAN",
                       "Q50NZ",
                       "Q50NZ",
                       "MRN",
                       "MRN",
                       "TMM",
                       "TMM"),
    Display_X_Line_0=c("No",
                       "squa",
                       "squa",
                       "Median",
                       "Median",
                       "Q50",
                       "Q50",
                       "MRN",
                       "MRN",
                       "TMM",
                       "TMM"),
    PARM_EVAL_SEG_OBJ=c(rep(TRUE, times=NUM_GROUPS)),
    PARM_SEG_OBJ_NAME=c(rep("NucleiSCA", times=NUM_GROUPS)),
    PARM_EVAL_SEG_OBJ_NAME=c(rep("NucleiSCA", times=NUM_GROUPS)),
    PARM_FILTER_MIN_DAPI_RNDtoRND_CORR= c(rep(filterCorr, times=NUM_GROUPS)),
    PARM_NORMALIZE_EVALUATE_PIPELINE=c(rep(value_PARM_NORMALIZE_EVALUATE_PIPELINE,times=NUM_GROUPS)),
    PARM_METRIC=c(CV_METRIC_RawLogSpace,
                  CV_METRIC_NormAbsSpace,
                  CV_METRIC_NormLogSpace,
                  CV_METRIC_NormAbsSpace,
                  CV_METRIC_NormLogSpace,
                  CV_METRIC_NormAbsSpace,
                  CV_METRIC_NormLogSpace,
                  CV_METRIC_NormAbsSpace,
                  CV_METRIC_NormLogSpace,
                  CV_METRIC_NormAbsSpace,
                  CV_METRIC_NormLogSpace),
    
    Display_X_Line_1=c("Correct",
                       " Absolute",
                       "Log",    
                       "Absolute",
                       "Log",
                       "Absolute",
                       "Log",
                       "Absolute",
                       "Log",
                       "Absolute",
                       "Log"),
    
    Display_X_Line_2= c(rep("", times=NUM_GROUPS))) 
  
  maxVirtualSlides = c(2,3,10,14)
  pdTrials = unique(datFLINOresults[,"PARM_EVAL_DAPI_RNDS"])
  pdTrials = pdTrials[(str_count(pdTrials, pattern = ":") + 1) %in% maxVirtualSlides]
  pdResults = datFLINOresults
  pdResults = pdResults[is.na(pdResults[,"PARM_EVALUATE_POSITION"]), ]
  pdResults = pdResults[pdResults[,"PARM_EVAL_DAPI_RNDS"] %in% pdTrials,]
  
  
  
  #par(fig=c(0,0.333+xOverlay,0.5+2*yOverlay,1), new=FALSE)
  #c(x1, x2, y1, y2)
  # make labels and margins smaller
  #par(mai=c(0.1,0.1,0.2,0.1))
  par(mar = c(3, 5, 1, 1)) #  bottom, left, top and right  c(5.1, 4.1, 4.1, 2.1).
  
  setDefaultPlotParm()
  DISPLY_GRID = FALSE
  DRAW_BOX_FRAME = FALSE
  acex <<- 1
  tcex <<-  0.9
  lcex <<- 0.9
  t2cex <<- 0.9
  ploty = 1.3
  pys = c(1.25,1.25,1.25, rep(0.3,times=NUM_GROUPS-3))
  ADD_HORIZONTAL_BAR_AT_GIVEN_VALUE  = TRUE
  HORIZONTAL_BAR_VALUE = 0.0738
  pdTable = plotVirtualSlideData(pdResults, pdTrials, pdPlotColumns[1:NUM_GROUPS,], 0 , ploty, "")
  abline(h=0)
  mtext("Error across virtual slides (Coefficient of Variation)", line=3.5, side =2, cex=t2cex)
  mtext("for dapi segmented nuclear object intensities", line=2.5, side =2, cex=t2cex)
  for(gi in c(1:NUM_GROUPS)) {
    text(format(pdTable[5,gi],digits = 3),  x = barCenters[gi],y =pys[gi], cex = tcex, xpd=NA)
  }
  
  if (PLOT_TO != SCREEN) {
    dev.off()
  }
}


###############################################################
###############################################################

if (GENERATE_PANEL_S4) {
  panelName = "Fig_S4"
  
  if (PLOT_TO == TIFF) {
    #jpeg(file.path(pathToFigures,paste0("Dapi",".slide",".jpeg")), quality = 100, width = 720, height = 480)
    myFile = file.path(pathToFigures,paste0(panelName,".tif"))
    tiff(myFile, height = 6, width = 7.5, units = 'in', compression="lzw",  type = "windows", res = 300)
    #par(cex=1.3)
  } else if (PLOT_TO == EPS) {
    myFile = file.path(pathToFigures,paste0(panelName,".eps"))
    setEPS()
    postscript(myFile,  height = 6, width = 7.5)
  }
  value_PARM_NORMALIZE_EVALUATE_PIPELINE = "NORMALIZE_ALL_FOVS_EVALUATE_ALL_FOVS"
  value_PARM_NORM_METHOD = "TMM"
  value_PARM_METRIC = CV_METRIC_NormLogSpace
  NUM_GROUPS = 5
  pdPlotColumns = data.frame(
    PARM_NORM_METHOD=c(rep(value_PARM_NORM_METHOD, times=NUM_GROUPS)),
    PARM_EVAL_SEG_OBJ=c(rep(TRUE, times=NUM_GROUPS)),
    PARM_SEG_OBJ_NAME=c(rep("NucleiSCA", times=NUM_GROUPS)),
    PARM_EVAL_SEG_OBJ_NAME=c(rep("NucleiSCA", times=NUM_GROUPS)),
    PARM_FILTER_MIN_DAPI_RNDtoRND_CORR=c(0,
                                         25,
                                         50,
                                         75,
                                         90),
    Display_X_Line_0=c("No Filter",
                       "<25",
                       "<50",
                       "<75",
                       "<90"),
    PARM_NORMALIZE_EVALUATE_PIPELINE=rep(value_PARM_NORMALIZE_EVALUATE_PIPELINE,times=NUM_GROUPS),
    PARM_METRIC=rep("EvalSegObj_SegObjErrCV_NormLogSpace",times=NUM_GROUPS),
    Display_X_Line_1= rep("", times=NUM_GROUPS), 
    Display_X_Line_2= rep("", times=NUM_GROUPS)) 
  
  
  
  maxVirtualSlides = c(2,3,10,14)
  
  pdTrials = unique(datFLINOresults[,"PARM_EVAL_DAPI_RNDS"])
  pdTrials = pdTrials[(str_count(pdTrials, pattern = ":") + 1) %in% maxVirtualSlides]
  pdResults = datFLINOresults
  pdResults = pdResults[is.na(pdResults[,"PARM_EVALUATE_POSITION"]), ]
  pdResults = pdResults[pdResults[,"PARM_EVAL_DAPI_RNDS"] %in% pdTrials,]
  
  par(fig=c(0.42,1,0.5,1.0), new=FALSE)
  #par(fig=c(0.42,1,0,0.5), new=FALSE)
  par(mar = c(2, 5, 1, 1)) #  bottom, left, top and right  c(5.1, 4.1, 4.1, 2.1).
  
  plot_max_value = 0.2
  plot_display_value = 0.17
  plot_display_text = 0.185
  setDefaultPlotParm()
  pdTable = plotVirtualSlideData(pdResults, pdTrials, pdPlotColumns, 0 , plot_max_value, "")
  mtext("Nuclei Evaluation Object Error (CV)", line=3, side =2, cex=t2cex)
  for(gi in c(1:NUM_GROUPS)) {
    text(format(pdTable[5,gi],digits = 3),  x = barCenters[gi],y =plot_display_value, cex = 0.8, xpd=NA)
  }
  text("CV of Evaluation Objects",  x = barCenters[3],y =plot_display_text, cex = 0.8, xpd=NA)
  abline(v=(barCenters[5]+barCenters[6])/2)
  mtext("Filtering tolerance level", side=1, line = 1, cex = 0.9)

  # correctionMethodTxt = paste0("Correction Method: ", value_PARM_NORM_METHOD, " log-norm")
  # mtext(correctionMethodTxt, side=1, line = 2, cex = 0.9)
  # mtext(value_PARM_NORMALIZE_EVALUATE_PIPELINE, side=1, line = 3, cex = 0.9)
  # correctionMethodTxt = paste0(correctionMethodTxt, " ", value_PARM_NORMALIZE_EVALUATE_PIPELINE, " ", CV_METRIC_TYPE)
  # cat("\n", correctionMethodTxt,"\n")
  # print(pdTable)
  
  
  
  value_PARM_NORMALIZE_EVALUATE_PIPELINE = "NORMALIZE_ALL_FOVS_EVALUATE_ALL_FOVS"
  value_PARM_NORM_METHOD = "TMM"
  value_PARM_METRIC = CV_METRIC_NormLogSpace
  NUM_GROUPS = 10
  pdPlotColumns = data.frame(
    PARM_NORM_METHOD=c(rep(value_PARM_NORM_METHOD, times=NUM_GROUPS)),
    PARM_EVAL_SEG_OBJ=c(rep(TRUE, times=NUM_GROUPS)),
    PARM_SEG_OBJ_NAME=c(rep("NucleiSCA", times=NUM_GROUPS)),
    PARM_EVAL_SEG_OBJ_NAME=c(rep("NucleiSCA", times=NUM_GROUPS)),
    PARM_FILTER_MIN_DAPI_RNDtoRND_CORR=c(0,
                                         25,
                                         50,
                                         75,
                                         90,
                                         0,
                                         25,
                                         50,
                                         75,
                                         90),
    Display_X_Line_0=c("No Filter",
                       "<25",
                       "<50",
                       "<75",
                       "<90",
                       "No Filter",
                       "<25",
                       "<50",
                       "<75",
                       "<90"),
    PARM_NORMALIZE_EVALUATE_PIPELINE=rep(value_PARM_NORMALIZE_EVALUATE_PIPELINE,times=NUM_GROUPS),
    PARM_METRIC=c(rep("Number_NormSegObj",times=NUM_GROUPS/2),rep("Number_EvalSegObj",times=NUM_GROUPS/2)),
    Display_X_Line_1= rep("", times=NUM_GROUPS), 
    Display_X_Line_2= rep("", times=NUM_GROUPS)) 
  
  maxVirtualSlides = c(2,3,10,14)
  pdTrials = unique(datFLINOresults[,"PARM_EVAL_DAPI_RNDS"])
  pdTrials = pdTrials[(str_count(pdTrials, pattern = ":") + 1) %in% maxVirtualSlides]
  pdResults = datFLINOresults
  pdResults = pdResults[is.na(pdResults[,"PARM_EVALUATE_POSITION"]), ]
  pdResults = pdResults[pdResults[,"PARM_EVAL_DAPI_RNDS"] %in% pdTrials,]

  
  par(fig=c(0,1,0,0.5), new=TRUE)
  #par(fig=c(0,1,0.5,1.0), new=TRUE)
  par(mar = c(3, 5, 1, 1)) #  bottom, left, top and right  c(5.1, 4.1, 4.1, 2.1).
  plot_max_value = 400000
  plot_display_value = 340000
  plot_display_text = 370000
  
  setDefaultPlotParm()
  pdTable = plotVirtualSlideData(pdResults, pdTrials, pdPlotColumns, 0 , plot_max_value, "")
  mtext("Number of nuclei objects", line=3.5, side =2, cex=t2cex)
  for(gi in c(1:NUM_GROUPS)) {
    text(format(pdTable[5,gi],digits = 3),  x = barCenters[gi],y =plot_display_value, cex = 0.8, xpd=NA)
  }
  text("Number of Normalization Objects", adj=0, x = barCenters[1],y =plot_display_text, cex = 0.8, xpd=NA)
  text("Number of Evaluation Objects",  x = barCenters[8],y =plot_display_text, cex = 0.8, xpd=NA)
  abline(v=(barCenters[5]+barCenters[6])/2)
  mtext("Filtering tolerance level", side=1, line = 1, cex = 0.9)
  # mtext(value_PARM_NORMALIZE_EVALUATE_PIPELINE, side=1, line = 3, cex = 0.9)
  # correctionMethodTxt = paste0(value_PARM_NORMALIZE_EVALUATE_PIPELINE, " NUMBER OF OBJECTS")
  # cat("\n", correctionMethodTxt,"\n")
  # print(pdTable)
  
  if (PLOT_TO != SCREEN) {
    dev.off()
  }
  
}


###############################################################
###############################################################
#PLOT_TO = SCREEN
if (GENERATE_PANEL_S5) {
  panelName = "Fig_S5"
  
  if (PLOT_TO == TIFF) {
    #jpeg(file.path(pathToFigures,paste0("Dapi",".slide",".jpeg")), quality = 100, width = 720, height = 480)
    myFile = file.path(pathToFigures,paste0(panelName,".tif"))
    tiff(myFile, height = 5, width = 7.5, units = 'in', compression="lzw",  type = "windows", res = 300)
    #par(cex=1.3)
  } else if (PLOT_TO == EPS) {
    myFile = file.path(pathToFigures,paste0(panelName,".eps"))
    setEPS()
    postscript(myFile,  height = 5, width = 7.5)
  }
  
  
  value_PARM_NORMALIZE_EVALUATE_PIPELINE = "NORMALIZE_ALL_FOVS_EVALUATE_ALL_FOVS"
  value_PARM_NORM_METHOD = "MEDIAN" # c("MEDIAN","Q2NORM", "Q3NORM" , "Q625NORM",  "Q875NORM", "Q9375NORM" ,"UPPER_QUARTILE","MRN","TMM")) { 
  filterCorr = 0
  Target_GRID_SIZE = "Grid32"
  value_PARM_METRIC = CV_METRIC_NormLogSpace
  NUM_GROUPS = 10
  pdPlotColumns = data.frame(
    PARM_SEG_OBJ_NAME= c(rep(Target_GRID_SIZE, times=NUM_GROUPS)),
    PARM_EVAL_SEG_OBJ= c(rep(TRUE, times=NUM_GROUPS)),
    PARM_NORM_METHOD=c("Q50NZ",
                       "Q625NZ",
                       "Q75NZ",
                       "Q875NZ",
                       "Q90NZ",
                       "Q95NZ",
                       "Q99NZ",
                       "QMAXNZ",
                       "MRN",
                       "TMM"),
    Display_X_Line_0=c("Q50",
                       "Q62.5",
                       "Q75",
                       "Q87.5",
                       "Q90",
                       "Q95",
                       "Q99",
                       "Q100",
                       "MRN",
                       "TMM"),
    PARM_EVAL_SEG_OBJ_NAME= c(rep("NucleiSCA", times=NUM_GROUPS)),
    PARM_FILTER_MIN_DAPI_RNDtoRND_CORR= c(rep(filterCorr, times=NUM_GROUPS)),
    PARM_NORMALIZE_EVALUATE_PIPELINE=c(rep(value_PARM_NORMALIZE_EVALUATE_PIPELINE, times=NUM_GROUPS)),
    PARM_METRIC=c(rep(value_PARM_METRIC, times=NUM_GROUPS)),
    Display_X_Line_1= c(rep("", times=NUM_GROUPS)),
    Display_X_Line_2= c(rep("", times=NUM_GROUPS)))     
  
  maxVirtualSlides = c(2,3,10,14)
  pdTrials = unique(datFLINOresults[,"PARM_EVAL_DAPI_RNDS"])
  pdTrials = pdTrials[(str_count(pdTrials, pattern = ":") + 1) %in% maxVirtualSlides]
  pdResults = datFLINOresults
  pdResults = pdResults[is.na(pdResults[,"PARM_EVALUATE_POSITION"]), ]
  pdResults = pdResults[pdResults[,"PARM_EVAL_DAPI_RNDS"] %in% pdTrials,]
  
  midpt = 1 - 0.515
  par(fig=c(0,1,0,midpt), new=FALSE)
  par(mar = c(2, 5, 1, 1)) #  bottom, left, top and right  c(5.1, 4.1, 4.1, 2.1).
  
  # par(fig=c(0,0.5,0,1), new=FALSE)
  # par(mar = c(3, 4.5, 1, 0)) #  bottom, left, top and right  c(5.1, 4.1, 4.1, 2.1).
  
  setDefaultPlotParm()
  ADD_HORIZONTAL_BAR_AT_GIVEN_VALUE  = TRUE
  HORIZONTAL_BAR_VALUE = 0.0738
  pdTable = plotVirtualSlideData(pdResults, pdTrials, pdPlotColumns, 0 , 0.3, "")
  mtext("Nuclei Evaluation Error (CV)", line=3, side =2, cex=t2cex)
  text(paste0("Normalization ", get_display_PARM_METRIC(value_PARM_METRIC), " space ",  " - grid size ", gsub("Grid","",Target_GRID_SIZE)), adj=0, x = barCenters[1] / 2,y =0.26, cex = 0.8, xpd=NA)
  for(gi in c(1:NUM_GROUPS)) {
    text(format(pdTable[5,gi],digits = 3),  x = barCenters[gi],y =0.22, cex = 0.8, xpd=NA)
  }
  mtext("Normalization Method", side=1, line = 1, cex = 0.9)
  
  
  
  
  value_PARM_NORM_METHOD = "UPPER_QUARTILE" # c("MEDIAN","UPPER_QUARTILE","MRN","TMM")) { 
  value_PARM_NORM_METHOD_LABEL = "Q75" # c("MEDIAN","UPPER_QUARTILE","MRN","TMM")) { 
  filterCorr = 0
  value_PARM_METRIC = CV_METRIC_NormLogSpace
  NUM_GROUPS = 9
  pdPlotColumns = data.frame(
    PARM_NORM_METHOD= c(rep(value_PARM_NORM_METHOD, times=NUM_GROUPS)),
    PARM_EVAL_SEG_OBJ= c(rep(TRUE, times=NUM_GROUPS)),
    PARM_SEG_OBJ_NAME=c("Grid2560",
                        "Grid1024",
                        "Grid512",
                        "Grid256",
                        "Grid128",
                        "Grid64",
                        "Grid32",
                        "Grid16",
                        "NucleiSCA"),
    Display_X_Line_0=c("2560",
                       "1024",
                       "512",
                       "256",
                       "128",
                       "64",
                       "32",
                       "16",
                       "Nuclei"),
    PARM_EVAL_SEG_OBJ_NAME= c(rep("NucleiSCA", times=NUM_GROUPS)),
    PARM_FILTER_MIN_DAPI_RNDtoRND_CORR= c(rep(filterCorr, times=NUM_GROUPS)),
    PARM_NORMALIZE_EVALUATE_PIPELINE=c(rep(value_PARM_NORMALIZE_EVALUATE_PIPELINE, times=NUM_GROUPS)),
    PARM_METRIC=c(rep(value_PARM_METRIC, times=NUM_GROUPS)),
    Display_X_Line_1= c(rep("", times=NUM_GROUPS)),
    Display_X_Line_2= c(rep("", times=NUM_GROUPS)))     
  
  maxVirtualSlides = c(2,3,10,14)
  pdTrials = unique(datFLINOresults[,"PARM_EVAL_DAPI_RNDS"])
  pdTrials = pdTrials[(str_count(pdTrials, pattern = ":") + 1) %in% maxVirtualSlides]
  pdResults = datFLINOresults
  pdResults = pdResults[is.na(pdResults[,"PARM_EVALUATE_POSITION"]), ]
  pdResults = pdResults[pdResults[,"PARM_EVAL_DAPI_RNDS"] %in% pdTrials,]
  
  par(fig=c(0,1,midpt,1), new=TRUE)
  par(mar = c(3, 5, 1, 1)) #  bottom, left, top and right  c(5.1, 4.1, 4.1, 2.1).
  # par(fig=c(0.5,1,0,1), new=TRUE)
  # par(mar = c(3, 3.5, 1, 1)) #  bottom, left, top and right  c(5.1, 4.1, 4.1, 2.1).
  
  setDefaultPlotParm()
  ADD_HORIZONTAL_BAR_AT_GIVEN_VALUE  = TRUE
  HORIZONTAL_BAR_VALUE = 0.0738
  pdTable = plotVirtualSlideData(pdResults, pdTrials, pdPlotColumns, 0 , 0.2, "")
  mtext("Nuclei Evaluation Error (CV)", line=3, side =2, cex=t2cex)
  text(paste0(value_PARM_NORM_METHOD_LABEL," normalization ", get_display_PARM_METRIC(value_PARM_METRIC), " space"), adj=0, x = barCenters[1] / 2,y =0.18, cex = 0.8, xpd=NA)
  for(gi in c(1:NUM_GROUPS)) {
    text(format(pdTable[5,gi],digits = 3),  x = barCenters[gi],y =0.16, cex = 0.8, xpd=NA)
  }
  mtext("Grid Size (Pixels)", side=1, line = 1, cex = 0.9)
  mtext("Objects", at = barCenters[NUM_GROUPS], side=1, line = 1, cex = 0.9)
  #text("Objects",  x = barCenters[NUM_GROUPS],y =0.16, cex = 0.8, xpd=NA)
  
  # filterTxt = paste0("Filter <", filterCorr, " Corr")
  # if (filterCorr == 0) {
  #   filterTxt ="No Filtering"
  # }
  # if (value_PARM_METRIC == CV_METRIC_NormLogSpace) {
  #   displayMetric = "log-norm"
  # } else if (value_PARM_METRIC == CV_METRIC_NormAbsSpace) {
  #   displayMetric = "abs-norm"
  # }
  # correctionMethodTxt = paste0("Correction Method: ", value_PARM_NORM_METHOD, " ", displayMetric, " ", filterTxt)
  # mtext(correctionMethodTxt, side=1, line = 2, cex = 0.9)
  # mtext(value_PARM_NORMALIZE_EVALUATE_PIPELINE, side=1, line = 3, cex = 0.9)
  # mtext(CV_METRIC_TYPE, side=1, line = 4, cex = 0.9)
  # correctionMethodTxt = paste0(correctionMethodTxt, " ", value_PARM_NORMALIZE_EVALUATE_PIPELINE, " ", CV_METRIC_TYPE)
  # cat("\n", correctionMethodTxt,"\n")
  # print(pdTable)
  
  if (PLOT_TO != SCREEN) {
    dev.off()
  }
}



###############################################################
###############################################################
#PLOT_TO = SCREEN
if (GENERATE_PANEL_S6a) {
  panelName = "Fig_S6a"
  
  if (PLOT_TO == TIFF) {
    #jpeg(file.path(pathToFigures,paste0("Dapi",".slide",".jpeg")), quality = 100, width = 720, height = 480)
    myFile = file.path(pathToFigures,paste0(panelName,".tif"))
    tiff(myFile, height = 5, width = 7.5, units = 'in', compression="lzw",  type = "windows", res = 300)
    #par(cex=1.3)
  } else if (PLOT_TO == EPS) {
    myFile = file.path(pathToFigures,paste0(panelName,".eps"))
    setEPS()
    postscript(myFile,  height = 5, width = 7.5)
  }
  
  value_PARM_NORMALIZE_EVALUATE_PIPELINE = "NORMALIZE_ALL_FOVS_EVALUATE_ALL_FOVS"
  value_PARM_NORM_METHOD = "TMM" # c("MEDIAN","UPPER_QUARTILE","MRN","TMM")) { 
  filterCorr = 75
  value_PARM_METRIC = CV_METRIC_NormLogSpace
  NUM_GROUPS = 9
  pdPlotColumns = data.frame(
    PARM_NORM_METHOD= c(rep(value_PARM_NORM_METHOD, times=NUM_GROUPS)),
    PARM_EVAL_SEG_OBJ= c(rep(TRUE, times=NUM_GROUPS)),
    PARM_SEG_OBJ_NAME=c("Grid2560",
                        "Grid1024",
                        "Grid512",
                        "Grid256",
                        "Grid128",
                        "Grid64",
                        "Grid32",
                        "Grid16",
                        "NucleiSCA"),
    Display_X_Line_0=c("2560",
                       "1024",
                       "512",
                       "256",
                       "128",
                       "64",
                       "32",
                       "16",
                       "NucObj"),
    PARM_EVAL_SEG_OBJ_NAME= c(rep("NucleiSCA", times=NUM_GROUPS)),
    PARM_FILTER_MIN_DAPI_RNDtoRND_CORR= c(rep(filterCorr, times=NUM_GROUPS)),
    PARM_NORMALIZE_EVALUATE_PIPELINE=c(rep(value_PARM_NORMALIZE_EVALUATE_PIPELINE, times=NUM_GROUPS)),
    PARM_METRIC=c(rep(value_PARM_METRIC, times=NUM_GROUPS)),
    Display_X_Line_1= c(rep("", times=NUM_GROUPS)),
    Display_X_Line_2= c(rep("", times=NUM_GROUPS)))     
  
  maxVirtualSlides = c(2,3,10,14)
  pdTrials = unique(datFLINOresults[,"PARM_EVAL_DAPI_RNDS"])
  pdTrials = pdTrials[(str_count(pdTrials, pattern = ":") + 1) %in% maxVirtualSlides]
  pdResults = datFLINOresults
  pdResults = pdResults[is.na(pdResults[,"PARM_EVALUATE_POSITION"]), ]
  pdResults = pdResults[pdResults[,"PARM_EVAL_DAPI_RNDS"] %in% pdTrials,]
  
  midpt = 0.515
  par(fig=c(0,1,midpt,1), new=FALSE)
  par(mar = c(2, 5, 1, 1)) #  bottom, left, top and right  c(5.1, 4.1, 4.1, 2.1).
  
  # par(fig=c(0,0.5,0,1), new=FALSE)
  # par(mar = c(3, 4.5, 1, 0)) #  bottom, left, top and right  c(5.1, 4.1, 4.1, 2.1).
  
  setDefaultPlotParm()
  ADD_HORIZONTAL_BAR_AT_GIVEN_VALUE  = TRUE
  HORIZONTAL_BAR_VALUE = 0.0738
  pdTable = plotVirtualSlideData(pdResults, pdTrials, pdPlotColumns, 0 , 0.2, "")
  mtext("Nuclei Evaluation Error (CV)", line=3, side =2, cex=t2cex)
  text(paste0(value_PARM_NORM_METHOD," normalization ", get_display_PARM_METRIC(value_PARM_METRIC), " space ",  "with dapi-corr filtering <", filterCorr), adj=0, x = barCenters[1] / 2,y =0.18, cex = 0.8, xpd=NA)
  for(gi in c(1:NUM_GROUPS)) {
    text(format(pdTable[5,gi],digits = 3),  x = barCenters[gi],y =0.16, cex = 0.8, xpd=NA)
  }
  mtext("Grid Size (Pixels)", side=1, line = 1, cex = 0.9)
  
  
  
  
  value_PARM_NORM_METHOD = "UPPER_QUARTILE" # c("MEDIAN","UPPER_QUARTILE","MRN","TMM")) { 
  filterCorr = 75
  value_PARM_METRIC = CV_METRIC_NormLogSpace
  NUM_GROUPS = 9
  pdPlotColumns = data.frame(
    PARM_NORM_METHOD= c(rep(value_PARM_NORM_METHOD, times=NUM_GROUPS)),
    PARM_EVAL_SEG_OBJ= c(rep(TRUE, times=NUM_GROUPS)),
    PARM_SEG_OBJ_NAME=c("Grid2560",
                        "Grid1024",
                        "Grid512",
                        "Grid256",
                        "Grid128",
                        "Grid64",
                        "Grid32",
                        "Grid16",
                        "NucleiSCA"),
    Display_X_Line_0=c("2560",
                       "1024",
                       "512",
                       "256",
                       "128",
                       "64",
                       "32",
                       "16",
                       "NucObj"),
    PARM_EVAL_SEG_OBJ_NAME= c(rep("NucleiSCA", times=NUM_GROUPS)),
    PARM_FILTER_MIN_DAPI_RNDtoRND_CORR= c(rep(filterCorr, times=NUM_GROUPS)),
    PARM_NORMALIZE_EVALUATE_PIPELINE=c(rep(value_PARM_NORMALIZE_EVALUATE_PIPELINE, times=NUM_GROUPS)),
    PARM_METRIC=c(rep(value_PARM_METRIC, times=NUM_GROUPS)),
    Display_X_Line_1= c(rep("", times=NUM_GROUPS)),
    Display_X_Line_2= c(rep("", times=NUM_GROUPS)))     
  
  maxVirtualSlides = c(2,3,10,14)
  pdTrials = unique(datFLINOresults[,"PARM_EVAL_DAPI_RNDS"])
  pdTrials = pdTrials[(str_count(pdTrials, pattern = ":") + 1) %in% maxVirtualSlides]
  pdResults = datFLINOresults
  pdResults = pdResults[is.na(pdResults[,"PARM_EVALUATE_POSITION"]), ]
  pdResults = pdResults[pdResults[,"PARM_EVAL_DAPI_RNDS"] %in% pdTrials,]
  
  par(fig=c(0,1,0,midpt), new=TRUE)
  par(mar = c(3, 5, 1, 1)) #  bottom, left, top and right  c(5.1, 4.1, 4.1, 2.1).
  # par(fig=c(0.5,1,0,1), new=TRUE)
  # par(mar = c(3, 3.5, 1, 1)) #  bottom, left, top and right  c(5.1, 4.1, 4.1, 2.1).
  
  setDefaultPlotParm()
  ADD_HORIZONTAL_BAR_AT_GIVEN_VALUE  = TRUE
  HORIZONTAL_BAR_VALUE = 0.0738
  pdTable = plotVirtualSlideData(pdResults, pdTrials, pdPlotColumns, 0 , 0.2, "")
  mtext("Nuclei Evaluation Error (CV)", line=3, side =2, cex=t2cex)
  text(paste0(value_PARM_NORM_METHOD," normalization ", get_display_PARM_METRIC(value_PARM_METRIC), " space ",  "with dapi-corr filtering <", filterCorr), adj=0, x = barCenters[1] / 2,y =0.18, cex = 0.8, xpd=NA)
  for(gi in c(1:NUM_GROUPS)) {
    text(format(pdTable[5,gi],digits = 3),  x = barCenters[gi],y =0.16, cex = 0.8, xpd=NA)
  }
  mtext("Grid Size (Pixels)", side=1, line = 1, cex = 0.9)
  
  # filterTxt = paste0("Filter <", filterCorr, " Corr")
  # if (filterCorr == 0) {
  #   filterTxt ="No Filtering"
  # }
  # if (value_PARM_METRIC == CV_METRIC_NormLogSpace) {
  #   displayMetric = "log-norm"
  # } else if (value_PARM_METRIC == CV_METRIC_NormAbsSpace) {
  #   displayMetric = "abs-norm"
  # }
  # correctionMethodTxt = paste0("Correction Method: ", value_PARM_NORM_METHOD, " ", displayMetric, " ", filterTxt)
  # mtext(correctionMethodTxt, side=1, line = 2, cex = 0.9)
  # mtext(value_PARM_NORMALIZE_EVALUATE_PIPELINE, side=1, line = 3, cex = 0.9)
  # mtext(CV_METRIC_TYPE, side=1, line = 4, cex = 0.9)
  # correctionMethodTxt = paste0(correctionMethodTxt, " ", value_PARM_NORMALIZE_EVALUATE_PIPELINE, " ", CV_METRIC_TYPE)
  # cat("\n", correctionMethodTxt,"\n")
  # print(pdTable)
  
  if (PLOT_TO != SCREEN) {
    dev.off()
  }
}



###############################################################
###############################################################

#PLOT_TO = SCREEN
if (GENERATE_PANEL_S6b) {
  panelName = "Fig_S6b"
  
  if (PLOT_TO == TIFF) {
    #jpeg(file.path(pathToFigures,paste0("Dapi",".slide",".jpeg")), quality = 100, width = 720, height = 480)
    myFile = file.path(pathToFigures,paste0(panelName,".tif"))
    tiff(myFile, height = 5, width = 7.5, units = 'in', compression="lzw",  type = "windows", res = 300)
    #par(cex=1.3)
  } else if (PLOT_TO == EPS) {
    myFile = file.path(pathToFigures,paste0(panelName,".eps"))
    setEPS()
    postscript(myFile,  height = 5, width = 7.5)
  }
  
  value_PARM_NORMALIZE_EVALUATE_PIPELINE = "NORMALIZE_ALL_FOVS_EVALUATE_ALL_FOVS"
  value_PARM_NORM_METHOD = "MRN" # c("MEDIAN","UPPER_QUARTILE","MRN","TMM")) { 
  filterCorr = 0
  value_PARM_METRIC = CV_METRIC_NormLogSpace
  NUM_GROUPS = 9
  pdPlotColumns = data.frame(
    PARM_NORM_METHOD= c(rep(value_PARM_NORM_METHOD, times=NUM_GROUPS)),
    PARM_EVAL_SEG_OBJ= c(rep(TRUE, times=NUM_GROUPS)),
    PARM_SEG_OBJ_NAME=c("Grid2560",
                        "Grid1024",
                        "Grid512",
                        "Grid256",
                        "Grid128",
                        "Grid64",
                        "Grid32",
                        "Grid16",
                        "NucleiSCA"),
    Display_X_Line_0=c("2560",
                       "1024",
                       "512",
                       "256",
                       "128",
                       "64",
                       "32",
                       "16",
                       "NucObj"),
    PARM_EVAL_SEG_OBJ_NAME= c(rep("NucleiSCA", times=NUM_GROUPS)),
    PARM_FILTER_MIN_DAPI_RNDtoRND_CORR= c(rep(filterCorr, times=NUM_GROUPS)),
    PARM_NORMALIZE_EVALUATE_PIPELINE=c(rep(value_PARM_NORMALIZE_EVALUATE_PIPELINE, times=NUM_GROUPS)),
    PARM_METRIC=c(rep(value_PARM_METRIC, times=NUM_GROUPS)),
    Display_X_Line_1= c(rep("", times=NUM_GROUPS)),
    Display_X_Line_2= c(rep("", times=NUM_GROUPS)))     
  
  maxVirtualSlides = c(2,3,10,14)
  pdTrials = unique(datFLINOresults[,"PARM_EVAL_DAPI_RNDS"])
  pdTrials = pdTrials[(str_count(pdTrials, pattern = ":") + 1) %in% maxVirtualSlides]
  pdResults = datFLINOresults
  pdResults = pdResults[is.na(pdResults[,"PARM_EVALUATE_POSITION"]), ]
  pdResults = pdResults[pdResults[,"PARM_EVAL_DAPI_RNDS"] %in% pdTrials,]
  
  midpt = 0.515
  par(fig=c(0,1,midpt,1), new=FALSE)
  par(mar = c(2, 5, 1, 1)) #  bottom, left, top and right  c(5.1, 4.1, 4.1, 2.1).
  
  # par(fig=c(0,0.5,0,1), new=FALSE)
  # par(mar = c(3, 4.5, 1, 0)) #  bottom, left, top and right  c(5.1, 4.1, 4.1, 2.1).
  
  setDefaultPlotParm()
  ADD_HORIZONTAL_BAR_AT_GIVEN_VALUE  = TRUE
  HORIZONTAL_BAR_VALUE = 0.0738
  pdTable = plotVirtualSlideData(pdResults, pdTrials, pdPlotColumns, 0 , 0.2, "")
  mtext("Nuclei Evaluation Error (CV)", line=3, side =2, cex=t2cex)
  text(paste0(value_PARM_NORM_METHOD," normalization ", get_display_PARM_METRIC(value_PARM_METRIC), " space"), adj=0, x = barCenters[1] / 2,y =0.18, cex = 0.8, xpd=NA)
  for(gi in c(1:NUM_GROUPS)) {
    text(format(pdTable[5,gi],digits = 3),  x = barCenters[gi],y =0.16, cex = 0.8, xpd=NA)
  }
  mtext("Grid Size (Pixels)", side=1, line = 1, cex = 0.9)
  
  
  
  
  value_PARM_NORM_METHOD = "MEDIAN" # c("MEDIAN","UPPER_QUARTILE","MRN","TMM")) { 
  filterCorr = 0
  value_PARM_METRIC = CV_METRIC_NormLogSpace
  NUM_GROUPS = 9
  pdPlotColumns = data.frame(
    PARM_NORM_METHOD= c(rep(value_PARM_NORM_METHOD, times=NUM_GROUPS)),
    PARM_EVAL_SEG_OBJ= c(rep(TRUE, times=NUM_GROUPS)),
    PARM_SEG_OBJ_NAME=c("Grid2560",
                        "Grid1024",
                        "Grid512",
                        "Grid256",
                        "Grid128",
                        "Grid64",
                        "Grid32",
                        "Grid16",
                        "NucleiSCA"),
    Display_X_Line_0=c("2560",
                       "1024",
                       "512",
                       "256",
                       "128",
                       "64",
                       "32",
                       "16",
                       "NucObj"),
    PARM_EVAL_SEG_OBJ_NAME= c(rep("NucleiSCA", times=NUM_GROUPS)),
    PARM_FILTER_MIN_DAPI_RNDtoRND_CORR= c(rep(filterCorr, times=NUM_GROUPS)),
    PARM_NORMALIZE_EVALUATE_PIPELINE=c(rep(value_PARM_NORMALIZE_EVALUATE_PIPELINE, times=NUM_GROUPS)),
    PARM_METRIC=c(rep(value_PARM_METRIC, times=NUM_GROUPS)),
    Display_X_Line_1= c(rep("", times=NUM_GROUPS)),
    Display_X_Line_2= c(rep("", times=NUM_GROUPS)))     
  
  maxVirtualSlides = c(2,3,10,14)
  pdTrials = unique(datFLINOresults[,"PARM_EVAL_DAPI_RNDS"])
  pdTrials = pdTrials[(str_count(pdTrials, pattern = ":") + 1) %in% maxVirtualSlides]
  pdResults = datFLINOresults
  pdResults = pdResults[is.na(pdResults[,"PARM_EVALUATE_POSITION"]), ]
  pdResults = pdResults[pdResults[,"PARM_EVAL_DAPI_RNDS"] %in% pdTrials,]
  
  par(fig=c(0,1,0,midpt), new=TRUE)
  par(mar = c(3, 5, 1, 1)) #  bottom, left, top and right  c(5.1, 4.1, 4.1, 2.1).
  # par(fig=c(0.5,1,0,1), new=TRUE)
  # par(mar = c(3, 3.5, 1, 1)) #  bottom, left, top and right  c(5.1, 4.1, 4.1, 2.1).
  
  setDefaultPlotParm()
  ADD_HORIZONTAL_BAR_AT_GIVEN_VALUE  = TRUE
  HORIZONTAL_BAR_VALUE = 0.0738
  pdTable = plotVirtualSlideData(pdResults, pdTrials, pdPlotColumns, 0 , 0.4, "")
  mtext("Nuclei Evaluation Error (CV)", line=3, side =2, cex=t2cex)
  text(paste0(value_PARM_NORM_METHOD," normalization ", get_display_PARM_METRIC(value_PARM_METRIC), " space"), adj=0, x = barCenters[1] / 2,y =0.36, cex = 0.8, xpd=NA)
  for(gi in c(1:NUM_GROUPS)) {
    text(format(pdTable[5,gi],digits = 3),  x = barCenters[gi],y =0.32, cex = 0.8, xpd=NA)
  }
  mtext("Grid Size (Pixels)", side=1, line = 1, cex = 0.9)
  
  # filterTxt = paste0("Filter <", filterCorr, " Corr")
  # if (filterCorr == 0) {
  #   filterTxt ="No Filtering"
  # }
  # if (value_PARM_METRIC == CV_METRIC_NormLogSpace) {
  #   displayMetric = "log-norm"
  # } else if (value_PARM_METRIC == CV_METRIC_NormAbsSpace) {
  #   displayMetric = "abs-norm"
  # }
  # correctionMethodTxt = paste0("Correction Method: ", value_PARM_NORM_METHOD, " ", displayMetric, " ", filterTxt)
  # mtext(correctionMethodTxt, side=1, line = 2, cex = 0.9)
  # mtext(value_PARM_NORMALIZE_EVALUATE_PIPELINE, side=1, line = 3, cex = 0.9)
  # mtext(CV_METRIC_TYPE, side=1, line = 4, cex = 0.9)
  # correctionMethodTxt = paste0(correctionMethodTxt, " ", value_PARM_NORMALIZE_EVALUATE_PIPELINE, " ", CV_METRIC_TYPE)
  # cat("\n", correctionMethodTxt,"\n")
  # print(pdTable)
  
  if (PLOT_TO != SCREEN) {
    dev.off()
  }
}



###############################################################
###############################################################

#PLOT_TO = SCREEN
if (GENERATE_PANEL_S6c) {
  panelName = "Fig_S6c"
  
  if (PLOT_TO == TIFF) {
    #jpeg(file.path(pathToFigures,paste0("Dapi",".slide",".jpeg")), quality = 100, width = 720, height = 480)
    myFile = file.path(pathToFigures,paste0(panelName,".tif"))
    tiff(myFile, height = 5, width = 7.5, units = 'in', compression="lzw",  type = "windows", res = 300)
    #par(cex=1.3)
  } else if (PLOT_TO == EPS) {
    myFile = file.path(pathToFigures,paste0(panelName,".eps"))
    setEPS()
    postscript(myFile,  height = 5, width = 7.5)
  }
  
  value_PARM_NORMALIZE_EVALUATE_PIPELINE = "NORMALIZE_ALL_FOVS_EVALUATE_ALL_FOVS"
  value_PARM_NORM_METHOD = "MEDIAN" # c("MEDIAN","Q2NORM", "Q3NORM" , "Q625NORM",  "Q875NORM", "Q9375NORM" ,"UPPER_QUARTILE","MRN","TMM")) { 
  filterCorr = 0
  value_PARM_METRIC = CV_METRIC_NormLogSpace
  NUM_GROUPS = 9
  pdPlotColumns = data.frame(
    PARM_NORM_METHOD= c(rep(value_PARM_NORM_METHOD, times=NUM_GROUPS)),
    PARM_EVAL_SEG_OBJ= c(rep(TRUE, times=NUM_GROUPS)),
    PARM_SEG_OBJ_NAME=c("Grid2560",
                        "Grid1024",
                        "Grid512",
                        "Grid256",
                        "Grid128",
                        "Grid64",
                        "Grid32",
                        "Grid16",
                        "NucleiSCA"),
    Display_X_Line_0=c("2560",
                       "1024",
                       "512",
                       "256",
                       "128",
                       "64",
                       "32",
                       "16",
                       "NucObj"),
    PARM_EVAL_SEG_OBJ_NAME= c(rep("NucleiSCA", times=NUM_GROUPS)),
    PARM_FILTER_MIN_DAPI_RNDtoRND_CORR= c(rep(filterCorr, times=NUM_GROUPS)),
    PARM_NORMALIZE_EVALUATE_PIPELINE=c(rep(value_PARM_NORMALIZE_EVALUATE_PIPELINE, times=NUM_GROUPS)),
    PARM_METRIC=c(rep(value_PARM_METRIC, times=NUM_GROUPS)),
    Display_X_Line_1= c(rep("", times=NUM_GROUPS)),
    Display_X_Line_2= c(rep("", times=NUM_GROUPS)))     
  
  maxVirtualSlides = c(2,3,10,14)
  pdTrials = unique(datFLINOresults[,"PARM_EVAL_DAPI_RNDS"])
  pdTrials = pdTrials[(str_count(pdTrials, pattern = ":") + 1) %in% maxVirtualSlides]
  pdResults = datFLINOresults
  pdResults = pdResults[is.na(pdResults[,"PARM_EVALUATE_POSITION"]), ]
  pdResults = pdResults[pdResults[,"PARM_EVAL_DAPI_RNDS"] %in% pdTrials,]
  
  midpt = 0.515
  par(fig=c(0,1,midpt,1), new=FALSE)
  par(mar = c(2, 5, 1, 1)) #  bottom, left, top and right  c(5.1, 4.1, 4.1, 2.1).
  
  # par(fig=c(0,0.5,0,1), new=FALSE)
  # par(mar = c(3, 4.5, 1, 0)) #  bottom, left, top and right  c(5.1, 4.1, 4.1, 2.1).
  
  setDefaultPlotParm()
  ADD_HORIZONTAL_BAR_AT_GIVEN_VALUE  = TRUE
  HORIZONTAL_BAR_VALUE = 0.0738
  pdTable = plotVirtualSlideData(pdResults, pdTrials, pdPlotColumns, 0 , 0.4, "")
  mtext("Nuclei Evaluation Error (CV)", line=3, side =2, cex=t2cex)
  text(paste0(value_PARM_NORM_METHOD," normalization ", get_display_PARM_METRIC(value_PARM_METRIC), " space"), adj=0, x = barCenters[1] / 2,y =0.36, cex = 0.8, xpd=NA)
  for(gi in c(1:NUM_GROUPS)) {
    text(format(pdTable[5,gi],digits = 3),  x = barCenters[gi],y =0.32, cex = 0.8, xpd=NA)
  }
  mtext("Grid Size (Pixels)", side=1, line = 1, cex = 0.9)
  
  
  
  
  value_PARM_NORM_METHOD = "Q3NORM" # c("MEDIAN","Q2NORM", "Q3NORM" , "Q625NORM",  "Q875NORM", "Q9375NORM","UPPER_QUARTILE","MRN","TMM")) { 
  filterCorr = 0
  value_PARM_METRIC = CV_METRIC_NormLogSpace
  NUM_GROUPS = 9
  pdPlotColumns = data.frame(
    PARM_NORM_METHOD= c(rep(value_PARM_NORM_METHOD, times=NUM_GROUPS)),
    PARM_EVAL_SEG_OBJ= c(rep(TRUE, times=NUM_GROUPS)),
    PARM_SEG_OBJ_NAME=c("Grid2560",
                        "Grid1024",
                        "Grid512",
                        "Grid256",
                        "Grid128",
                        "Grid64",
                        "Grid32",
                        "Grid16",
                        "NucleiSCA"),
    Display_X_Line_0=c("2560",
                       "1024",
                       "512",
                       "256",
                       "128",
                       "64",
                       "32",
                       "16",
                       "NucObj"),
    PARM_EVAL_SEG_OBJ_NAME= c(rep("NucleiSCA", times=NUM_GROUPS)),
    PARM_FILTER_MIN_DAPI_RNDtoRND_CORR= c(rep(filterCorr, times=NUM_GROUPS)),
    PARM_NORMALIZE_EVALUATE_PIPELINE=c(rep(value_PARM_NORMALIZE_EVALUATE_PIPELINE, times=NUM_GROUPS)),
    PARM_METRIC=c(rep(value_PARM_METRIC, times=NUM_GROUPS)),
    Display_X_Line_1= c(rep("", times=NUM_GROUPS)),
    Display_X_Line_2= c(rep("", times=NUM_GROUPS)))     
  
  maxVirtualSlides = c(2,3,10,14)
  pdTrials = unique(datFLINOresults[,"PARM_EVAL_DAPI_RNDS"])
  pdTrials = pdTrials[(str_count(pdTrials, pattern = ":") + 1) %in% maxVirtualSlides]
  pdResults = datFLINOresults
  pdResults = pdResults[is.na(pdResults[,"PARM_EVALUATE_POSITION"]), ]
  pdResults = pdResults[pdResults[,"PARM_EVAL_DAPI_RNDS"] %in% pdTrials,]
  
  par(fig=c(0,1,0,midpt), new=TRUE)
  par(mar = c(3, 5, 1, 1)) #  bottom, left, top and right  c(5.1, 4.1, 4.1, 2.1).
  # par(fig=c(0.5,1,0,1), new=TRUE)
  # par(mar = c(3, 3.5, 1, 1)) #  bottom, left, top and right  c(5.1, 4.1, 4.1, 2.1).
  
  setDefaultPlotParm()
  ADD_HORIZONTAL_BAR_AT_GIVEN_VALUE  = TRUE
  HORIZONTAL_BAR_VALUE = 0.0738
  pdTable = plotVirtualSlideData(pdResults, pdTrials, pdPlotColumns, 0 , 0.4, "")
  mtext("Nuclei Evaluation Error (CV)", line=3, side =2, cex=t2cex)
  text(paste0(value_PARM_NORM_METHOD," normalization ", get_display_PARM_METRIC(value_PARM_METRIC), " space"), adj=0, x = barCenters[1] / 2,y =0.36, cex = 0.8, xpd=NA)
  for(gi in c(1:NUM_GROUPS)) {
    text(format(pdTable[5,gi],digits = 3),  x = barCenters[gi],y =0.32, cex = 0.8, xpd=NA)
  }
  mtext("Grid Size (Pixels)", side=1, line = 1, cex = 0.9)
  
  # filterTxt = paste0("Filter <", filterCorr, " Corr")
  # if (filterCorr == 0) {
  #   filterTxt ="No Filtering"
  # }
  # if (value_PARM_METRIC == CV_METRIC_NormLogSpace) {
  #   displayMetric = "log-norm"
  # } else if (value_PARM_METRIC == CV_METRIC_NormAbsSpace) {
  #   displayMetric = "abs-norm"
  # }
  # correctionMethodTxt = paste0("Correction Method: ", value_PARM_NORM_METHOD, " ", displayMetric, " ", filterTxt)
  # mtext(correctionMethodTxt, side=1, line = 2, cex = 0.9)
  # mtext(value_PARM_NORMALIZE_EVALUATE_PIPELINE, side=1, line = 3, cex = 0.9)
  # mtext(CV_METRIC_TYPE, side=1, line = 4, cex = 0.9)
  # correctionMethodTxt = paste0(correctionMethodTxt, " ", value_PARM_NORMALIZE_EVALUATE_PIPELINE, " ", CV_METRIC_TYPE)
  # cat("\n", correctionMethodTxt,"\n")
  # print(pdTable)
  
  if (PLOT_TO != SCREEN) {
    dev.off()
  }
}



###############################################################
###############################################################

#PLOT_TO = SCREEN
if (GENERATE_PANEL_S6d) {
  panelName = "Fig_S6d"
  
  if (PLOT_TO == TIFF) {
    #jpeg(file.path(pathToFigures,paste0("Dapi",".slide",".jpeg")), quality = 100, width = 720, height = 480)
    myFile = file.path(pathToFigures,paste0(panelName,".tif"))
    tiff(myFile, height = 5, width = 7.5, units = 'in', compression="lzw",  type = "windows", res = 300)
    #par(cex=1.3)
  } else if (PLOT_TO == EPS) {
    myFile = file.path(pathToFigures,paste0(panelName,".eps"))
    setEPS()
    postscript(myFile,  height = 5, width = 7.5)
  }
  
  value_PARM_NORMALIZE_EVALUATE_PIPELINE = "NORMALIZE_ALL_FOVS_EVALUATE_ALL_FOVS"
  value_PARM_NORM_METHOD = "MEDIAN" # c("MEDIAN","Q2NORM", "Q3NORM" , "Q625NORM",  "Q875NORM", "Q9375NORM" ,"UPPER_QUARTILE","MRN","TMM")) { 
  filterCorr = 0
  Target_GRID_SIZE = "Grid32"
  value_PARM_METRIC = CV_METRIC_NormLogSpace
  NUM_GROUPS = 11
  pdPlotColumns = data.frame(
    PARM_SEG_OBJ_NAME= c(rep(Target_GRID_SIZE, times=NUM_GROUPS)),
    PARM_EVAL_SEG_OBJ= c(rep(TRUE, times=NUM_GROUPS)),
    PARM_NORM_METHOD=c("MEDIAN",
                        "Q625NORM",
                        "Q3NORM",
                        "Q875NORM",
                       "Q90NORM",
                       "Q95NORM",
                       "Q99NORM",
                       "MAX",
                        "UPPER_QUARTILE",
                        "MRN",
                        "TMM"),
    Display_X_Line_0=c("Median",
                       "N62.5%",
                       "N75%",
                       "N87.5%",
                       "N90%",
                       "N95%",
                       "N99%",
                       "N100%",
                       "uqua",
                       "MRN",
                       "TMM"),
    PARM_EVAL_SEG_OBJ_NAME= c(rep("NucleiSCA", times=NUM_GROUPS)),
    PARM_FILTER_MIN_DAPI_RNDtoRND_CORR= c(rep(filterCorr, times=NUM_GROUPS)),
    PARM_NORMALIZE_EVALUATE_PIPELINE=c(rep(value_PARM_NORMALIZE_EVALUATE_PIPELINE, times=NUM_GROUPS)),
    PARM_METRIC=c(rep(value_PARM_METRIC, times=NUM_GROUPS)),
    Display_X_Line_1= c(rep("", times=NUM_GROUPS)),
    Display_X_Line_2= c(rep("", times=NUM_GROUPS)))     
  
  maxVirtualSlides = c(2,3,10,14)
  pdTrials = unique(datFLINOresults[,"PARM_EVAL_DAPI_RNDS"])
  pdTrials = pdTrials[(str_count(pdTrials, pattern = ":") + 1) %in% maxVirtualSlides]
  pdResults = datFLINOresults
  pdResults = pdResults[is.na(pdResults[,"PARM_EVALUATE_POSITION"]), ]
  pdResults = pdResults[pdResults[,"PARM_EVAL_DAPI_RNDS"] %in% pdTrials,]
  
  midpt = 0.515
  par(fig=c(0,1,midpt,1), new=FALSE)
  par(mar = c(2, 5, 1, 1)) #  bottom, left, top and right  c(5.1, 4.1, 4.1, 2.1).
  
  # par(fig=c(0,0.5,0,1), new=FALSE)
  # par(mar = c(3, 4.5, 1, 0)) #  bottom, left, top and right  c(5.1, 4.1, 4.1, 2.1).
  
  setDefaultPlotParm()
  ADD_HORIZONTAL_BAR_AT_GIVEN_VALUE  = TRUE
  HORIZONTAL_BAR_VALUE = 0.0738
  pdTable = plotVirtualSlideData(pdResults, pdTrials, pdPlotColumns, 0 , 0.4, "")
  mtext("Nuclei Evaluation Error (CV)", line=3, side =2, cex=t2cex)
  text(paste0("Normalization ", get_display_PARM_METRIC(value_PARM_METRIC), " space ",  " - grid size ", gsub("Grid","",Target_GRID_SIZE)), adj=0, x = barCenters[1] / 2,y =0.36, cex = 0.8, xpd=NA)
  for(gi in c(1:NUM_GROUPS)) {
    text(format(pdTable[5,gi],digits = 3),  x = barCenters[gi],y =0.32, cex = 0.8, xpd=NA)
  }
  mtext("Normalization Method", side=1, line = 1, cex = 0.9)
  
  
  
  
  value_PARM_NORMALIZE_EVALUATE_PIPELINE = "NORMALIZE_ALL_FOVS_EVALUATE_ALL_FOVS"
  value_PARM_NORM_METHOD = "MEDIAN" # c("MEDIAN","Q2NORM", "Q3NORM" , "Q625NORM",  "Q875NORM", "Q9375NORM" ,"UPPER_QUARTILE","MRN","TMM")) { 
  filterCorr = 0
  Target_GRID_SIZE = "Grid32"
  value_PARM_METRIC = CV_METRIC_NormAbsSpace
  NUM_GROUPS = 11
  pdPlotColumns = data.frame(
    PARM_SEG_OBJ_NAME= c(rep(Target_GRID_SIZE, times=NUM_GROUPS)),
    PARM_EVAL_SEG_OBJ= c(rep(TRUE, times=NUM_GROUPS)),
    PARM_NORM_METHOD=c("MEDIAN",
                       "Q625NORM",
                       "Q3NORM",
                       "Q875NORM",
                       "Q90NORM",
                       "Q95NORM",
                       "Q99NORM",
                       "MAX",
                       "UPPER_QUARTILE",
                       "MRN",
                       "TMM"),
    Display_X_Line_0=c("Median",
                       "N62.5%",
                       "N75%",
                       "N87.5%",
                       "N90%",
                       "N95%",
                       "N99%",
                       "N100%",
                       "uqua",
                       "MRN",
                       "TMM"),
    PARM_EVAL_SEG_OBJ_NAME= c(rep("NucleiSCA", times=NUM_GROUPS)),
    PARM_FILTER_MIN_DAPI_RNDtoRND_CORR= c(rep(filterCorr, times=NUM_GROUPS)),
    PARM_NORMALIZE_EVALUATE_PIPELINE=c(rep(value_PARM_NORMALIZE_EVALUATE_PIPELINE, times=NUM_GROUPS)),
    PARM_METRIC=c(rep(value_PARM_METRIC, times=NUM_GROUPS)),
    Display_X_Line_1= c(rep("", times=NUM_GROUPS)),
    Display_X_Line_2= c(rep("", times=NUM_GROUPS)))     
  
  maxVirtualSlides = c(2,3,10,14)
  pdTrials = unique(datFLINOresults[,"PARM_EVAL_DAPI_RNDS"])
  pdTrials = pdTrials[(str_count(pdTrials, pattern = ":") + 1) %in% maxVirtualSlides]
  pdResults = datFLINOresults
  pdResults = pdResults[is.na(pdResults[,"PARM_EVALUATE_POSITION"]), ]
  pdResults = pdResults[pdResults[,"PARM_EVAL_DAPI_RNDS"] %in% pdTrials,]
  
  par(fig=c(0,1,0,midpt), new=TRUE)
  par(mar = c(3, 5, 1, 1)) #  bottom, left, top and right  c(5.1, 4.1, 4.1, 2.1).
  # par(fig=c(0.5,1,0,1), new=TRUE)
  # par(mar = c(3, 3.5, 1, 1)) #  bottom, left, top and right  c(5.1, 4.1, 4.1, 2.1).
  
  setDefaultPlotParm()
  ADD_HORIZONTAL_BAR_AT_GIVEN_VALUE  = TRUE
  HORIZONTAL_BAR_VALUE = 0.0738
  pdTable = plotVirtualSlideData(pdResults, pdTrials, pdPlotColumns, 0 , 0.8, "")
  mtext("Nuclei Evaluation Error (CV)", line=3, side =2, cex=t2cex)
  text(paste0("Normalization ", get_display_PARM_METRIC(value_PARM_METRIC), " space ",  " - grid size ", gsub("Grid","",Target_GRID_SIZE)), adj=0, x = barCenters[1] / 2,y =0.72, cex = 0.8, xpd=NA)
  for(gi in c(1:NUM_GROUPS)) {
    text(format(pdTable[5,gi],digits = 3),  x = barCenters[gi],y =0.64, cex = 0.8, xpd=NA)
  }
  mtext("Normalization Method", side=1, line = 1, cex = 0.9)
  
  
  
  if (PLOT_TO != SCREEN) {
    dev.off()
  }
}



###############################################################
###############################################################

#PLOT_TO = SCREEN
if (GENERATE_PANEL_S6e) {
  panelName = "Fig_S6e"
  
  if (PLOT_TO == TIFF) {
    #jpeg(file.path(pathToFigures,paste0("Dapi",".slide",".jpeg")), quality = 100, width = 720, height = 480)
    myFile = file.path(pathToFigures,paste0(panelName,".tif"))
    tiff(myFile, height = 5, width = 7.5, units = 'in', compression="lzw",  type = "windows", res = 300)
    #par(cex=1.3)
  } else if (PLOT_TO == EPS) {
    myFile = file.path(pathToFigures,paste0(panelName,".eps"))
    setEPS()
    postscript(myFile,  height = 5, width = 7.5)
  }
  
  value_PARM_NORMALIZE_EVALUATE_PIPELINE = "NORMALIZE_ALL_FOVS_EVALUATE_ALL_FOVS"
  value_PARM_NORM_METHOD = "MEDIAN" # c("MEDIAN","Q2NORM", "Q3NORM" , "Q625NORM",  "Q875NORM", "Q9375NORM" ,"UPPER_QUARTILE","MRN","TMM")) { 
  filterCorr = 0
  Target_GRID_SIZE = "Grid32"
  value_PARM_METRIC = CV_METRIC_NormLogSpace
  NUM_GROUPS = 11
  pdPlotColumns = data.frame(
    PARM_SEG_OBJ_NAME= c(rep(Target_GRID_SIZE, times=NUM_GROUPS)),
    PARM_EVAL_SEG_OBJ= c(rep(TRUE, times=NUM_GROUPS)),
    PARM_NORM_METHOD=c("Q50NZ",
                       "Q625NZ",
                       "Q75NZ",
                       "Q875NZ",
                       "Q90NZ",
                       "Q95NZ",
                       "Q99NZ",
                       "QMAXNZ",
                       "UPPER_QUARTILE",
                       "MRN",
                       "TMM"),
    Display_X_Line_0=c("50%",
                       "62.5%",
                       "75%",
                       "87.5%",
                       "90%",
                       "95%",
                       "99%",
                       "100%",
                       "uqua",
                       "MRN",
                       "TMM"),
    PARM_EVAL_SEG_OBJ_NAME= c(rep("NucleiSCA", times=NUM_GROUPS)),
    PARM_FILTER_MIN_DAPI_RNDtoRND_CORR= c(rep(filterCorr, times=NUM_GROUPS)),
    PARM_NORMALIZE_EVALUATE_PIPELINE=c(rep(value_PARM_NORMALIZE_EVALUATE_PIPELINE, times=NUM_GROUPS)),
    PARM_METRIC=c(rep(value_PARM_METRIC, times=NUM_GROUPS)),
    Display_X_Line_1= c(rep("", times=NUM_GROUPS)),
    Display_X_Line_2= c(rep("", times=NUM_GROUPS)))     
  
  maxVirtualSlides = c(2,3,10,14)
  pdTrials = unique(datFLINOresults[,"PARM_EVAL_DAPI_RNDS"])
  pdTrials = pdTrials[(str_count(pdTrials, pattern = ":") + 1) %in% maxVirtualSlides]
  pdResults = datFLINOresults
  pdResults = pdResults[is.na(pdResults[,"PARM_EVALUATE_POSITION"]), ]
  pdResults = pdResults[pdResults[,"PARM_EVAL_DAPI_RNDS"] %in% pdTrials,]
  
  midpt = 0.515
  par(fig=c(0,1,midpt,1), new=FALSE)
  par(mar = c(2, 5, 1, 1)) #  bottom, left, top and right  c(5.1, 4.1, 4.1, 2.1).
  
  # par(fig=c(0,0.5,0,1), new=FALSE)
  # par(mar = c(3, 4.5, 1, 0)) #  bottom, left, top and right  c(5.1, 4.1, 4.1, 2.1).
  
  setDefaultPlotParm()
  ADD_HORIZONTAL_BAR_AT_GIVEN_VALUE  = TRUE
  HORIZONTAL_BAR_VALUE = 0.0738
  pdTable = plotVirtualSlideData(pdResults, pdTrials, pdPlotColumns, 0 , 0.3, "")
  mtext("Nuclei Evaluation Error (CV)", line=3, side =2, cex=t2cex)
  text(paste0("Normalization ", get_display_PARM_METRIC(value_PARM_METRIC), " space ",  " - grid size ", gsub("Grid","",Target_GRID_SIZE)), adj=0, x = barCenters[1] / 2,y =0.26, cex = 0.8, xpd=NA)
  for(gi in c(1:NUM_GROUPS)) {
    text(format(pdTable[5,gi],digits = 3),  x = barCenters[gi],y =0.22, cex = 0.8, xpd=NA)
  }
  mtext("Normalization Method", side=1, line = 1, cex = 0.9)
  
  
  
  
  value_PARM_NORMALIZE_EVALUATE_PIPELINE = "NORMALIZE_ALL_FOVS_EVALUATE_ALL_FOVS"
  value_PARM_NORM_METHOD = "MEDIAN" # c("MEDIAN","Q2NORM", "Q3NORM" , "Q625NORM",  "Q875NORM", "Q9375NORM" ,"UPPER_QUARTILE","MRN","TMM")) { 
  filterCorr = 0
  Target_GRID_SIZE = "Grid32"
  value_PARM_METRIC = CV_METRIC_NormAbsSpace
  NUM_GROUPS = 11
  pdPlotColumns = data.frame(
    PARM_SEG_OBJ_NAME= c(rep(Target_GRID_SIZE, times=NUM_GROUPS)),
    PARM_EVAL_SEG_OBJ= c(rep(TRUE, times=NUM_GROUPS)),
    PARM_NORM_METHOD=c("Q50NZ",
                       "Q625NZ",
                       "Q75NZ",
                       "Q875NZ",
                       "Q90NZ",
                       "Q95NZ",
                       "Q99NZ",
                       "QMAXNZ",
                       "UPPER_QUARTILE",
                       "MRN",
                       "TMM"),
    Display_X_Line_0=c("50%",
                       "62.5%",
                       "75%",
                       "87.5%",
                       "90%",
                       "95%",
                       "99%",
                       "100%",
                       "uqua",
                       "MRN",
                       "TMM"),
    PARM_EVAL_SEG_OBJ_NAME= c(rep("NucleiSCA", times=NUM_GROUPS)),
    PARM_FILTER_MIN_DAPI_RNDtoRND_CORR= c(rep(filterCorr, times=NUM_GROUPS)),
    PARM_NORMALIZE_EVALUATE_PIPELINE=c(rep(value_PARM_NORMALIZE_EVALUATE_PIPELINE, times=NUM_GROUPS)),
    PARM_METRIC=c(rep(value_PARM_METRIC, times=NUM_GROUPS)),
    Display_X_Line_1= c(rep("", times=NUM_GROUPS)),
    Display_X_Line_2= c(rep("", times=NUM_GROUPS)))     
  
  maxVirtualSlides = c(2,3,10,14)
  pdTrials = unique(datFLINOresults[,"PARM_EVAL_DAPI_RNDS"])
  pdTrials = pdTrials[(str_count(pdTrials, pattern = ":") + 1) %in% maxVirtualSlides]
  pdResults = datFLINOresults
  pdResults = pdResults[is.na(pdResults[,"PARM_EVALUATE_POSITION"]), ]
  pdResults = pdResults[pdResults[,"PARM_EVAL_DAPI_RNDS"] %in% pdTrials,]
  
  par(fig=c(0,1,0,midpt), new=TRUE)
  par(mar = c(3, 5, 1, 1)) #  bottom, left, top and right  c(5.1, 4.1, 4.1, 2.1).
  # par(fig=c(0.5,1,0,1), new=TRUE)
  # par(mar = c(3, 3.5, 1, 1)) #  bottom, left, top and right  c(5.1, 4.1, 4.1, 2.1).
  
  setDefaultPlotParm()
  ADD_HORIZONTAL_BAR_AT_GIVEN_VALUE  = TRUE
  HORIZONTAL_BAR_VALUE = 0.0738
  pdTable = plotVirtualSlideData(pdResults, pdTrials, pdPlotColumns, 0 , 0.3, "")
  mtext("Nuclei Evaluation Error (CV)", line=3, side =2, cex=t2cex)
  text(paste0("Normalization ", get_display_PARM_METRIC(value_PARM_METRIC), " space ",  " - grid size ", gsub("Grid","",Target_GRID_SIZE)), adj=0, x = barCenters[1] / 2,y =0.26, cex = 0.8, xpd=NA)
  for(gi in c(1:NUM_GROUPS)) {
    text(format(pdTable[5,gi],digits = 3),  x = barCenters[gi],y =0.22, cex = 0.8, xpd=NA)
  }
  mtext("Normalization Method", side=1, line = 1, cex = 0.9)
  
  
  
  if (PLOT_TO != SCREEN) {
    dev.off()
  }
}



###############################################################
###############################################################
#PLOT_TO = SCREEN

if (GENERATE_PANEL_3) {
  panelName = "Fig_3"
  
  if (PLOT_TO == TIFF) {
    #jpeg(file.path(pathToFigures,paste0("Dapi",".slide",".jpeg")), quality = 100, width = 720, height = 480)
    myFile = file.path(pathToFigures,paste0(panelName,".tif"))
    tiff(myFile, height = 3.75, width = 4, units = 'in', compression="lzw",  type = "windows", res = 300)
    #par(cex=1.3)
  } else if (PLOT_TO == EPS) {
    myFile = file.path(pathToFigures,paste0(panelName,".eps"))
    setEPS()
    postscript(myFile,  height = 3.75, width = 45)
  }
  
  
  # par(fig=c(0,0.5,0,1), new=FALSE)
  # par(mar = c(3, 4.5, 1, 0)) #  bottom, left, top and right  c(5.1, 4.1, 4.1, 2.1).
  # 
  par(mar = c(3, 4.5, 1, 1))
  
  
  value_PARM_NORM_METHOD = "UPPER_QUARTILE"
  value_PARM_NORM_METHOD_LABEL = "Q75"
  Target_GRID_SIZE  = "Grid32" 
  filterCorr = 0
  value_PARM_METRIC = CV_METRIC_NormLogSpace
  NUM_GROUPS = 7
  pdPlotColumns = data.frame(
    PARM_NORM_METHOD= c(rep(value_PARM_NORM_METHOD, times=NUM_GROUPS)),
    PARM_EVAL_SEG_OBJ= c(rep(TRUE, times=NUM_GROUPS)),
    PARM_SEG_OBJ_NAME= c("NucleiSCA",rep(Target_GRID_SIZE, times=(NUM_GROUPS-1))),
    PARM_EVAL_SEG_OBJ_NAME= c(rep("NucleiSCA", times=NUM_GROUPS)),
    PARM_FILTER_MIN_DAPI_RNDtoRND_CORR= c(rep(filterCorr, times=NUM_GROUPS)),
    PARM_NORMALIZE_EVALUATE_PIPELINE=c("NORMALIZE_ALL_FOVS_EVALUATE_ALL_FOVS",
                                       "NORMALIZE_1_FOVS_EVALUATE_ALL_FOVS",
                                       "NORMALIZE_2_FOVS_EVALUATE_ALL_FOVS",
                                       "NORMALIZE_3_FOVS_EVALUATE_ALL_FOVS",
                                       "NORMALIZE_5_FOVS_EVALUATE_ALL_FOVS",
                                       "NORMALIZE_10_FOVS_EVALUATE_ALL_FOVS",
                                       "NORMALIZE_20_FOVS_EVALUATE_ALL_FOVS"),
    Display_X_Line_0=c("None",
                       "1",
                       "2",
                       "3",
                       "5",
                       "10",
                       "20"),
    PARM_METRIC= c(CV_METRIC_RawLogSpace,rep(CV_METRIC_NormLogSpace, times=(NUM_GROUPS-1))),
    Display_X_Line_1= c(rep("", times=NUM_GROUPS)),
    Display_X_Line_2= c(rep("", times=NUM_GROUPS)))    
  
  # par(mfrow=c(1,1))
  
  pdTrials = unique(datFLINOresults[,"PARM_EVAL_DAPI_RNDS"])
  #maxVirtualSlides = c(2,3,10,14)
  maxVirtualSlides = c(2,3)
  pdTrials = unique(datFLINOresults[,"PARM_EVAL_DAPI_RNDS"])
  pdTrials = pdTrials[(str_count(pdTrials, pattern = ":") + 1) %in% maxVirtualSlides]
  pdResults = datFLINOresults
  pdResults = pdResults[is.na(pdResults[,"PARM_EVALUATE_POSITION"]), ]
  pdResults = pdResults[pdResults[,"PARM_EVAL_DAPI_RNDS"] %in% pdTrials,]
  
  # par(fig=c(0.5,1,0,1), new=TRUE)
  # par(mar = c(3, 3.5, 1, 1)) #  bottom, left, top and right  c(5.1, 4.1, 4.1, 2.1).
  
  setDefaultPlotParm()
  ADD_HORIZONTAL_BAR_AT_GIVEN_VALUE  = TRUE
  HORIZONTAL_BAR_VALUE = 0.0738
  pdTable = plotVirtualSlideData(pdResults, pdTrials, pdPlotColumns,0, 2, "")
  mtext("Nuclei Evaluation Error (CV)", line=3, side =2, cex=t2cex)
  text(paste0(value_PARM_NORM_METHOD_LABEL," ", get_display_PARM_METRIC(value_PARM_METRIC), " space ","- grid size ", gsub("Grid","",Target_GRID_SIZE)), adj=0, x = barCenters[1] / 2,y =1.9, cex = 0.8, xpd=NA)
  for(gi in c(1:NUM_GROUPS)) {
    text(format(pdTable[5,gi],digits = 3),  x = barCenters[gi],y =1.75, cex = 0.8, xpd=NA)
  }
  mtext("Number of Control Samples Per Slide", side=1, line = 1, cex = 0.9)
  
  
  
  # filterTxt = paste0(Target_GRID_SIZE," Filter <", filterCorr, " Corr")
  # if (filterCorr == 0) {
  #   filterTxt = paste0(Target_GRID_SIZE," No Filtering")
  #   
  # }
  # correctionMethodTxt = paste0("Correction Method: ", value_PARM_NORM_METHOD, " log-norm ", filterTxt)
  # mtext(correctionMethodTxt, side=1, line = 2, cex = 0.9)
  # cat("\n", correctionMethodTxt,"\n")
  # print(pdTable)
  
  if (PLOT_TO != SCREEN) {
    dev.off()
  }
}

###############################################################
###############################################################





load(file.path(pathToData,"FOLFOX_MarkerData_Grid32.RData"), verbose = TRUE)
slides = unique(dat.dapi_SEG[,"slide"])  


selRoundNames = c("RG003", "RG005", "RG007", "RG008", "RG010", "RG012", "RG014", "RG015", "RG017", "RG018", "RG020", "RG021", "RG022", "RG024")
selRounds = as.numeric(gsub("RG","",selRoundNames))

# plotSlideExposureTime(slideDat, "Dapi", 150, "Dapi_exp_time")
plotSlideExposureTime = function(slideDat, channel, maxY, channelExposureLabel, pLegendLoc = "topleft") {
  slides = unique(slideDat[,"slide"])
  slides = slides[order(slides)]
  pcol = c("red","green","blue","black","orange")
  psym = c(20,21,22,23,24,25)
  db = slideDat[,c("Round",channelExposureLabel)]
  par(mar = c(5, 4, 1, 1)) #  bottom, left, top and right  c(5.1, 4.1, 4.1, 2.1).
  myXlab = paste0("Round - ", channel ," channel")
  myYlab="Exposure Time"
  plot(db[,"Round"],db[,1],pch=NA,xlim=c(0,25), ylim=c(0,maxY),ylab="",xlab="", panel.first = grid())
  for(i in 1:length(slides)) {
    db = slideDat[slideDat[,"slide"] == slides[i],c("Round",channelExposureLabel)]
    points(db[,"Round"],db[,channelExposureLabel],type='o',col=pcol[i],pch=psym[i])
  }
  legend(pLegendLoc,getSlideNames(slides),col=pcol,lty=1,pch = psym)
}



plotSlideDAPIRounds = function(slideDat, channel, maxY, roundInfo, pLegendLoc = "topleft") {
  
  db = slideDat
  
  pcol = c("red","green","blue","black","orange")
  psym = c(20,21,22,23,24,25)
  par(mar = c(5, 4, 1, 1)) #  bottom, left, top and right  c(5.1, 4.1, 4.1, 2.1).
  myXlab = paste0("Round - ", channel ," channel")
  myYlab="Slide Median (FOV Median (Grid Mean pixel intensity))"
  plot(roundInfo[,"Round"],db[,1],pch=NA,xlim=c(0,25), ylim=c(0,maxY),ylab="",xlab="", panel.first = grid())
  for(i in 1:ncol(db)) {
    points(roundInfo[,"Round"],db[,i],type='o',col=pcol[i],pch=psym[i])
  }
  legend(pLegendLoc,getSlideNames(colnames(db)),col=pcol,lty=1,pch = psym)
}



plotSlideRounds = function(slideDat, channel, maxY, roundInfo, pLegendLoc = "topleft") {
  
  
  db = slideDat[row.names(slideDat) %in% roundInfo[,"MarkerMeasureLabel"], ]
  db = db[match(roundInfo[,"MarkerMeasureLabel"],row.names(db)), ]
  
  pcol = c("red","green","blue","black","orange")
  psym = c(20,21,22,23,24,25)
  par(mar = c(5, 4, 1, 1)) #  bottom, left, top and right  c(5.1, 4.1, 4.1, 2.1).
  myXlab = paste0("Round - ", channel ," channel")
  myYlab="Slide Median (FOV Median (Grid Mean pixel intensity))"
  plot(roundInfo[,"Round"],db[,1],pch=NA,xlim=c(0,25), ylim=c(0,maxY),ylab="",xlab="", panel.first = grid())
  for(i in 1:ncol(db)) {
    points(roundInfo[,"Round"],db[,i],type='o',col=pcol[i],pch=psym[i])
  }
  for(i in 1:nrow(roundInfo)) {
    v = as.character(str_split_fixed(roundInfo[i,1],"\\.",4)[4])
    text(as.numeric(roundInfo[i,2]),1.1*as.numeric(db[i,1]),v, adj=0)
  }
  legend(pLegendLoc,getSlideNames(colnames(db)),col=pcol,lty=1,pch = psym)
}




plotSlideFOVDistributions = function(dat, round, channel, marker) {
  par(mar = c(5, 4, 1, 1)) #  bottom, left, top and right  c(5.1, 4.1, 4.1, 2.1).
  myXlab= paste0(marker, " Round ", round , " ", channel ," Channel")
  myYlab=paste0(marker, " (FOV Median (Grid Mean pixel intensity))")
  boxplot(marker~slide, data=dat, names=getSlideNames(slides), outpch = NA, xlab="",ylab="")
  for(si in 1:length(slides)) {
    stripchart(dat[dat[,"slide"] == slides[si],2],at =si, vertical = TRUE, method = "jitter", 
               pch = 20, col = "blue",
               add = TRUE,cex=1) 
  }
}




if (GENERATE_PANEL_S1) {
  
  ############ Figure S1 ############################################
  if (PLOT_TO == TIFF) {
    #jpeg(file.path(pathToFigures,paste0("Fig_S1_","Dapi",".slide",".jpeg")), quality = 100, width = 720, height = 480)
    myFile = file.path(pathToFigures,paste0("Fig_S1_","Dapi",".slide",".tif"))
    tiff(myFile, height = 4.75, width = 7.5, units = 'in', compression="lzw",  type = "windows", res = 300)
    par(cex=1.3)
  } else if (PLOT_TO == EPS) {
    myFile = file.path(pathToFigures,paste0(panelName,".eps"))
    setEPS()
    postscript(myFile,  height = 4.75, width = 7.5)
  }
  
  roundInfo = dapiRounds[dapiRounds[,"Round"] %in% selRounds, ]
  slideDat = tdat.marker_SEG.Slide[,c(1:3)]
  slideDat = slideDat[row.names(slideDat) %in% roundInfo[,"MarkerMeasureLabel"], ]
  plotSlideDAPIRounds(slideDat, "dapi", 2500, roundInfo)
  mtext("Imaging Round",side=1,line=3,cex = par("cex"))
  mtext("DAPI Staining Intensity",side=2,line=3,cex = par("cex"))
  if (PLOT_TO != SCREEN) {
    dev.off()
  }
}
####################################################################

if (GENERATE_PANEL_S2) {
  
  ############ Figure S2 ############################################
  if (PLOT_TO == TIFF) {
    #jpeg(file.path(pathToFigures,paste0("Fig_S2_","Exposure_Time." ,"Dapi", ".channel",".jpeg")), quality = 100, width = 720, height = 480)
    myFile = file.path(pathToFigures,paste0("Fig_S2_","Exposure_Time." ,"Dapi", ".channel",".tif"))
    tiff(myFile, height = 4.75, width = 7.5, units = 'in', compression="lzw",  type = "windows", res = 300)
    par(cex=1.3)
  } else if (PLOT_TO == EPS) {
    myFile = file.path(pathToFigures,paste0(panelName,".eps"))
    setEPS()
    postscript(myFile,  height = 4.75, width = 7.5)
  }
  slideDat = slide_round_channel_exposuretimes[slide_round_channel_exposuretimes[,"Round"] %in% selRounds, ]
  slideDat = slideDat[slideDat[,"slide"] %in% c("S18030274", "S18030275", "S18030276"), ]
  plotSlideExposureTime(slideDat, "Dapi", 150, "Dapi_exp_time")
  mtext("Imaging Round",side=1,line=3,cex = par("cex"))
  mtext("Dapi Channel Exposure Time",side=2,line=3,cex = par("cex"))
  if (PLOT_TO != SCREEN) {
    dev.off()
  }
}




###############################################################
###############################################################

###############################################################
###############################################################

# Build Figure 4

# Now Perform Grid-Based Normalization as an example

slides = c("S18030274", "S18030275", "S18030276")
cellObjectDataFile = file.path(pathToData, "FOLFOX_MarkerData_NucleiSCA.RData")
gridObjectDataFile = file.path(pathToData,  "FOLFOX_MarkerData_Grid32.RData")

applyGridBased_FLINO(cellObjectDataFile = cellObjectDataFile, gridObjectDataFile = gridObjectDataFile, slides = slides)

##################################################################
##################################################################


if (GENERATE_PANEL_4) {
  panelName = "Fig_4"
  
  fig_width = 7    # inches
  fig_height = 6   # inches
  
  fig_width = 6.3    # inches
  fig_height = 5.4   # inches
  
  if (PLOT_TO == TIFF) {
    #jpeg(file.path(pathToFigures,paste0("Dapi",".slide",".jpeg")), quality = 100, width = 720, height = 480)
    myFile = file.path(pathToFigures,paste0(panelName,".tif"))
    tiff(myFile, height = fig_height, width = fig_width, units = 'in', compression="lzw",  type = "windows", res = 300)
    #par(cex=1.3)
  } else if (PLOT_TO == EPS) {
    myFile = file.path(pathToFigures,paste0(panelName,".eps"))
    setEPS()
    postscript(myFile,  height = fig_height, width = fig_width)
  }
  
  
  n = layout(mat = rbind(c(9, 11, 10, 10), 
                         c(1, 11, 3, 4), 
                         c(1, 11, 5, 6), 
                         c(2, 11, 5, 6), 
                         c(2, 11, 7, 8)),
             widths =c(1.3, 0.2, 1, 1), 
             heights =  c(0.05, 1, 0.5, 0.5 , 1), respect = FALSE)
  
  # layout.show(n)
  
  cellPositions = c(81,82,83,84)
  positions = unique(nuclei.BMdat.FOV[,"position"])
  # Position - Cell Lines
  # 81 – HCT116_XIAP KO
  # 82 - JURKAT
  # 83 - HeLa
  # 84 – MCF7
  
  
  setDefaultPlotParm()
  DISPLAY_DATA_LINE_SEGMENTS = TRUE
  symSize = 1
  pd_bm = "Mean.AF000.BAX.BAX" 
  
  pd_ylim = NA
  aggStr = "FOV"
  
  dat = nuclei.BMdat.gridNorm.FOV
  pd_dat = dat[dat[,"position"] %in% positions, c("slide","position",pd_bm)]
  colnames(pd_dat) = c("slide","position","marker")
  pd_dat[,"marker"] = log2(pd_dat[,"marker"]+1)
  pd_ylim1 = c(min( pd_dat[,"marker"],na.rm = TRUE),max( pd_dat[,"marker"],na.rm = TRUE))
  
  dat = nuclei.BMdat.FOV
  pd_dat = dat[dat[,"position"] %in% positions, c("slide","position",pd_bm)]
  colnames(pd_dat) = c("slide","position","marker")
  pd_dat[,"marker"] = log2(pd_dat[,"marker"]+1)
  pd_ylim3 = c(min( pd_dat[,"marker"],na.rm = TRUE),max( pd_dat[,"marker"],na.rm = TRUE))
  
  pd_ylim = c(max(6,min(c(pd_ylim1[1],pd_ylim3[1]))),max(c(pd_ylim1[2],pd_ylim3[2])))
  

  
  dat = nuclei.BMdat.FOV
  pd_dat = dat[dat[,"position"] %in% positions, c("slide","position",pd_bm)]
  colnames(pd_dat) = c("slide","position","marker")
  pd_dat[,"marker"] = log2(pd_dat[,"marker"]+1)
  par(mar = c(3, 5, 3, 0)) #  bottom, left, top and right 
  plotBMdata(pd_dat, pd_bm, pd_ylim, cellPositions)
  mtext("BAX Image Intensity (Log2)",side=2,line=3,cex = 1)
  mtext("Uncorrected",side=3,line=0,cex = 1.2, font=2)
  text(x=3.15,y=pd_dat[pd_dat[,"slide"] == "S18030276" & pd_dat[,"position"] == 81, "marker"], adj=0,"HCT116", xpd=NA, cex=1.2, font=2)
  text(x=3.15,y=pd_dat[pd_dat[,"slide"] == "S18030276" & pd_dat[,"position"] == 82, "marker"], adj=0,"JURKAT", xpd=NA, cex=1.2, font=2)
  text(x=3.45,y=pd_dat[pd_dat[,"slide"] == "S18030276" & pd_dat[,"position"] == 83, "marker"], adj=0,"HeLa", xpd=NA, cex=1.4, font=2)
  text(x=3.15,y=pd_dat[pd_dat[,"slide"] == "S18030276" & pd_dat[,"position"] == 84, "marker"], adj=0,"MCF7", xpd=NA, cex=1.2, font=2)
  
  
  dat = nuclei.BMdat.gridNorm.FOV
  pd_dat = dat[dat[,"position"] %in% positions, c("slide","position",pd_bm)]
  colnames(pd_dat) = c("slide","position","marker")
  pd_dat[,"marker"] = log2(pd_dat[,"marker"]+1)
  par(mar = c(3, 5, 3, 0)) #  bottom, left, top and right 
  plotBMdata(pd_dat, pd_bm, pd_ylim, cellPositions)
  mtext("BAX Image Intensity (Log2)",side=2,line=3,cex = 1)
  mtext("Normalized",side=3,line=0,cex = 1.2, font=2)
  text(x=3.15,y=pd_dat[pd_dat[,"slide"] == "S18030276" & pd_dat[,"position"] == 81, "marker"], adj=0,"HCT116", xpd=NA, cex=1.2, font=2)
  text(x=3.15,y=pd_dat[pd_dat[,"slide"] == "S18030276" & pd_dat[,"position"] == 82, "marker"], adj=0,"JURKAT", xpd=NA, cex=1.2, font=2)
  text(x=3.45,y=pd_dat[pd_dat[,"slide"] == "S18030276" & pd_dat[,"position"] == 83, "marker"], adj=0,"HeLa", xpd=NA, cex=1.4, font=2)
  text(x=3.15,y=pd_dat[pd_dat[,"slide"] == "S18030276" & pd_dat[,"position"] == 84, "marker"], adj=0,"MCF7", xpd=NA, cex=1.2, font=2)
  
  
  for(fn in c("BAX_HeLa_Slide_A1_Uncorrected.png",
          "BAX_HeLa_Slide_A1_Normalized.png",
           "BAX_HeLa_Slide_A2_Uncorrected.png",
           "BAX_HeLa_Slide_A2_Normalized.png",
           "BAX_HeLa_Slide_A3_Uncorrected.png",
           "BAX_HeLa_Slide_A3_Normalized.png" )) {
    par(mar = c(0, 0, 0, 0)) #  bottom, left, top and right 
    image_fn =  file.path(pathToImages,fn)
    pp <- readPNG(image_fn)
    x0 = 0
    y0 = 0
    plot(c(0, ncol(pp) * 1), c(0, nrow(pp) * 1), frame.plot = FALSE  ,type = "n", xaxt='n', yaxt='n', xlab = "", ylab = "") #panel.first = grid())
    rasterImage(pp, x0, y0, x0 + ncol(pp), y0 + nrow(pp))
    
    if (length(grep("_Uncorrected.png",fn)) > 0) {
      slideName = gsub("BAX_HeLa_Slide_","",gsub("_Uncorrected.png","",gsub("_Normalized.png","",fn)))
      text(10,50, paste0("Slide ", slideName),col='white', adj=0,cex=1.5, font=2)
      
    }
    
    if (fn == "BAX_HeLa_Slide_A1_Uncorrected.png") {
      text(20,nrow(pp)-75, "Uncorrected\nHeLa cell Images",col='white', adj=0,cex=1.5, font=2)
    }
    
    if (fn == "BAX_HeLa_Slide_A1_Normalized.png") {
      text(20,nrow(pp)-50, "Normalized",col='white', adj=0,cex=1.5, font=2)
    }

  }


  par(mar=c(0,0,0,0)) # par(mar = c(bottom, left, top, right))
  plot.new()
  #rect(0,0,1,1)
  mtext("A", side=3, at=0.05, outer=FALSE, line=-2, cex = 1.5, font=2, xpd = NA)
  #mtext("BAX statining intensity of four cell\nlines across three physical slides", side=3, adj=0, at=0.15, outer=FALSE, line=-3, cex = 1, font=1, xpd = NA)
   
  
  par(mar=c(0,0,0,0)) # par(mar = c(bottom, left, top, right))
  plot.new()
  #rect(0,0,1,1)
  mtext("B", side=3, at=-0.1, outer=FALSE, line=-2, cex = 1.5, font=2, xpd = NA)
  # mtext("BAX images of HeLa cells across three slides", side=3, adj=0, at=0.05, outer=FALSE, line=-2, cex = 1, font=1, xpd = NA)
  # 
  # mtext("Uncorrected",side=3, at=0.2, line=-5,cex = 1.1, font=2)
  # mtext("Normalized",side=3, at=0.7,line=-5,cex = 1.1, font=2)
  

  if (PLOT_TO != SCREEN) {
    dev.off()
  }
}

#par(mfrow=c(1,1))



####################################################################
####################################################################







