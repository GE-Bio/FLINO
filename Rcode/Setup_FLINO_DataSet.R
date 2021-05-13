options(stringsAsFactors=FALSE)

pathToData = "Data"

REMOVE_BUILD_FILES = FALSE

if (!file.exists(file.path(pathToData,"FOLFOX_MarkerData_Grid32.RData"))) {
  
  cat("Loading...","\n")
  cat("FOLFOX_MarkerData_Grid32_S18030274.RData","\n")
  load(file.path(pathToData,"FOLFOX_MarkerData_Grid32_S18030274.RData"), verbose = FALSE)
  dat = dat.marker_SEG
  cat("FOLFOX_MarkerData_Grid32_S18030275.RData","\n")
  load(file.path(pathToData,"FOLFOX_MarkerData_Grid32_S18030275.RData"), verbose = FALSE)
  dat = rbind(dat,dat.marker_SEG)
  cat("FOLFOX_MarkerData_Grid32_S18030276.RData","\n")
  load(file.path(pathToData,"FOLFOX_MarkerData_Grid32_S18030276.RData"), verbose = FALSE)
  dat.marker_SEG = rbind(dat,dat.marker_SEG)
  rm(dat)
  fn.out = "FOLFOX_MarkerData_Grid32.RData"
  cat("Saving as...","\n")
  cat(fn.out,"\n")
  save(   cy3rounds,
          cy5rounds,
          dapiRounds,
          dat.marker_SEG,
          dat.marker_SEG.FOV,
          dat.marker_SEG.Slide,
          tdat.marker_SEG.Slide,
          slide_position_patient_id,
          slide_round_channel_exposuretimes,
          rounds,
          slides,
          file=file.path(pathToData,fn.out))

  if (REMOVE_BUILD_FILES & file.exists(file.path(pathToData,"FOLFOX_MarkerData_Grid32.RData"))) {
    #Delete file if it exists
    file.remove("FOLFOX_MarkerData_Grid32_S18030274.RData")
    file.remove("FOLFOX_MarkerData_Grid32_S18030275.RData")
    file.remove("FOLFOX_MarkerData_Grid32_S18030276.RData")
  }
  
}

if (!file.exists(file.path(pathToData,"FOLFOX_MarkerData_NucleiSCA.RData"))) {
    
  cat("Loading...","\n")
  cat("FOLFOX_MarkerData_NucleiSCA_S18030274.RData","\n")
  load(file.path(pathToData,"FOLFOX_MarkerData_NucleiSCA_S18030274.RData"), verbose = FALSE)
  dat = dat.marker_SEG
  cat("FOLFOX_MarkerData_NucleiSCA_S18030275.RData","\n")
  load(file.path(pathToData,"FOLFOX_MarkerData_NucleiSCA_S18030275.RData"), verbose = FALSE)
  dat = rbind(dat,dat.marker_SEG)
  cat("FOLFOX_MarkerData_NucleiSCA_S18030276.RData","\n")
  load(file.path(pathToData,"FOLFOX_MarkerData_NucleiSCA_S18030276.RData"), verbose = FALSE)
  dat.marker_SEG = rbind(dat,dat.marker_SEG)
  rm(dat)
  fn.out = "FOLFOX_MarkerData_NucleiSCA.RData"
  cat("Saving as...","\n")
  cat(fn.out,"\n")
  save(   cy3rounds,
          cy5rounds,
          dapiRounds,
          dat.marker_SEG,
          dat.marker_SEG.FOV,
          dat.marker_SEG.Slide,
          tdat.marker_SEG.Slide,
          slide_position_patient_id,
          slide_round_channel_exposuretimes,
          rounds,
          slides,
          file=file.path(pathToData,fn.out))
  
  if (REMOVE_BUILD_FILES & file.exists(file.path(pathToData,"FOLFOX_MarkerData_NucleiSCA.RData"))) {
    #Delete file if it exists
    file.remove("FOLFOX_MarkerData_NucleiSCA_S18030274.RData")
    file.remove("FOLFOX_MarkerData_NucleiSCA_S18030275.RData")
    file.remove("FOLFOX_MarkerData_NucleiSCA_S18030276.RData")
  }
}

