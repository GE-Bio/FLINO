####################################################
##### Run_FLINO_Evaluator ##########################
####################################################

if (!"analysisRun" %in% objects()) {

  # Example run that can be used as a template
  analysisRun = "eRuns_Grid256_Q75NZ_14VS.txt"

}
  
pathToData = "Data"
pathToRcode= "Rcode" 
pathToResults = "Results"
pathToeRuns = "eRuns"

dir.create(file.path(pathToResults), showWarnings = FALSE)

source(file.path(pathToRcode,"Setup_FLINO_DataSet.R"))

source(file.path(pathToRcode,"FLINO_Functions.R"))

fnInputRuns = file.path(pathToeRuns, analysisRun) 
fnResults = file.path(pathToResults, paste0("results_", analysisRun))
  


if (!"displayFits" %in% objects()) {
  displayFits = FALSE
}

if ("RND_SEED_VALUE" %in% objects()) {
  cat("Setting RND_SEED_VALUE to: ", RND_SEED_VALUE,"\n")
  set.seed(RND_SEED_VALUE)
}

if (!"scaleIndexes" %in% objects()) {
  scaleIndexes = c(2,1)
}


eRuns = read.delim(file=fnInputRuns, header = TRUE, sep = "\t", blank.lines.skip = TRUE, na.strings = c("NA","NaN"))
eRuns = eRuns[eRuns[,"PARM_RUN"] == TRUE, ]



#########################################################

library("plyr")


computeGridRndToRndVariation = function(db) {

    d1 = aggregate(x = db[,bm], by = list(db[,"position"],db[,"Cell.ID"]), FUN = "cvRemoveNAs")
    colnames(d1) = c("position","Cell.ID",paste0("cv.",bm))
  
    d2 = aggregate(x = db[,bm], by = list(db[,"position"],db[,"Cell.ID"]), FUN = "countNonNAs")
    colnames(d2) = c("position","Cell.ID",paste0("n.",bm))
  
    #d <- merge(d1,d2,by=c("position","Cell.ID"))
    d <- join(d1,d2,by=c("position","Cell.ID"))
  
    return(d)
  
}

#########################################################

getNumToNormalize = function(x) {
  x = gsub("NORMALIZE_","",x, ignore.case = TRUE)
  x = gsub("_OBJS_EVALUATE_ALL_FOVS","",x, ignore.case = TRUE)
  x = gsub("_FOVS_EVALUATE_ALL_FOVS","",x, ignore.case = TRUE)
  if (x == "ALL") {
    return(0)
  } else {
    x = as.numeric(x)
    return(x)
  }
}

#########################################################

results = data.frame()
printResults = function(elistNum, x, results) {


  
  xNames = names(x)
  ri = nrow(results) + 1
  for(cn in c("PARM_EVALUATE_SLIDE",
             # "PARM_EVALUATE_POSITION",
              "MinNumPerSlide_NormSegObj",
              "MeanNumPerSlide_NormSegObj",
              "MedianNumPerSlide_NormSegObj",
              "MaxNumPerSlide_NormSegObj",
              "SegObjErrCV_RawAbsSpace",
              "SegObjErrCV_NormAbsSpace",
              "SegObjErrCV_RawLogSpace",
              "SegObjErrCV_NormLogSpace",
              "Number_NormSegObj",
              "PARM_FILTER_MIN_DAPI_RNDtoRND_CORR",
              "PARM_NORM_METHOD",
              "PARM_SEG_OBJ_NAME",
              "PARM_NORMALIZE_EVALUATE_PIPELINE",
              "PARM_EVAL_SEG_OBJ_NAME",
              "EvalSegObj_SegObjErrCV_RawAbsSpace",
              "EvalSegObj_SegObjErrCV_NormAbsSpace",
              "EvalSegObj_SegObjErrCV_RawLogSpace",
              "EvalSegObj_SegObjErrCV_NormLogSpace",
              "Number_EvalSegObj")) {
    if (cn %in% xNames) {
      results[ri,cn] = x[[cn]]
    } else {
      results[ri,cn] = NA
    }
  }
  return(results)
}


EvalSegObj_PARM_EVAL_SEG_OBJ_NAME = ""
EvalSegObj_PARM_EVALUATE_SLIDE = ""

NormSegObj_PARM_SEG_OBJ_NAME = ""
NormSegObj_PARM_EVALUATE_SLIDE = ""

for (elistNum in 1:getNumEvaluationLists()) {

  
  cat("Performing evaluation run ",elistNum, " of ", getNumEvaluationLists(),"\n")
  
  elist = getEvaluationList(elistNum)
  
  
  ####*******************************************************************
  ####******* START of load and initialize EvalSegObj  ******************
  ####*******************************************************************
  if (elist[["PARM_EVAL_SEG_OBJ"]]) {

    if (elist[["PARM_EVAL_SEG_OBJ_NAME"]] != EvalSegObj_PARM_EVAL_SEG_OBJ_NAME ||
        elist[["PARM_EVALUATE_SLIDE"]] != EvalSegObj_PARM_EVALUATE_SLIDE) {
      
      
      cat("Loading...","\n")
      cat(file.path(pathToData,elist[["PARM_EVAL_SEG_OBJ_FILE_NAME"]]),"\n")
      load(file.path(pathToData,elist[["PARM_EVAL_SEG_OBJ_FILE_NAME"]]), verbose = FALSE)
      corrFN = gsub(".RData","_Corr.RData", file.path(pathToData,elist[["PARM_EVAL_SEG_OBJ_FILE_NAME"]]))
      if (file.exists(corrFN)) {
        load(corrFN, verbose = FALSE)
      }
      
      cat("Initializing EvalSegObj...","\n")
      
      roundInfo = dapiRounds[dapiRounds[,"Round"] %in% unique(c(cy3rounds[,"Round"],cy5rounds[,"Round"])),]
      roundInfo = roundInfo[order(roundInfo[,"Round"]), ]
      row.names(roundInfo) = NULL
      
      bm = roundInfo[,"MarkerMeasureLabel"]
      dat0_AllSlides = get_dat0(dat.dapi_SEG, dat.Corr_SEG, dapiRounds, dapiRounds, bm) #, minCorrValue, minIntensityValue)
      
      roundInfo[,"SlideLabel"] = gsub(".dapi.dapi","",gsub("Mean.","",roundInfo[,"MarkerMeasureLabel"]), ignore.case = TRUE)
      
      dat0_slide = dat0_AllSlides[dat0_AllSlides[,"slide"] == elist[["PARM_EVALUATE_SLIDE"]], ]

      cns = colnames(dat0_slide)
      dat0_AllPos = data.frame()
      for(ri in 1:nrow(roundInfo)) {
        temp = dat0_slide[,c(cns[!startsWith(cns,"Mean.")],gsub(".dapi.dapi",".Corr.Corr",roundInfo[ri,"MarkerMeasureLabel"]),roundInfo[ri,"MarkerMeasureLabel"])]
        colnames(temp)[(ncol(temp)-1)] = "dapi_RndToRnd_Corr"
        colnames(temp)[ncol(temp)] = "Mean.AF000.dapi.dapi"
        temp[,"slide"] = roundInfo[ri,"SlideLabel"]
        if (nrow(dat0_AllPos) == 0) {
          dat0_AllPos = temp
        } else {
          dat0_AllPos = rbind(dat0_AllPos,temp)
        }
      }
      rm(temp)
      dat0_AllPos[,"Mean.AF000.dapi.dapi_NormInLogSpace"] = dat0_AllPos[,"Mean.AF000.dapi.dapi"]
      
      
      
      
      rm(cy3rounds)
      rm(cy5rounds)
      rm(dapiRounds)
      rm(dat.Corr_SEG)
      rm(dat.dapi_SEG)
      rm(dat.dapi_SEG.FOV)
      rm(dat.dapi_SEG.Slide)
      rm(tdat.dapi_SEG.Slide)
      rm(slide_position_patient_id)
      rm(slide_round_channel_exposuretimes)
      rm(rounds)
      rm(slides)
      
      rm(roundInfo)
      rm(dat0_AllSlides)
      rm(dat0_slide)
      
      EvalSegObj_dat0_AllPos = dat0_AllPos
      EvalSegObj_PARM_SEG_OBJ_NAME = elist[["PARM_SEG_OBJ_NAME"]]
      EvalSegObj_PARM_EVALUATE_SLIDE =  elist[["PARM_EVALUATE_SLIDE"]]
      
      rm(dat0_AllPos)
      
    } else {
      cat("Reinitializing EvalSegObj...","\n")
    }
    
    dat0 = EvalSegObj_dat0_AllPos
    
    if (!is.null(elist[["PARM_EVALUATE_POSITION"]]) & !is.na(elist[["PARM_EVALUATE_POSITION"]][1])) {
      dat0 = dat0[dat0[,"position"] %in% elist[["PARM_EVALUATE_POSITION"]], ]
    }
    
    ####################################################################
    ####################################################################
    # remove virtual slides that will not be used (i.e. dapi rounds)
    dat0 = dat0[dat0[,"slide"] %in% elist[["PARM_EVAL_DAPI_RNDS"]], ]
    # this is important to prevent problems
    if (elist[["PARM_EVAL_MIN_NUM_DAPI_RNDS"]] > length(elist[["PARM_EVAL_DAPI_RNDS"]])) {
      elist[["PARM_EVAL_MIN_NUM_DAPI_RNDS"]] = length(elist[["PARM_EVAL_DAPI_RNDS"]])
    }
    ####################################################################
    ####################################################################
    
    # **********************************************************************************************
    # **********************************************************************************************
    vi = which(!is.na(dat0[,"Mean.AF000.dapi.dapi"]) & as.numeric(dat0[,"Mean.AF000.dapi.dapi"]) >=  elist[["PARM_EVAL_MIN_DAPI_INTENSITY"]]  
               & !is.na(dat0[,"dapi_RndToRnd_Corr"])  & as.numeric(dat0[,"dapi_RndToRnd_Corr"]) >= elist[["PARM_EVAL_MIN_DAPI_RNDtoRND_CORR"]] 
               & !is.na(dat0[,"Area"])  & as.numeric(dat0[,"Area"]) >= elist[["PARM_EVAL_MIN_OBJ_AREA"]]  & as.numeric(dat0[,"Area"]) <= elist[["PARM_EVAL_MAX_OBJ_AREA"]])
    
    if (length(vi) < 1) {
      cat("xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx","\n")
      cat("RUN ABORTED: Filtering of EVAL_SEG_OBJ caused vi==0 skipping run","\n")
      cat("xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx","\n")
      next
    }
    
    # > quantile(EvalSegObj_dat0[,"Area"],probs = c(0,0.05,0.1,0.25,0.5,0.75,0.9,0.95,1))
    # 0%    5%   10%   25%   50%   75%   90%   95%  100% 
    # 62    73    85   123   203   320   471   593 14500 
    
    bm = c("Mean.AF000.dapi.dapi","Mean.AF000.dapi.dapi_NormInLogSpace")
    BMdat=dat0[,bm]
    BMdat[vi,]= apply(BMdat[vi,],2,pmax,1) #forces each element to be >= to pixel minimum value of 1
    # **********************************************************************************************
    # **********************************************************************************************
    
    
    EvalSegObj_dat0 = dat0
    EvalSegObj_BMdat = BMdat
    EvalSegObj_vi = vi
    
    
    rm(dat0)
    rm(BMdat)
    rm(vi)
  }
  ####*******************************************************************
  ####******* END of load and initialize EvalSegObj  ********************
  ####*******************************************************************
  
  if (elist[["PARM_SEG_OBJ_NAME"]] != NormSegObj_PARM_SEG_OBJ_NAME ||
      elist[["PARM_EVALUATE_SLIDE"]] != NormSegObj_PARM_EVALUATE_SLIDE) {
    
    
    cat("Loading...","\n")
    cat(file.path(pathToData,elist[["PARM_SEG_OBJ_FILE_NAME"]]),"\n")
    load(file.path(pathToData,elist[["PARM_SEG_OBJ_FILE_NAME"]]), verbose = FALSE)
    corrFN = gsub(".RData","_Corr.RData", file.path(pathToData,elist[["PARM_SEG_OBJ_FILE_NAME"]]))
    if (file.exists(corrFN)) {
      load(corrFN, verbose = FALSE)
    }
    cat("Initializing...","\n")
  
    roundInfo = dapiRounds[dapiRounds[,"Round"] %in% unique(c(cy3rounds[,"Round"],cy5rounds[,"Round"])),]
    roundInfo = roundInfo[order(roundInfo[,"Round"]), ]
    row.names(roundInfo) = NULL
    
    bm = roundInfo[,"MarkerMeasureLabel"]

    dat0_AllSlides = get_dat0(dat.dapi_SEG, dat.Corr_SEG, dapiRounds, dapiRounds, bm)
    
    roundInfo[,"SlideLabel"] = gsub(".dapi.dapi","",gsub("Mean.","",roundInfo[,"MarkerMeasureLabel"]), ignore.case = TRUE)

    dat0_slide = dat0_AllSlides[dat0_AllSlides[,"slide"] == elist[["PARM_EVALUATE_SLIDE"]], ]
    
    cns = colnames(dat0_slide)
    dat0_AllPos = data.frame()
    for(ri in 1:nrow(roundInfo)) {
      temp = dat0_slide[,c(cns[!startsWith(cns,"Mean.")],gsub(".dapi.dapi",".Corr.Corr",roundInfo[ri,"MarkerMeasureLabel"]),roundInfo[ri,"MarkerMeasureLabel"])]
      colnames(temp)[(ncol(temp)-1)] = "dapi_RndToRnd_Corr"
      colnames(temp)[ncol(temp)] = "Mean.AF000.dapi.dapi"
      temp[,"slide"] = roundInfo[ri,"SlideLabel"]
      if (nrow(dat0_AllPos) == 0) {
        dat0_AllPos = temp
      } else {
        dat0_AllPos = rbind(dat0_AllPos,temp)
      }
    }
    rm(temp)
    dat0_AllPos[,"Mean.AF000.dapi.dapi_NormInLogSpace"] = dat0_AllPos[,"Mean.AF000.dapi.dapi"]
    
    
    rm(cy3rounds)
    rm(cy5rounds)
    rm(dapiRounds)
    rm(dat.Corr_SEG)
    rm(dat.dapi_SEG)
    rm(dat.dapi_SEG.FOV)
    rm(dat.dapi_SEG.Slide)
    rm(tdat.dapi_SEG.Slide)
    rm(slide_position_patient_id)
    rm(slide_round_channel_exposuretimes)
    rm(rounds)
    rm(slides)
    
    rm(roundInfo)
    rm(dat0_AllSlides)
    rm(dat0_slide)
    
    NormSegObj_dat0_AllPos = dat0_AllPos
    NormSegObj_PARM_SEG_OBJ_NAME = elist[["PARM_SEG_OBJ_NAME"]]
    NormSegObj_PARM_EVALUATE_SLIDE =  elist[["PARM_EVALUATE_SLIDE"]]
    
    rm(dat0_AllPos)


  } else {
    cat("Reinitializing...","\n")
  }
  

  dat0 = NormSegObj_dat0_AllPos

  if (!is.null(elist[["PARM_EVALUATE_POSITION"]]) & !is.na(elist[["PARM_EVALUATE_POSITION"]][1])) {
    dat0 = dat0[dat0[,"position"] %in% elist[["PARM_EVALUATE_POSITION"]], ]
  }
  
  ####################################################################
  ####################################################################
  # remove virtual slides that will not be used (i.e. dapi rounds)
  dat0 = dat0[dat0[,"slide"] %in% elist[["PARM_EVAL_DAPI_RNDS"]], ]
  # this is important to prevent problems
  if (elist[["PARM_EVAL_MIN_NUM_DAPI_RNDS"]] > length(elist[["PARM_EVAL_DAPI_RNDS"]])) {
    elist[["PARM_EVAL_MIN_NUM_DAPI_RNDS"]] = length(elist[["PARM_EVAL_DAPI_RNDS"]])
  }
  ####################################################################
  ####################################################################

  
  # **********************************************************************************************
  # **********************************************************************************************
  
  vi = which(!is.na(dat0[,"Mean.AF000.dapi.dapi"]) & as.numeric(dat0[,"Mean.AF000.dapi.dapi"]) >=  elist[["PARM_FILTER_MIN_DAPI_INTENSITY"]]  
             & !is.na(dat0[,"dapi_RndToRnd_Corr"])  & as.numeric(dat0[,"dapi_RndToRnd_Corr"]) >= elist[["PARM_FILTER_MIN_DAPI_RNDtoRND_CORR"]] 
             & !is.na(dat0[,"Area"])  & as.numeric(dat0[,"Area"]) >= elist[["PARM_FILTER_MIN_OBJ_AREA"]]  & as.numeric(dat0[,"Area"]) <= elist[["PARM_FILTER_MAX_OBJ_AREA"]])
  
  if (length(vi) < 1) {
    cat("xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx","\n")
    cat("RUN ABORTED: Filtering caused vi==0 skipping run","\n")
    cat("xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx","\n")
    next
  }
  
  bm = c("Mean.AF000.dapi.dapi","Mean.AF000.dapi.dapi_NormInLogSpace")
  BMdat=dat0[,bm]
  BMdat[vi,]= apply(BMdat[vi,],2,pmax,1) #forces each element to be >= to pixel minimum value of 1
  # **********************************************************************************************
  # **********************************************************************************************
 
  
  if (elist[["PARM_NORMALIZE_EVALUATE_PIPELINE"]] == "NORMALIZE_SAMPLE_OBJ_EVALUATE_ALL_FOVS") {
    
    #randomly select N OBJ's for each virtual slide with no repeats
    
    vs = unique(dat0[vi,"slide"])
    nvs = length(vs)
    # only consider Objects that are found on all virtual slides
    d = unique(dat0[vi,c("slide","position","Cell.ID")])               
    d = aggregate(d[,"slide"],by=list(d[,"position"], d[,"Cell.ID"]),FUN=length)
    colnames(d) = c("position","Cell.ID","x")
    virtual_objects = d[d[,"x"] >= nvs, c("position","Cell.ID") ]
    virtual_objects[,"VALID_VIRTUAL_OBJ"] = TRUE
    #dat0 = merge(dat0,virtual_objects,by=c("position","Cell.ID"),all.x=TRUE)
    dat0 = join(dat0,virtual_objects,by=c("position","Cell.ID"))
    
    numVirtualOBJs = floor(nrow(virtual_objects) / nvs)
    
    sampledSlideObjs  = cbind(data.frame(slide=rep(vs,times=numVirtualOBJs)), virtual_objects[sample(c(1:nrow(virtual_objects)), size=numVirtualOBJs * nvs, replace =F),])
    sampledSlideObjs[,"USE_TO_NORMALIZE"] = TRUE
    #dat0 = merge(dat0,sampledSlideObjs[,c("slide","position","Cell.ID","USE_TO_NORMALIZE")],by=c("slide","position","Cell.ID"),all.x=TRUE)
    dat0 = join(dat0,sampledSlideObjs[,c("slide","position","Cell.ID","USE_TO_NORMALIZE")],by=c("slide","position","Cell.ID"))
    
    viNorm = which(!is.na(dat0[,"USE_TO_NORMALIZE"]) & dat0[,"USE_TO_NORMALIZE"] == TRUE)
    
  } else {
    if (length(grep("_OBJS_", elist[["PARM_NORMALIZE_EVALUATE_PIPELINE"]])) > 0) {
      cat("Select OBJs for Normalization...","\n")
      numVirtualOBJs = getNumToNormalize(elist[["PARM_NORMALIZE_EVALUATE_PIPELINE"]])
      numVirtualFOVs = 0
    } else {
      cat("Select FOVs for Normalization...","\n")
      numVirtualOBJs = 0
      numVirtualFOVs = getNumToNormalize(elist[["PARM_NORMALIZE_EVALUATE_PIPELINE"]])
    }
    
    
    
    if (numVirtualFOVs == 0 & numVirtualOBJs == 0) {
      
      dat0[vi,"VALID_VIRTUAL_POSITION"] = TRUE
      dat0[vi,"USE_TO_NORMALIZE"] = TRUE
      
      viNorm = vi
      
    } else if (numVirtualFOVs > 0) {
      #randomly select N FOV's for each virtual slide with no repeats
      
      vs = unique(dat0[vi,"slide"])
      nvs = length(vs)
      # only consider positions that are found on all virtual slides
      d = unique(dat0[vi,c("slide","position")])               
      d = aggregate(d[,"slide"],by=list(d[,"position"]),FUN=length)
      virtual_positions = d[d[,"x"] >= nvs, 1]
      if (length(virtual_positions) < numVirtualFOVs  * nvs) {
        cat("xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx","\n")
        cat("RUN ABORTED: virtual_positions < numVirtualFOVs  * nvs : ", "\t", length(virtual_positions), "\t",numVirtualFOVs, "\t", nvs,"\n")
        cat("xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx","\n")
        next
      }
      dat0[vi,"VALID_VIRTUAL_POSITION"] = FALSE
      dat0[!is.na(dat0[,"VALID_VIRTUAL_POSITION"]) & dat0[,"position"] %in% virtual_positions,"VALID_VIRTUAL_POSITION"] = TRUE

      
      sampledSlidePos = data.frame(slide=rep(vs,times=numVirtualFOVs),position=sample(virtual_positions, size=numVirtualFOVs * nvs, replace =F))
      sampledSlidePos = sampledSlidePos[order(sampledSlidePos[,"slide"], as.numeric(sampledSlidePos[,"position"])),]
      

      
      dat0[vi,"USE_TO_NORMALIZE"] = FALSE
      for(slide in unique(sampledSlidePos[,"slide"])) {
        dat0[!is.na(dat0[,"USE_TO_NORMALIZE"]) & dat0[,"slide"] == slide & dat0[,"position"] %in% sampledSlidePos[sampledSlidePos[,"slide"] == slide, "position"],"USE_TO_NORMALIZE"] = TRUE
      }
      viNorm = which(!is.na(dat0[,"USE_TO_NORMALIZE"]) & dat0[,"USE_TO_NORMALIZE"] == TRUE)
      
    } else if (numVirtualOBJs > 0) {
      #randomly select N OBJ's for each virtual slide with no repeats
      
      vs = unique(dat0[vi,"slide"])
      nvs = length(vs)
      # only consider Objects that are found on all virtual slides
      d = unique(dat0[vi,c("slide","position","Cell.ID")])               
      d = aggregate(d[,"slide"],by=list(d[,"position"], d[,"Cell.ID"]),FUN=length)
      colnames(d) = c("position","Cell.ID","x")
      virtual_objects = d[d[,"x"] >= nvs, c("position","Cell.ID") ]
      virtual_objects[,"VALID_VIRTUAL_OBJ"] = TRUE
      #dat0 = merge(dat0,virtual_objects,by=c("position","Cell.ID"),all.x=TRUE)
      dat0 = join(dat0,virtual_objects,by=c("position","Cell.ID"))
      
      
      sampledSlideObjs  = cbind(data.frame(slide=rep(vs,times=numVirtualOBJs)), virtual_objects[sample(c(1:nrow(virtual_objects)), size=numVirtualOBJs * nvs, replace =F),])
      sampledSlideObjs[,"USE_TO_NORMALIZE"] = TRUE
      #dat0 = merge(dat0,sampledSlideObjs[,c("slide","position","Cell.ID","USE_TO_NORMALIZE")],by=c("slide","position","Cell.ID"),all.x=TRUE)
      dat0 = join(dat0,sampledSlideObjs[,c("slide","position","Cell.ID","USE_TO_NORMALIZE")],by=c("slide","position","Cell.ID"))
      
      viNorm = which(!is.na(dat0[,"USE_TO_NORMALIZE"]) & dat0[,"USE_TO_NORMALIZE"] == TRUE)
      
    }
    
  }
  
  numVirtualSlides = length(unique(elist[["PARM_EVAL_DAPI_RNDS"]]))
  if ((length(viNorm) < numVirtualSlides) || (length(unique(dat0[viNorm,"slide"])) < numVirtualSlides)) {
    cat("xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx","\n")
    cat("RUN ABORTED: viNorm does not include all virtual slides","\t",length(viNorm),"\t",numVirtualSlides,"\n")
    cat("xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx","\n")
    next
  }
  
  # Across all virtual slides, what is the min, max, mean, and median number of objects used for normalization. 
  dt2 = aggregate(dat0[viNorm,"slide"],by=list(dat0[viNorm ,"slide"]),FUN=length)
  elist[["MinNumPerSlide_NormSegObj"]] = min(dt2[,"x"])
  elist[["MeanNumPerSlide_NormSegObj"]] = mean(dt2[,"x"])
  elist[["MedianNumPerSlide_NormSegObj"]] = median(dt2[,"x"])
  elist[["MaxNumPerSlide_NormSegObj"]] = max(dt2[,"x"])

  cat("Performing Normalization...","\n")
  BMdat[,"Mean.AF000.dapi.dapi_NormInLogSpace"] = log2(1 + BMdat[,"Mean.AF000.dapi.dapi_NormInLogSpace"])
  BMdat.norm = BMdat
  BMdat.norm[,bm] = NA
  BMdat.norm[viNorm,] = get_datNorm(dat0[viNorm,], BMdat[viNorm,], method=elist[["PARM_NORM_METHOD"]], SCALE_NORM_TO_RAW_MEDIAN=elist[["PARM_SCALE_NORM_TO_RAW_MEDIAN"]], UNEVEN_MATRIX_LENGTH_METHOD = elist[["PARM_UNEVEN_MATRIX_LENGTH_METHOD"]])

  
  cat("Evaluating PARM_EVAL_SEG_OBJ Normalization...","\n")
  if (elist[["PARM_EVAL_SEG_OBJ"]]) {
    

    SlideNormFuncList_AbsSpace = getNormalizationFunctions(dat0, BMdat, BMdat.norm, bm[1], DISPLAY_FITS = displayFits)
    SlideNormFuncList_LogSpace = getNormalizationFunctions(dat0, BMdat, BMdat.norm, bm[2], DISPLAY_FITS = displayFits)

    EvalSegObj_BMdat[,"Mean.AF000.dapi.dapi_NormInLogSpace"] = log2(1 + EvalSegObj_BMdat[,"Mean.AF000.dapi.dapi_NormInLogSpace"])
    EvalSegObj_BMdat.norm = EvalSegObj_BMdat
    EvalSegObj_BMdat.norm[,bm] = NA

    EvalSegObj_BMdat.norm = applyNormalizationFunctions(EvalSegObj_dat0, EvalSegObj_BMdat, EvalSegObj_BMdat.norm, bm[1], SlideNormFuncList_AbsSpace)
    EvalSegObj_BMdat.norm = applyNormalizationFunctions(EvalSegObj_dat0, EvalSegObj_BMdat, EvalSegObj_BMdat.norm, bm[2], SlideNormFuncList_LogSpace)

    EvalSegObj_BMdat[,"Mean.AF000.dapi.dapi_NormInLogSpace"] = (2^EvalSegObj_BMdat[,"Mean.AF000.dapi.dapi_NormInLogSpace"] - 1)
    EvalSegObj_BMdat.norm[,"Mean.AF000.dapi.dapi_NormInLogSpace"] = (2^EvalSegObj_BMdat.norm[,"Mean.AF000.dapi.dapi_NormInLogSpace"] - 1)

    EvalSegObj_dbr = cbind(EvalSegObj_dat0[EvalSegObj_vi,c("slide","position", "Cell.ID")],EvalSegObj_BMdat[EvalSegObj_vi,])
    EvalSegObj_dbn = cbind(EvalSegObj_dat0[EvalSegObj_vi,c("slide","position", "Cell.ID")],EvalSegObj_BMdat.norm[EvalSegObj_vi,])

    cat("Performing.... computeGridRndToRndVariation on Evaluation objects","\n")
    EvalSegObj_dbr_RndToRndVar = computeGridRndToRndVariation(EvalSegObj_dbr)
    EvalSegObj_dbn_RndToRndVar = computeGridRndToRndVariation(EvalSegObj_dbn)
    # cat(system.time(EvalSegObj_dbr_RndToRndVar <- computeGridRndToRndVariation(EvalSegObj_dbr)),"\n")
    # cat(system.time(EvalSegObj_dbn_RndToRndVar <- computeGridRndToRndVariation(EvalSegObj_dbn)),"\n")
    cat("Completed computeGridRndToRndVariation on Evaluation objects","\n")
    
    for(scaleIndex in scaleIndexes) {
      
      if (scaleIndex == 1) {
        normScale = "Abs"
      } else if (scaleIndex == 2) {
        normScale = "Log"
      }

      
      EvalSegObj_dbrh = EvalSegObj_dbr_RndToRndVar[EvalSegObj_dbr_RndToRndVar[,paste0("n.",bm[scaleIndex])] >= elist[["PARM_EVAL_MIN_NUM_DAPI_RNDS"]] ,paste0("cv.",bm[scaleIndex])]
      EvalSegObj_dbnh = EvalSegObj_dbn_RndToRndVar[EvalSegObj_dbn_RndToRndVar[,paste0("n.",bm[scaleIndex])] >= elist[["PARM_EVAL_MIN_NUM_DAPI_RNDS"]] ,paste0("cv.",bm[scaleIndex])]

      elist[[paste0("EvalSegObj_SegObjErrCV_","Raw", normScale,"Space")]] = mean(EvalSegObj_dbrh)
      elist[[paste0("EvalSegObj_SegObjErrCV_","Norm", normScale,"Space")]] = mean(EvalSegObj_dbnh)

      # The number of objects in which they are all present (not filtered) 
      # on each of the virtual slides.
      elist[["Number_EvalSegObj"]] = length(EvalSegObj_dbrh)


      
    }
  }
  
  cat("Evaluating PARM_SEG_OBJ_NAME Normalization...","\n")
  
  BMdat.norm = extrapolateNormalization(dat0, BMdat, BMdat.norm, DISPLAY_FITS = displayFits)
  
  BMdat[,"Mean.AF000.dapi.dapi_NormInLogSpace"] = (2^BMdat[,"Mean.AF000.dapi.dapi_NormInLogSpace"] - 1)
  BMdat.norm[,"Mean.AF000.dapi.dapi_NormInLogSpace"] = (2^BMdat.norm[,"Mean.AF000.dapi.dapi_NormInLogSpace"] - 1)
  
  ve = which(!is.na(dat0[,"Mean.AF000.dapi.dapi"]) & as.numeric(dat0[,"Mean.AF000.dapi.dapi"]) >=  elist[["PARM_FILTER_MIN_DAPI_INTENSITY"]]  
             & !is.na(dat0[,"dapi_RndToRnd_Corr"])  & as.numeric(dat0[,"dapi_RndToRnd_Corr"]) >= elist[["PARM_FILTER_MIN_DAPI_RNDtoRND_CORR"]] 
             & !is.na(dat0[,"dapi_RndToRnd_Corr"])  & as.numeric(dat0[,"dapi_RndToRnd_Corr"]) >= elist[["PARM_NORM_VE_MIN_CORR"]] 
             & !is.na(dat0[,"Area"])  & as.numeric(dat0[,"Area"]) >= elist[["PARM_FILTER_MIN_OBJ_AREA"]]  & as.numeric(dat0[,"Area"]) <= elist[["PARM_FILTER_MAX_OBJ_AREA"]])
  
  if (length(ve) == 0) {
    cat("xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx","\n")
    cat("RUN ABORTED: Number of evaluation rows is zero (ve parameter is empty).","\n")
    cat("dat0 number of rows:","\t",nrow(dat0),"\n")
    cat("dat0 number of nonNA rows:","\t",nrow(dat0[!is.na(dat0[,"Mean.AF000.dapi.dapi"]), ]),"\n")
    cat("PARM_FILTER_MIN_DAPI_RNDtoRND_CORR:","\t",elist[["PARM_FILTER_MIN_DAPI_RNDtoRND_CORR"]],"\n")  
    cat("PARM_NORM_VE_MIN_CORR:","\t",elist[["PARM_NORM_VE_MIN_CORR"]],"\n")
    cat("PARM_FILTER_MIN_DAPI_INTENSITY:","\t",elist[["PARM_FILTER_MIN_DAPI_INTENSITY"]],"\n")
    cat("PARM_FILTER_MIN_OBJ_AREA:","\t",elist[["PARM_FILTER_MIN_OBJ_AREA"]],"\n")
    cat("PARM_FILTER_MAX_OBJ_AREA:","\t",elist[["PARM_FILTER_MAX_OBJ_AREA"]],"\n")
    
    cat("xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx","\n")
    next
  }
  
  
  dbr = cbind(dat0[ve,c("slide","position", "Cell.ID")],BMdat[ve,])
  dbn = cbind(dat0[ve,c("slide","position", "Cell.ID")],BMdat.norm[ve,])

  dapiRndToRndsSlides = unique(dbr[,"slide"])
  dapiRndToRndsSlides = dapiRndToRndsSlides[order(dapiRndToRndsSlides)]
  
  cat("Performing.... computeGridRndToRndVariation on Normalization objects","\n")
  dbr_RndToRndVar = computeGridRndToRndVariation(dbr)
  dbn_RndToRndVar = computeGridRndToRndVariation(dbn)
  cat("Completed computeGridRndToRndVariation on Normalization objects","\n")

  for(scaleIndex in scaleIndexes) {
    
    #scaleIndex = 1
    binWidth = 0.02
    minHistRange = 0
    maxHistRange = 1
    p1col = rgb(0.2,0.2,0.2,1/4)
    p2col = rgb(0,0,1,1/4)
    if (scaleIndex == 1) {
      normScale = "Abs"
    } else if (scaleIndex == 2) {
      normScale = "Log"
    }

    dbrh = dbr_RndToRndVar[dbr_RndToRndVar[,paste0("n.",bm[scaleIndex])] >= elist[["PARM_EVAL_MIN_NUM_DAPI_RNDS"]] ,paste0("cv.",bm[scaleIndex])]
    dbnh = dbn_RndToRndVar[dbn_RndToRndVar[,paste0("n.",bm[scaleIndex])] >= elist[["PARM_EVAL_MIN_NUM_DAPI_RNDS"]] ,paste0("cv.",bm[scaleIndex])]

    elist[[paste0("SegObjErrCV_","Raw", normScale,"Space")]] = mean(dbrh)
    elist[[paste0("SegObjErrCV_","Norm", normScale,"Space")]] = mean(dbnh)
 

    # The number of objects in which they are all present (not filtered) on each of the virtual slides.
    # This number will be equal to OR LESS THAN the minimum number of objects used for 
    # normalization on the virtual slides ("MinNumPerSlide_NormSegObj").
    elist[["Number_NormSegObj"]] = length(dbrh)


  }
  
  cat("Save Evaluation Results...","\n")
  
  if (elist[["SAVE_EVAL_DATA"]]) {
    save(file=file.path(pathToResults, elist[["SAVE_EVAL_DATA_FILENAME"]]),elist)
  }
  
  
  results = printResults(elistNum, elist, results)
  
  cns1 = colnames(eRuns)
  cns2 = colnames(results)
  resultsOutput = cbind(results,eRuns[c(1:nrow(results)),cns1[!cns1 %in% cns2]])
  write.table(resultsOutput, file=fnResults, row.names = FALSE, quote=FALSE, na = "NA", sep = "\t")
  
  
}

cat("Completed Evaluation Runs!","\n")


