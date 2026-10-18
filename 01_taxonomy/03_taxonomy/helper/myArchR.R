addIterativeLSI_with_seurat <- function(
    ArchRProj = NULL, 
    useMatrix = "TileMatrix",
    name = "IterativeLSI",
    iterations = 2,
    clusterParams = list(
      resolution = c(2), 
      sampleCells = 10000,
      maxClusters = 6,
      n.start = 10
    ),
    firstSelection = "top",
    depthCol = "nFrags",
    varFeatures = 25000,
    dimsToUse = 1:30,
    LSIMethod = 2,
    scaleDims = TRUE,
    corCutOff = 0.75,
    binarize = TRUE,
    outlierQuantiles = c(0.02, 0.98),
    filterBias = TRUE,
    sampleCellsPre = 10000,
    projectCellsPre = FALSE,
    sampleCellsFinal = NULL,
    selectionMethod = "var",
    scaleTo = 10000,
    totalFeatures = 500000,
    filterQuantile = 0.995,
    excludeChr = c(),
    saveIterations = TRUE,
    UMAPParams = list(
      n_neighbors = 40, 
      min_dist = 0.4, 
      metric = "cosine", 
      verbose = FALSE, 
      fast_sgd = TRUE
    ),
    nPlot = 10000,
    outDir = getOutputDirectory(ArchRProj),
    threads = getArchRThreads(),
    seed = 1,
    verbose = TRUE,
    force = FALSE,
    inputPeak=NULL,
    precluster=NULL,
    logFile = createLogFile("addIterativeLSI")
){
  
  
  if(verbose) message("Checking Inputs...")
  ArchR:::.validInput(input = ArchRProj, name = "ArchRProj", valid = c("ArchRProj"))
  ArchR:::.validInput(input = useMatrix, name = "useMatrix", valid = c("character"))
  ArchR:::.validInput(input = name, name = "name", valid = c("character"))
  ArchR:::.validInput(input = iterations, name = "iterations", valid = c("integer"))
  ArchR:::.validInput(input = clusterParams, name = "clusterParams", valid = c("list"))
  ArchR:::.validInput(input = varFeatures, name = "varFeatures", valid = c("integer"))
  ArchR:::.validInput(input = dimsToUse, name = "dimsToUse", valid = c("integer"))
  ArchR:::.validInput(input = LSIMethod, name = "LSIMethod", valid = c("integer", "character"))
  ArchR:::.validInput(input = scaleDims, name = "scaleDims", valid = c("boolean"))
  ArchR:::.validInput(input = corCutOff, name = "corCutOff", valid = c("numeric"))
  ArchR:::.validInput(input = binarize, name = "binarize", valid = c("boolean"))
  ArchR:::.validInput(input = outlierQuantiles, name = "outlierQuantiles", valid = c("numeric", "null"))
  ArchR:::.validInput(input = filterBias, name = "filterBias", valid = c("boolean"))
  ArchR:::.validInput(input = sampleCellsPre, name = "sampleCellsPre", valid = c("integer", "null"))
  ArchR:::.validInput(input = sampleCellsFinal, name = "sampleCellsFinal", valid = c("integer", "null"))
  ArchR:::.validInput(input = selectionMethod, name = "selectionMethod", valid = c("character"))
  ArchR:::.validInput(input = scaleTo, name = "scaleTo", valid = c("numeric"))
  ArchR:::.validInput(input = totalFeatures, name = "totalFeatures", valid = c("integer"))
  ArchR:::.validInput(input = filterQuantile, name = "filterQuantile", valid = c("numeric"))
  ArchR:::.validInput(input = excludeChr, name = "excludeChr", valid = c("character", "null"))
  ArchR:::.validInput(input = saveIterations, name = "saveIterations", valid = c("boolean"))
  ArchR:::.validInput(input = UMAPParams, name = "UMAPParams", valid = c("list"))
  ArchR:::.validInput(input = nPlot, name = "nPlot", valid = c("integer"))
  ArchR:::.validInput(input = outDir, name = "outDir", valid = c("character"))
  ArchR:::.validInput(input = threads, name = "threads", valid = c("integer"))
  ArchR:::.validInput(input = seed, name = "seed", valid = c("integer"))
  ArchR:::.validInput(input = verbose, name = "verbose", valid = c("boolean"))
  ArchR:::.validInput(input = force, name = "force", valid = c("boolean"))
  ArchR:::.validInput(input = logFile, name = "logFile", valid = c("character"))
  
  
  if(varFeatures < 1000){
    stop("Please provide more than 1000 varFeatures!")
  }
  
  ArchR:::.startLogging(logFile = logFile)
  ArchR:::.logThis(mget(names(formals()),sys.frame(sys.nframe())), "IterativeLSI Input-Parameters", logFile=logFile)
  
  ArchR:::.requirePackage("Matrix", source = "cran")
  tstart <- Sys.time()
  
  if(!is.null(ArchRProj@reducedDims[[name]])){
    if(!force){
      stop("Error name in reducedDims Already Exists! Set force = TRUE or pick a different name!")
    }
  }
  
  #Set Seed
  set.seed(seed)
  outDir <- file.path(outDir, name)
  dir.create(outDir, showWarnings = FALSE, recursive = TRUE)
  
  #All the Cell Names
  cellNames <- rownames(getCellColData(ArchRProj))
  if(!is.null(sampleCellsPre)){
    if(length(cellNames) < sampleCellsPre){
      sampleCellsPre <- NULL
    }
  }
  if(!is.null(sampleCellsFinal)){
    if(length(cellNames) <= sampleCellsFinal){
      sampleCellsFinal <- NULL
    }
  }
  
  #MatrixFiles
  ArrowFiles <- getSampleColData(ArchRProj)[,"ArrowFiles"]
  
  #Check if Matrix is supported and check type
  if(tolower(useMatrix) == "tilematrix"){
    useMatrix <- "TileMatrix"
    tileSizes <- lapply(ArrowFiles, function(x){
      h5read(x, "TileMatrix/Info/Params/")$tileSize[1]
    }) %>% unlist
    if(length(unique(tileSizes)) != 1){
      stop("Error not all TileMatrices are the same tileSize!")
    }
    tileSize <- unique(tileSizes)
  }else if(tolower(useMatrix) == "peakmatrix"){
    useMatrix <- "PeakMatrix"
    tileSize <- NA
  }else{
    tileSize <- NA
  }
  
  units <- unique(unlist(lapply(ArrowFiles, function(x) h5read(x, paste0(useMatrix, "/Info/Units")))))
  if(length(units) != 1){
    stop("Units of matrices are not identical!")
  }
  if(grepl("log",units,ignore.case=TRUE)){
    stop("Cannot use log transformed values for iterativeLSI!")
  }
  
  tstart <- Sys.time()  
  ############################################################################################################################
  # Organize Information for LSI
  ############################################################################################################################
  chrToRun <- ArchR:::.availableSeqnames(ArrowFiles, subGroup = useMatrix)
  
  if(tolower(firstSelection) == "top"){
    
  #  if(!binarize){
  #    matClass <- h5read(ArrowFiles[1], paste0(useMatrix,"/Info/Class"))
  #    if(matClass != "Sparse.Binary.Matrix"){
  #      stop("Input matrix is not binarized and binarize != TRUE. Please use binarized data if using top selection for first iteration! Set binarize = TRUE!")
  #    }
  #  }
    
    #Compute Row Sums Across All Samples
    ArchR:::.logDiffTime("Computing Total Across All Features", tstart, addHeader = FALSE, verbose = verbose, logFile = logFile)
    if(useMatrix == "TileMatrix"){
      totalAcc <- ArchR:::.getRowSums(ArrowFiles = ArrowFiles, useMatrix = useMatrix, seqnames = chrToRun, addInfo = FALSE)
      totalAcc$start <- (totalAcc$idx - 1) * tileSize
    }else{
      totalAcc <- ArchR:::.getRowSums(ArrowFiles = ArrowFiles, useMatrix = useMatrix, seqnames = chrToRun, addInfo = TRUE)
    }
    
    #Filter Chromosomes
    if(length(excludeChr) > 0){
      totalAcc <- totalAcc[BiocGenerics::which(totalAcc$seqnames %bcni% excludeChr), , drop = FALSE]
    }
    if(!is.null(inputPeak)){
      totalAcc=totalAcc[overlapsAny(makeGRangesFromDataFrame(totalAcc),inputPeak),]
    }
    #Identify the top features to be used here
    ArchR:::.logDiffTime("Computing Top Features", tstart, addHeader = FALSE, verbose = verbose, logFile = logFile)
    nFeature <- varFeatures[1]
    rmTop <- floor((1-filterQuantile) * totalFeatures)
    topIdx <- head(order(totalAcc$rowSums, decreasing=TRUE), nFeature + rmTop)[-seq_len(rmTop)]
    topFeatures <- totalAcc[sort(topIdx),]
    
    gc()
    
  }else if(tolower(firstSelection) %in% c("var", "variable")){
    
    if(binarize){
      stop("Please do not binarize data if using variable selection for first iteration! Set binarize = FALSE!")
    }
    
    if(units %in% "BinarizedCounts"){
      stop("Cannot do variable selection with BinarizedCounts. Set firstSelection = Top!")
    }
    
    #Compute Row Sums Across All Samples
    ArchR:::.logDiffTime("Computing Variability Across All Features", tstart, addHeader = FALSE, verbose = verbose, logFile = logFile)
    if(useMatrix == "TileMatrix"){
      totalAcc <- ArchR:::.getRowVars(ArrowFiles = ArrowFiles, useMatrix = useMatrix, seqnames = chrToRun, useLog2 = TRUE)
      totalAcc$start <- (totalAcc$idx - 1) * tileSize
    }else{
      totalAcc <- ArchR:::.getRowVars(ArrowFiles = ArrowFiles, useMatrix = useMatrix, seqnames = chrToRun, useLog2 = TRUE)
    }
    
    #Filter Chromosomes
    if(length(excludeChr) > 0){
      totalAcc <- totalAcc[BiocGenerics::which(totalAcc$seqnames %bcni% excludeChr), , drop = FALSE]
    }
    if(!is.null(inputPeak)){
      totalAcc=totalAcc[overlapsAny(makeGRangesFromDataFrame(totalAcc),inputPeak),]
    }
    #Identify the top features to be used here
    ArchR:::.logDiffTime("Computing Variable Features", tstart, addHeader = FALSE, verbose = verbose, logFile = logFile)
    nFeature <- varFeatures[1]
    #if(nFeature > 0.5 * nrow(totalAcc)){
    #  stop("nFeature for variable selection must be less than 1/2 the total features!")
    #}
    topIdx <- head(order(totalAcc$combinedVars, decreasing=TRUE), nFeature)
    topFeatures <- totalAcc[sort(topIdx),]
    gc()
    
  }else{
    
    stop("firstSelect method must be Top or Var/Variable!")
    
  }
  totalAcc$IDX=seq_len(nrow(totalAcc))
  cellDepth <- tryCatch({
    df <- getCellColData(ArchRProj = ArchRProj, select = depthCol)
    v <- df[,1]
    names(v) <- rownames(df)
    v
  }, error = function(e){
    tryCatch({
      ArchR:::.getColSums(ArrowFiles = ArrowFiles, useMatrix = useMatrix, seqnames = chrToRun)[ArchRProj$cellNames]
    }, error = function(y){
      stop("Could not determine depth from depthCol or colSums!")
    })
  }
  )
  cellDepth <- log10(cellDepth + 1)
  ############################################################################################################################
  # LSI Iteration 1
  ############################################################################################################################
  if(!is.null(precluster)){
    variableFeatures <- identifyVarFeatures_from_Seurat(
      seuratobj = precluster,
      ArrowFiles = ArrowFiles,
      useMatrix = useMatrix,
      prevFeatures = totalAcc,
      scaleTo = scaleTo,
      totalAcc = totalAcc,
      totalFeatures = totalFeatures,
      firstSelection = firstSelection,
      selectionMethod = selectionMethod,
      varFeatures = varFeatures,
      tstart = tstart,
      threads = threads,
      verbose = verbose,
      logFile = logFile
    )
  }else{
    ArchR:::.logDiffTime(paste0("Running LSI (1 of ",iterations,") on Top Features"), tstart, addHeader = TRUE, verbose = verbose, logFile = logFile)
    j <- 1
    
    if(!is.null(clusterParams$sampleCells)){
      if(!is.na(clusterParams$sampleCells[j])){
        sampleJ <- clusterParams$sampleCells[j]
      }else if(!is.na(clusterParams$sampleCells[1])){
        sampleJ <- clusterParams$sampleCells[1]
      }else{
        sampleJ <- sampleCellsPre
      }
    }else{
      sampleJ <- sampleCellsPre 
    }
    ##########################################################
    #############################seurat 
    mat <- PartialMatrix(
      ArrowFiles = ArrowFiles, 
      featureDF = totalAcc,
      cellNames = cellNames, 
      cellDepth = cellDepth,
      useMatrix = useMatrix,
      sampleNames = getCellColData(ArchRProj)$Sample, 
      LSIMethod = LSIMethod, 
      scaleTo = scaleTo,
      dimsToUse = dimsToUse, 
      binarize = binarize, 
      outlierQuantiles = outlierQuantiles,
      sampleCells = sampleCellsPre,
      projectAll = F,
      threads = threads,
      useIndex = FALSE,
      seed = seed,
      tstart = tstart,
      verbose = verbose,
      logFile = logFile
    )
    rowc=rownames(mat)=paste0(totalAcc$seqnames,":",totalAcc$start,'-',totalAcc$end)
    colc=colnames(mat)
    library(Signac)
    library(Seurat)
    library(harmony)
    obj1=CreateSeuratObject(mat, project='scATAC', min.cells=0, min.features=0,meta.data = as.data.frame(ArchRProj@cellColData[colc,]))                      
    
    obj1 <- RunTFIDF(obj1)
    if(F){ ################# 
      npeaks <- colSums(obj1@assays$RNA@data)
      obj1[['RNA']]@data <- tcrossprod(x = obj1[['RNA']]@data, y = Diagonal(x = 1 / npeaks))
      obj1[['RNA']]@data@x <- log(obj1[['RNA']]@data@x* 1e4 + 1)
      colnames(obj1[['RNA']]@data)=colnames(mat)  
    }
    
    obj1 <- FindTopFeatures(obj1, min.cutoff = 'q0')
    obj1 <- RunSVD(obj1)
    obj1 <- RunHarmony(obj1, group.by.vars = 'Sample', 
                       project.dim = F,reduction='lsi',
                       reduction.save = "harmony")
    obj1 <- FindNeighbors(object =obj1, reduction = "harmony", dims = 2:30)
    obj1 <- FindClusters(object = obj1, verbose = FALSE, algorithm = 3,resolution = clusterParams$resolution[1])
    obj1 <- RunUMAP(object = obj1, reduction = "harmony", dims = 2:30)
    pdf('seurat.pdf')
    p1=DimPlot(object = obj1, label = TRUE) + NoLegend()
    p2=DimPlot(object = obj1, group.by = 'newID') + NoLegend()
    if('GEX_anno0_1'%in%colnames(obj1@meta.data)){
      p3=DimPlot(object = obj1, group.by = 'GEX_anno0_1',label = TRUE) + NoLegend()  
      print(p3)
    }
    
    p4=DimPlot(object = obj1, group.by = 'seurat_clusters',label = TRUE) + NoLegend()
    print(p1)
    print(p2)
    
    print(p4)
    dev.off()
    # save(obj1,file='seurat.Rdata')
    variableFeatures <- identifyVarFeatures_from_Seurat(
      seuratobj = obj1,
      ArrowFiles = ArrowFiles,
      useMatrix = useMatrix,
      prevFeatures = totalAcc,
      scaleTo = scaleTo,
      totalAcc = totalAcc,
      totalFeatures = totalFeatures,
      firstSelection = firstSelection,
      selectionMethod = selectionMethod,
      varFeatures = varFeatures,
      tstart = tstart,
      threads = threads,
      verbose = verbose,
      logFile = logFile
    )
    rm(obj1)
  }


  gc()
  #########################################
  
  outLSI <- LSIPartialMatrix(
    ArrowFiles = ArrowFiles, 
    featureDF = variableFeatures,
    useMatrix = useMatrix,
    cellNames = cellNames, 
    cellDepth = cellDepth,
    sampleNames = getCellColData(ArchRProj)$Sample, 
    LSIMethod = LSIMethod, 
    scaleTo = scaleTo, 
    dimsToUse = dimsToUse,
    binarize = binarize,
    outlierQuantiles = outlierQuantiles, 
    sampleCells = sampleCellsFinal,
    projectAll = T,
    threads = threads,
    useIndex = FALSE,
    seed = seed,
    tstart = tstart,
    verbose = verbose,
    logFile = logFile
  )
  outLSI$scaleDims <- scaleDims
  outLSI$useMatrix <- useMatrix
  outLSI$tileSize <- tileSize
  gc()
  #Organize Output
  ArchR:::.logDiffTime("Finished Running IterativeLSI", tstart, addHeader = FALSE, verbose = verbose, logFile = logFile)
  ArchRProj@reducedDims[[name]] <- outLSI
  
  return(ArchRProj)
  
}


LSIPartialMatrix <- function(
    ArrowFiles = NULL, 
    featureDF = NULL, 
    useMatrix = NULL,
    cellNames = NULL, 
    cellDepth = NULL,
    sampleNames = NULL, 
    dimsToUse = NULL, 
    binarize = TRUE, 
    outlierQuantiles = c(0.02, 0.98),
    LSIMethod = FALSE,
    scaleTo = 10^4,
    sampleCells = 5000, 
    projectAll = TRUE,
    threads = 1, 
    seed = 1,
    useIndex = FALSE, 
    tstart = NULL, 
    verbose = TRUE,
    logFile = NULL
){
  
  if(is.null(tstart)){
    tstart <- Sys.time()
  }
  
  errorList <- append(args, mget(names(formals()),sys.frame(sys.nframe())))
  
  outLSI2 <- tryCatch({
    
    if(is.null(sampleCells)){
      
      #Construct Matrix
      ArchR:::.logDiffTime("Creating Partial Matrix", tstart, addHeader = FALSE, verbose = verbose, logFile = logFile)
      
      mat <- ArchR:::.getPartialMatrix(
        ArrowFiles = ArrowFiles,
        featureDF = featureDF,
        useMatrix = useMatrix,
        cellNames = cellNames,
        doSampleCells = FALSE,
        threads = threads,
        verbose = FALSE,
        logFile = logFile
      )
      
      #Compute LSI
      ArchR:::.logDiffTime("Computing LSI", tstart, addHeader = FALSE, verbose = verbose, logFile = logFile)
      outLSI <- computeLSI(
        mat = mat, 
        LSIMethod = LSIMethod, 
        scaleTo = scaleTo,
        nDimensions = max(dimsToUse),
        binarize = binarize, 
        outlierQuantiles = outlierQuantiles,
        verbose = FALSE, 
        seed = seed,
        tstart = tstart,
        logFile = logFile
      )
      outLSI$LSIFeatures <- featureDF
      outLSI$corToDepth <- list(
        scaled = abs(cor(ArchR:::.scaleDims(outLSI[[1]]), cellDepth[rownames(outLSI[[1]])]))[,1],
        none = abs(cor(outLSI[[1]], cellDepth[rownames(outLSI[[1]])]))[,1]
      )
      
    }else{
      
      sampledCellNames <- sampleBySample(
        cellNames = cellNames,
        sampleNames = sampleNames,
        cellDepth = cellDepth,
        sampleCells = sampleCells,
        outlierQuantiles = outlierQuantiles,
        factor = 2      
      )
      ArchR:::.logDiffTime(sprintf("Sampling Cells (N = %s) for Estimated LSI", length(sampledCellNames)), tstart, addHeader = FALSE, verbose = verbose, logFile = logFile)
      
      #Construct Sampled Matrix
      ArchR:::.logDiffTime("Creating Sampled Partial Matrix", tstart, addHeader = FALSE, verbose = verbose, logFile = logFile)
      o <- h5closeAll()
      if(!projectAll){
        
        mat <- ArchR:::.getPartialMatrix(
          ArrowFiles = ArrowFiles,
          featureDF = featureDF,
          useMatrix = useMatrix,
          cellNames = sampledCellNames,
          doSampleCells = FALSE,
          threads = threads,
          verbose = FALSE,
          logFile = logFile
        )
        
        #Compute LSI
        ArchR:::.logDiffTime("Computing Estimated LSI (projectAll = FALSE)", tstart, addHeader = FALSE, verbose = verbose, logFile = logFile)
        outLSI <- computeLSI(
          mat = mat, 
          LSIMethod = LSIMethod, 
          scaleTo = scaleTo,
          nDimensions = max(dimsToUse),
          binarize = binarize, 
          outlierQuantiles = outlierQuantiles,
          seed = seed,
          tstart = tstart,
          logFile = logFile
        )
        outLSI$LSIFeatures <- featureDF
        outLSI$corToDepth <- list(
          scaled = abs(cor(ArchR:::.scaleDims(outLSI[[1]]), cellDepth[rownames(outLSI[[1]])]))[,1],
          none = abs(cor(outLSI[[1]], cellDepth[rownames(outLSI[[1]])]))[,1]
        )
        
      }else{
        
        tmpPath <- ArchR:::.tempfile(pattern = "tmp-LSI-PM")
        
        ArchR:::.logDiffTime(sprintf("Sampling Cells (N = %s) for Estimated LSI", length(sampledCellNames)), tstart, addHeader = FALSE, verbose = verbose, logFile = logFile)
        out <- ArchR:::.getPartialMatrix(
          ArrowFiles = ArrowFiles,
          featureDF = featureDF,
          useMatrix = useMatrix,
          cellNames = cellNames,
          doSampleCells = TRUE,
          sampledCellNames = sampledCellNames,
          tmpPath = tmpPath,
          useIndex = useIndex,
          threads = threads,
          verbose = FALSE,
          logFile = logFile
        )
        gc()
        
        #Perform LSI on Partial Sampled Matrix
        ArchR:::.logDiffTime("Computing Estimated LSI (projectAll = TRUE)", tstart, addHeader = FALSE, verbose = verbose, logFile = logFile)
        outLSI <- computeLSI(
          mat = out$mat, 
          LSIMethod = LSIMethod, 
          scaleTo = scaleTo,
          nDimensions = max(dimsToUse),
          binarize = binarize, 
          outlierQuantiles = outlierQuantiles,
          seed = seed,
          tstart = tstart,
          logFile = logFile
        )
        
        tmpMatFiles <- out[[2]]
        rm(out)
        gc()
        
        #Read In Matrices and Project into Manifold
        #Do Threads / 3 just in case of memory here (Needs testing QQQ)
        threads2 <- 1 #max(ceiling(threads / 3), 1)
        ArchR:::.logDiffTime("Projecting Matrices with LSI-Projection (Granja* et al 2019)", tstart, addHeader = FALSE, verbose = verbose, logFile = logFile)
        pLSI <- ArchR:::.safelapply(seq_along(tmpMatFiles), function(x){
          ArchR:::.logDiffTime(sprintf("Projecting Matrix (%s of %s) with LSI-Projection", x, length(tmpMatFiles)), tstart, addHeader = FALSE, verbose = FALSE, logFile = logFile)
          projectLSI(mat = readRDS(tmpMatFiles[x]), LSI = outLSI, verbose = FALSE, tstart = tstart, logFile = logFile)
        }, threads = threads2) %>% Reduce("rbind", .)
        
        #Remove Temporary Matrices
        rmf <- file.remove(tmpMatFiles)
        
        #Set To LSI the SVD Matrices
        outLSI$exlcude <- cellNames[which(cellNames %ni% rownames(pLSI))]
        outLSI$matSVD <- as.matrix(pLSI[cellNames[which(cellNames %in% rownames(pLSI))],])
        
      }
      
      outLSI$LSIFeatures <- featureDF
      outLSI$corToDepth <- list(
        scaled = abs(cor(ArchR:::.scaleDims(outLSI[[1]]), cellDepth[rownames(outLSI[[1]])]))[,1],
        none = abs(cor(outLSI[[1]], cellDepth[rownames(outLSI[[1]])]))[,1]
      )
      
    }
    
    outLSI
    
  }, error = function(e){
    
    errorList$outLSI <- if(exists("outLSI", inherits = FALSE)) outLSI else "Error with outLSI!"
    errorList$matSVD <- if(exists("outLSI", inherits = FALSE)) outLSI[[1]] else "Error with matSVD!"
    
    ArchR:::.logError(e, fn = ".LSIPartialMatrix", info = "", errorList = errorList, logFile = logFile)
    
  })
  
  return(outLSI2)
  
}


#########################################################################################
# Sampling Helpers
#########################################################################################

sampleBySample <- function(
    cellNames = NULL, 
    cellDepth = NULL, 
    sampleNames = NULL, 
    sampleCells = NULL, 
    outlierQuantiles = NULL, 
    factor = 2
){
  
  if(sampleCells < length(cellNames)){
    
    sampleN <- ceiling(sampleCells / length(sampleNames)* table(sampleNames) )
    splitCells <- split(cellNames, sampleNames)
    splitDepth <- split(cellDepth, sampleNames)
    
    sampledCellNames <- lapply(seq_along(splitCells), function(x){
      ArchR:::.filterSample(
        x = splitCells[[x]], 
        n = sampleN[names(splitCells)[x]], 
        vals = splitDepth[[x]], 
        outlierQuantiles = outlierQuantiles, 
        factor = factor
      )
    }) %>% unlist %>% sort
    
    return(sampledCellNames)
    
  }else{
    
    return(cellNames)
    
  }
  
}
PartialMatrix <- function(
    ArrowFiles = NULL, 
    featureDF = NULL, 
    useMatrix = NULL,
    cellNames = NULL, 
    cellDepth = NULL,
    sampleNames = NULL, 
    dimsToUse = NULL, 
    binarize = TRUE, 
    outlierQuantiles = c(0.02, 0.98),
    LSIMethod = FALSE,
    scaleTo = 10^4,
    sampleCells = 5000, 
    projectAll = TRUE,
    threads = 1, 
    seed = 1,
    useIndex = FALSE, 
    tstart = NULL, 
    verbose = TRUE,
    logFile = NULL
){
  if(is.null(sampleCells)){
    
    #Construct Matrix
    ArchR:::.logDiffTime("Creating Partial Matrix", tstart, addHeader = FALSE, verbose = verbose, logFile = logFile)
    
    mat <- ArchR:::.getPartialMatrix(
      ArrowFiles = ArrowFiles,
      featureDF = featureDF,
      useMatrix = useMatrix,
      cellNames = cellNames,
      doSampleCells = FALSE,
      threads = threads,
      verbose = FALSE,
      logFile = logFile
    )
  }else{
    
    sampledCellNames <- sampleBySample(
      cellNames = cellNames,
      sampleNames = sampleNames,
      cellDepth = cellDepth,
      sampleCells = sampleCells,
      outlierQuantiles = outlierQuantiles,
      factor = 2      
    )
    ArchR:::.logDiffTime(sprintf("Sampling Cells (N = %s) for Estimated LSI", length(sampledCellNames)), tstart, addHeader = FALSE, verbose = verbose, logFile = logFile)
    
    #Construct Sampled Matrix
    ArchR:::.logDiffTime("Creating Sampled Partial Matrix", tstart, addHeader = FALSE, verbose = verbose, logFile = logFile)
    o <- h5closeAll()
    mat <- ArchR:::.getPartialMatrix(
      ArrowFiles = ArrowFiles,
      featureDF = featureDF,
      useMatrix = useMatrix,
      cellNames = sampledCellNames,
      doSampleCells = FALSE,
      threads = threads,
      verbose = FALSE,
      logFile = logFile)
    
    
  }
  return(mat)
}



identifyVarFeatures_from_Seurat <- function(
    seuratobj = NULL,
    clusterDF = NULL,
    prevFeatures = NULL,
    ArrowFiles = NULL,
    useMatrix = NULL,
    totalAcc = NULL,
    scaleTo = NULL,
    firstSelection = NULL,
    totalFeatures = NULL,
    selectionMethod = NULL,
    varFeatures = NULL,
    tstart = NULL,
    threads = NULL,
    verbose = NULL,
    logFile = NULL
){
  
  errorList <- append(args, mget(names(formals()),sys.frame(sys.nframe())))
  
  variableFeatures2 <- tryCatch({
    
    nClust <- length(unique(seuratobj$seurat_clusters))
    
    if(nClust == 1){
      
      ArchR:::.logDiffTime("Identified 1 Cluster, Continuing with Previous Features!", tstart, addHeader = FALSE, verbose = verbose, logFile = logFile)
      return(prevFeatures)
      
    }else{
      
      #Random Sampling for Quick Estimation of Variance
      
      
      #idxSub <- seq_len(nrow(outLSI$matSVD))
      ArchR:::.logDiffTime("Creating Cluster Matrix on the total Group Features", tstart, addHeader = FALSE, verbose = verbose, logFile = logFile)
      groupList <- SimpleList(split(rownames(seuratobj@meta.data), seuratobj$seurat_clusters))
      
      if(tolower(firstSelection) == "top"){
        groupFeatures <- totalAcc[sort(head(order(totalAcc$rowSums, decreasing = TRUE), totalFeatures)),]
      }else if(tolower(firstSelection) %in% c("var", "variable")){
        groupFeatures <- totalAcc
      }
      
      groupMat <- ArchR:::.getGroupMatrix(
        ArrowFiles = ArrowFiles, 
        featureDF = groupFeatures,
        useMatrix = useMatrix, 
        threads = threads,
        groupList = groupList,
        useIndex = FALSE,
        verbose = FALSE
      )
      
      ArchR:::.logDiffTime("Computing Variable Features", tstart, addHeader = FALSE, verbose = verbose, logFile = logFile)
      if(length(varFeatures) > 1){
        nFeature <- varFeatures[j]
      }else{
        nFeature <- varFeatures
      }
      
      if(tolower(selectionMethod) == "var"){
        
        #Log-Normalize
        groupMat <- log2(t(t(groupMat) / colSums(groupMat)) * scaleTo + 1)
        var <- matrixStats::rowVars(groupMat)
        idx <- sort(head(order(var,decreasing=TRUE), nFeature))
        variableFeatures <- groupFeatures[idx,]
        
      }else if(tolower(selectionMethod) == "vmr"){
        
        #Variance-to-Mean Ratio
        vmr <- matrixStats::rowVars(groupMat) / rowMeans(groupMat)
        idx <- sort(head(order(vmr, decreasing=TRUE), nFeature))
        variableFeatures <- groupFeatures[idx,]
        
      }else{
        
        stop("Error selectionMethod is not valid requires var or vmr")
        
      }
      
    }
    
    variableFeatures
    
  }, error = function(e){
    
    ArchR:::.logError(e, fn = ".identifyVarFeatures", info = "", errorList = errorList, logFile = logFile)
    
  })
  
  return(variableFeatures2)
  
}
computeLSI <- function(
    mat = NULL, 
    LSIMethod = 1,
    scaleTo = 10^4,
    nDimensions = 50, 
    binarize = TRUE, 
    outlierQuantiles = c(0.02, 0.98),
    seed = 1, 
    verbose = FALSE, 
    tstart = NULL,
    logFile = NULL
){
  
  ArchR:::.logThis(append(args, mget(names(formals()),sys.frame(sys.nframe()))), "LSI Parameters", logFile=logFile)
  
  out2 <- tryCatch({
    
    set.seed(seed)
    
    if(is.null(tstart)){
      tstart <- Sys.time()
    }
    
    ArchR:::.logDiffTime(sprintf("Running LSI, Input Matrix = %s GB", round(object.size(mat)/10^9, 3)), tstart, addHeader = verbose, verbose = verbose, logFile = logFile)
    
    #TF IDF LSI adapted from flyATAC
    if(binarize){
      ArchR:::.logDiffTime("Binarizing Matrix", tstart, addHeader = FALSE, verbose = verbose, logFile = logFile)
      mat@x[mat@x > 0] <- 1 
    }
    
    #Compute Col Sums
    ArchR:::.logDiffTime("Computing Term Frequency", tstart, addHeader = FALSE, verbose = verbose, logFile = logFile)
    colSm <- Matrix::colSums(mat)
    if(any(colSm == 0)){
      exclude <- which(colSm==0)
      mat <- mat[,-exclude, drop = FALSE]
      colSm <- colSm[-exclude]
    }else{
      exclude <- c()
    }
    
    cn <- colnames(mat)
    filterOutliers <- 0
    if(!is.null(outlierQuantiles)){
      qCS <- quantile(colSm, probs = c(min(outlierQuantiles), max(outlierQuantiles)))
      idxOutlier <- which(colSm <= qCS[1] | colSm >= qCS[2])
      if(length(idxOutlier) > 0){
        ArchR:::.logDiffTime("Filtering Outliers Based On Counts", tstart, addHeader = FALSE, verbose = verbose, logFile = logFile)
        #.safeSaveRDS(mat, "temp.rds", compress = FALSE)
        matO <- mat[, idxOutlier, drop = FALSE]
        mat <- mat[, -idxOutlier, drop = FALSE]
        mat2 <- mat[, head(seq_len(ncol(mat)), 50), drop = FALSE] # A 2nd Matrix to Check Projection is Working
        colSm <- colSm[-idxOutlier]
        filterOutliers <- 1       
      }
    }
    
    #Clean up zero rows
    ArchR:::.logDiffTime("Removing 0 Sum Rows", tstart, addHeader = FALSE, verbose = verbose, logFile = logFile)
    rowSm <- Matrix::rowSums(mat)
    idx <- which(rowSm > 0)
    mat <- mat[idx, ]
    rowSm <- rowSm[idx]
    
    #TF - Normalize
    mat@x <- mat@x / rep.int(colSm, Matrix::diff(mat@p))
    
    if(LSIMethod == 1 | tolower(LSIMethod) == "tf-logidf"){
      
      #Adapted from Casanovich et al.
      
      #LogIDF
      ArchR:::.logDiffTime("Computing Inverse Document Frequency", tstart, addHeader = FALSE, verbose = verbose, logFile = logFile)
      idf   <- as(log(1 + ncol(mat) / rowSm), "sparseVector")
      
      #TF-LogIDF
      ArchR:::.logDiffTime("Computing TF-IDF Matrix", tstart, addHeader = FALSE, verbose = verbose, logFile = logFile)
      mat <- as(Matrix::Diagonal(x=as.vector(idf)), "sparseMatrix") %*% mat
      
    }else if(LSIMethod == 2 | tolower(LSIMethod) == "log(tf-idf)"){
      
      #Adapted from Stuart et al.
      
      #IDF
      ArchR:::.logDiffTime("Computing Inverse Document Frequency", tstart, addHeader = FALSE, verbose = verbose, logFile = logFile)
      idf   <- as(ncol(mat) / rowSm, "sparseVector")
      
      #TF-IDF
      ArchR:::.logDiffTime("Computing TF-IDF Matrix", tstart, addHeader = FALSE, verbose = verbose)
      mat <- as(Matrix::Diagonal(x=as.vector(idf)), "sparseMatrix") %*% mat
      
      #Log transform TF-IDF
      mat@x <- log(mat@x * scaleTo + 1)  
      
    }else if(LSIMethod == 3 | tolower(LSIMethod) == "logtf-logidf"){
      
      #LogTF
      mat@x <- log(mat@x + 1)
      
      #LogIDF
      ArchR:::.logDiffTime("Computing Inverse Document Frequency", tstart, addHeader = FALSE, verbose = verbose, logFile = logFile)
      idf   <- as(log(1 + ncol(mat) / rowSm), "sparseVector")
      
      #TF-IDF
      ArchR:::.logDiffTime("Computing TF-IDF Matrix", tstart, addHeader = FALSE, verbose = verbose, logFile = logFile)
      mat <- as(Matrix::Diagonal(x=as.vector(idf)), "sparseMatrix") %*% mat
      
      
    }else if (LSIMethod == 4 ){
      # log TF-log IDF * norm total + log1p
      mat@x <- log(mat@x + 1)
      
      #LogIDF
      ArchR:::.logDiffTime("Computing Inverse Document Frequency", tstart, addHeader = FALSE, verbose = verbose, logFile = logFile)
      idf   <- as(log(1 + ncol(mat) / rowSm), "sparseVector")
      
      #TF-IDF
      ArchR:::.logDiffTime("Computing TF-IDF Matrix", tstart, addHeader = FALSE, verbose = verbose, logFile = logFile)
      mat <- as(Matrix::Diagonal(x=as.vector(idf)), "sparseMatrix") %*% mat
      colSm <- Matrix::colSums(mat)
      mat@x <- mat@x / rep.int(colSm, Matrix::diff(mat@p))
      mat@x <- log(mat@x* scaleTo + 1)
    }else{
      
      stop("LSIMethod unrecognized please select valid method!")
      
    }
    
    gc()
    
    #Calc SVD then LSI
    ArchR:::.logDiffTime("Computing SVD using irlba", tstart, addHeader = FALSE, verbose = verbose, logFile = logFile)
    svd <- irlba::irlba(mat, nDimensions, nDimensions)
    svdDiag <- matrix(0, nrow=nDimensions, ncol=nDimensions)
    diag(svdDiag) <- svd$d
    matSVD <- t(svdDiag %*% t(svd$v))
    rownames(matSVD) <- colnames(mat)
    colnames(matSVD) <- paste0("LSI",seq_len(ncol(matSVD)))
    
    #Return Object
    out <- SimpleList(
      matSVD = matSVD, 
      rowSm = rowSm, 
      nCol = length(colSm),
      exclude = exclude, 
      idx = idx, 
      svd = svd, 
      binarize = binarize, 
      scaleTo = scaleTo,
      nDimensions = nDimensions,
      LSIMethod = LSIMethod,
      outliers = NA,
      date = Sys.Date(),
      seed = seed
    )
    
    if(filterOutliers == 1){
      ArchR:::.logDiffTime("Projecting Outliers with LSI-Projection (Granja* et al 2019)", tstart, addHeader = FALSE, verbose = verbose, logFile = logFile)
      #Quick Check LSI-Projection Works
      pCheck <- projectLSI(mat = mat2, LSI = out, verbose = verbose, logFile = logFile)
      pCheck2 <- out[[1]][rownames(pCheck), ]
      pCheck3 <- lapply(seq_len(ncol(pCheck)), function(x){
        cor(pCheck[,x], pCheck2[,x])
      }) %>% unlist
      if(min(pCheck3) < 0.95){
        ArchR:::.logThis(pCheck, "pCheck", logFile=logFile)
        ArchR:::.logThis(pCheck2, "pCheck2", logFile=logFile)
        ArchR:::.logThis(pCheck3, "pCheck3", logFile=logFile)
        warning("Warning with LSI-projection! Cor less than 0.95 of re-projection. Please report this to github with logFile!")
      }
      #Project LSI Outliers
      out$outliers <- colnames(matO)
      outlierLSI <- projectLSI(mat = matO, LSI = out, verbose = verbose, logFile = logFile)
      allLSI <- rbind(out[[1]], outlierLSI)
      allLSI <- allLSI[cn, , drop = FALSE] #Re-Order Correctly to original
      out[[1]] <- allLSI
    }
    ArchR:::.logDiffTime("Finished LSI (TF-IDF SVD) using irlba", tstart, addHeader = FALSE, verbose = verbose, logFile = logFile)
    
    rm(mat)
    gc()
    
    out
    
  }, error = function(e){
    
    errorList <- list(
      mat = mat,
      colSm = if(exists("colSm", inherits = FALSE)) colSm else "Error with colSm!",
      idf = if(exists("idf", inherits = FALSE)) idf else "Error with idf!",
      V = if(exists("V", inherits = FALSE)) V else "Error with V!",
      matSVD = if(exists("matSVD", inherits = FALSE)) matSVD else "Error with matSVD!",
      out = if(exists("out", inherits = FALSE)) out else "Error with out!"
    )
    
    ArchR:::.logError(e, fn = ".computeLSI", info = "", errorList = errorList, logFile = logFile)
    
  })
  
  out2
  
}
projectLSI <- function(
    mat = NULL, 
    LSI = NULL, 
    returnModel = FALSE, 
    verbose = FALSE, 
    tstart = NULL,
    logFile = NULL
){   
  
  ArchR:::.logThis(append(args, mget(names(formals()),sys.frame(sys.nframe()))), "LSI-Projection Parameters", logFile=logFile)
  
  out2 <- tryCatch({
    
    require(Matrix)
    set.seed(LSI$seed)
    
    if(is.null(tstart)){
      tstart <- Sys.time()
    }
    
    ArchR:::.logDiffTime(sprintf("Projecting LSI, Input Matrix = %s GB", round(object.size(mat)/10^9, 3)), tstart, addHeader = verbose, verbose = verbose, logFile = logFile)
    
    #Get Same Features
    ArchR:::.logDiffTime("Subsetting by Non-Zero features in inital Matrix", tstart, addHeader = FALSE, verbose = verbose, logFile = logFile)
    mat <- mat[LSI$idx,]
    
    #Binarize Matrix
    if(LSI$binarize){
      ArchR:::.logDiffTime("Binarizing Matrix", tstart, addHeader = FALSE, verbose = verbose, logFile = logFile)
      mat@x[mat@x > 0] <- 1       
    }
    
    #TF
    ArchR:::.logDiffTime("Computing Term Frequency", tstart, addHeader = FALSE, verbose = verbose, logFile = logFile)
    colSm <- Matrix::colSums(mat)
    if(any(colSm == 0)){
      exclude <- which(colSm==0)
      mat <- mat[,-exclude]
      colSm <- colSm[-exclude]
    }
    mat@x <- mat@x / rep.int(colSm, Matrix::diff(mat@p))
    
    if(LSI$LSIMethod == 1 | tolower(LSI$LSIMethod) == "tf-logidf"){
      
      #Adapted from Casanovich et al.
      
      #LogIDF
      ArchR:::.logDiffTime("Computing Inverse Document Frequency", tstart, addHeader = FALSE, verbose = verbose, logFile = logFile)
      idf   <- as(log(1 + LSI$nCol / LSI$rowSm), "sparseVector")
      
      #TF-LogIDF
      ArchR:::.logDiffTime("Computing TF-IDF Matrix", tstart, addHeader = FALSE, verbose = verbose, logFile = logFile)
      mat <- as(Matrix::Diagonal(x=as.vector(idf)), "sparseMatrix") %*% mat
      
    }else if(LSI$LSIMethod == 2 | tolower(LSI$LSIMethod) == "log(tf-idf)"){
      
      #Adapted from Stuart et al.
      
      #IDF
      ArchR:::.logDiffTime("Computing Inverse Document Frequency", tstart, addHeader = FALSE, verbose = verbose, logFile = logFile)
      idf   <- as(LSI$nCol / LSI$rowSm, "sparseVector")
      
      #TF-IDF
      ArchR:::.logDiffTime("Computing TF-IDF Matrix", tstart, addHeader = FALSE, verbose = verbose, logFile = logFile)
      mat <- as(Matrix::Diagonal(x=as.vector(idf)), "sparseMatrix") %*% mat
      
      #Log transform TF-IDF
      mat@x <- log(mat@x * LSI$scaleTo + 1)  
      
    }else if(LSI$LSIMethod == 3 | tolower(LSI$LSIMethod) == "logtf-logidf"){
      
      #LogTF
      mat@x <- log(mat@x + 1)
      
      #LogIDF
      ArchR:::.logDiffTime("Computing Inverse Document Frequency", tstart, addHeader = FALSE, verbose = verbose, logFile = logFile)
      idf   <- as(log(1 + LSI$nCol / LSI$rowSm), "sparseVector")
      
      #TF-IDF
      ArchR:::.logDiffTime("Computing TF-IDF Matrix", tstart, addHeader = FALSE, verbose = verbose, logFile = logFile)
      mat <- as(Matrix::Diagonal(x=as.vector(idf)), "sparseMatrix") %*% mat
      
      
    }else if (LSI$LSIMethod == 4 ){
      # log TF-log IDF * norm total + log1p
      #LogTF
      mat@x <- log(mat@x + 1)
      
      #LogIDF
      ArchR:::.logDiffTime("Computing Inverse Document Frequency", tstart, addHeader = FALSE, verbose = verbose, logFile = logFile)
      idf   <- as(log(1 + LSI$nCol / LSI$rowSm), "sparseVector")
      
      #TF-IDF
      ArchR:::.logDiffTime("Computing TF-IDF Matrix", tstart, addHeader = FALSE, verbose = verbose, logFile = logFile)
      mat <- as(Matrix::Diagonal(x=as.vector(idf)), "sparseMatrix") %*% mat
      
      colSm <- Matrix::colSums(mat)
      mat@x <- mat@x / rep.int(colSm, Matrix::diff(mat@p))
      mat@x <- log(mat@x* LSI$scaleTo + 1)
    }else{
      
      stop("LSIMethod unrecognized please select valid method!")
      
    }
    
    gc()
    
    #Clean Up Matrix
    idxNA <- Matrix::which(is.na(mat),arr.ind=TRUE)
    if(length(idxNA) > 0){
      ArchR:::.logDiffTime(sprintf("Zeroing %s NA elements", length(idxNA)), tstart, addHeader = FALSE, verbose = verbose, logFile = logFile)
      mat[idxNA] <- 0
    }
    
    #Calc V
    ArchR:::.logDiffTime("Calculating V Matrix", tstart, addHeader = FALSE, verbose = verbose, logFile = logFile)
    V <- Matrix::t(mat) %*% LSI$svd$u %*% Matrix::diag(1/LSI$svd$d)
    
    #LSI Diagonal
    ArchR:::.logDiffTime("Computing Projected Coordinates", tstart, addHeader = FALSE, verbose = verbose, logFile = logFile)
    svdDiag <- matrix(0, nrow=LSI$nDimensions, ncol=LSI$nDimensions)
    diag(svdDiag) <- LSI$svd$d
    matSVD <- Matrix::t(svdDiag %*% Matrix::t(V))
    matSVD <- as.matrix(matSVD)
    rownames(matSVD) <- colnames(mat)
    colnames(matSVD) <- paste0("LSI",seq_len(ncol(matSVD)))
    
    if(returnModel){
      ArchR:::.logDiffTime("Calculating Re-Projected Matrix", tstart, addHeader = FALSE, verbose = verbose, logFile = logFile)
      X <- LSI$svd$u %*% diag(LSI$svd$d) %*% t(V)
      out <- list(matSVD = matSVD, V = V, X = X)
    }else{
      out <- matSVD
    }
    
    out
    
    
  }, error = function(e){
    
    errorList <- list(
      mat = mat,
      colSm = if(exists("colSm", inherits = FALSE)) colSm else "Error with colSm!",
      idf = if(exists("idf", inherits = FALSE)) idf else "Error with idf!",
      V = if(exists("V", inherits = FALSE)) V else "Error with V!",
      matSVD = if(exists("matSVD", inherits = FALSE)) matSVD else "Error with matSVD!"
    )
    
    ArchR:::.logError(e, fn = ".projectLSI", info = "", errorList = errorList, logFile = logFile)
    
  })
  
  out2
  
}

myaddFeatureMatrix <- function(
    input = NULL,
    features = NULL,
    matrixName = "FeatureMatrix",
    ceiling = 10^9, 
    binarize = FALSE,
    verbose = TRUE,
    threads = getArchRThreads(),
    parallelParam = NULL,
    force = TRUE,
    logFile = createLogFile("addFeatureMatrix")
){
  
  ArchR:::.validInput(input = input, name = "input", valid = c("ArchRProj", "character"))
  ArchR:::.validInput(input = features, name = "features", valid = c("GRanges"))
  ArchR:::.validInput(input = matrixName, name = "matrixName", valid = c("character"))
  ArchR:::.validInput(input = ceiling, name = "ceiling", valid = c("integer"))
  ArchR:::.validInput(input = binarize, name = "binarize", valid = c("boolean"))
  ArchR:::.validInput(input = verbose, name = "verbose", valid = c("boolean"))
  ArchR:::.validInput(input = threads, name = "threads", valid = c("integer"))
  ArchR:::.validInput(input = parallelParam, name = "parallelParam", valid = c("parallelparam", "null"))
  ArchR:::.validInput(input = force, name = "force", valid = c("boolean"))
  ArchR:::.validInput(input = logFile, name = "logFile", valid = c("character"))
  
  matrixName <- ArchR:::.isProtectedArray(matrixName)
  
  ArchR:::.startLogging(logFile = logFile)
  
  if(inherits(input, "ArchRProject")){
    ArrowFiles <- getArrowFiles(input)
    allCells <- rownames(getCellColData(input))
    outDir <- getOutputDirectory(input)
  }else if(inherits(input, "character")){
    outDir <- ""
    ArrowFiles <- input
    allCells <- NULL
  }else{
    stop("Error Unrecognized Input!")
  }
  if(!all(file.exists(ArrowFiles))){
    stop("Error Input Arrow Files do not all exist!")
  }
  
  #Add Index To Features
  features <- sort(sortSeqlevels(features), ignore.strand = TRUE)
  features <- split(features, seqnames(features))
  features <- lapply(features, function(x){
    mcols(x)$idx <- seq_along(x)
    return(x)
  })
  features <- Reduce("c",features)
  
  ArchR:::.logThis(features, "features", logFile = logFile)
  
  #Add args to list
  args <- mget(names(formals()),sys.frame(sys.nframe()))
  args$ArrowFiles <- ArrowFiles
  args$allCells <- allCells
  args$X <- seq_along(ArrowFiles)
  args$features <- features
  args$FUN <- myaddFeatureMatrix_helper
  args$registryDir <- file.path(outDir, "CountFeaturesRegistry")
  
  #Remove input from args
  args$input <- NULL
  
  #Run With Parallel or lapply
  outList <- ArchR:::.batchlapply(args)
  
  ArchR:::.endLogging(logFile = logFile)
  
  if(inherits(input, "ArchRProject")){
    return(input)
  }else{
    return(unlist(outList))
  }
  
}

myaddPeakMatrix <- function(
    ArchRProj = NULL,
    ceiling = 4, 
    binarize = FALSE,
    verbose = TRUE,
    threads = getArchRThreads(),
    parallelParam = NULL,
    force = TRUE,
    matrixName =  "PeakMatrix",
    logFile = createLogFile("addPeakMatrix")
){
  
  ArchR:::.validInput(input = ArchRProj, name = "ArchRProj", valid = c("ArchRProj"))
  ArchR:::.validInput(input = ceiling, name = "ceiling", valid = c("numeric"))
  ArchR:::.validInput(input = binarize, name = "binarize", valid = c("boolean"))
  ArchR:::.validInput(input = verbose, name = "verbose", valid = c("boolean"))
  ArchR:::.validInput(input = threads, name = "threads", valid = c("integer"))
  ArchR:::.validInput(input = parallelParam, name = "parallelParam", valid = c("parallelparam", "null"))
  ArchR:::.validInput(input = force, name = "force", valid = c("boolean"))
  ArchR:::.validInput(input = logFile, name = "logFile", valid = c("character"))
  
  if(is.null(ArchRProj@peakSet)){
    stop("No peakSet found in ArchRProject!")
  }
  
  ArrowFiles <- getArrowFiles(ArchRProj)
  allCells <- rownames(getCellColData(ArchRProj))
  outDir <- getOutputDirectory(ArchRProj)
  
  if(!all(file.exists(ArrowFiles))){
    stop("Error Input Arrow Files do not all exist!")
  }
  
  ArchR:::.startLogging(logFile = logFile)
  ArchR:::.logThis(ArchRProj@peakSet, "peakSet", logFile = logFile)
  
  #Add args to list
  args <- mget(names(formals()),sys.frame(sys.nframe()))#as.list(match.call())
  args$ArrowFiles <- ArrowFiles
  args$allCells <- allCells
  args$matrixName = matrixName
  args$features <- ArchRProj@peakSet
  args$X <- seq_along(ArrowFiles)
  args$FUN <- myaddFeatureMatrix_helper
  args$registryDir <- file.path(outDir, "CountPeaksRegistry")
  
  #Remove project from args
  args$ArchRProj <- NULL
  
  #Run With Parallel or lapply
  outList <- ArchR:::.batchlapply(args)
  
  readsInPeaks <- lapply(outList, function(x) x$RIP) %>% unlist
  FRIP <- lapply(outList, function(x) x$FRIP) %>% unlist
  ArchRProj <- addCellColData(ArchRProj, data = readsInPeaks, name = "ReadsInPeaks", names(readsInPeaks), force = force)
  ArchRProj <- addCellColData(ArchRProj, data = FRIP, name = "FRIP", names(readsInPeaks), force = force)
  
  ArchR:::.endLogging(logFile = logFile)
  
  return(ArchRProj)
  
}

myaddFeatureMatrix_helper <- function(
    i = NULL,
    ArrowFiles = NULL, 
    features = NULL,
    cellNames = NULL, 
    allCells = NULL,
    matrixName = "PeakMatrix", 
    ceiling = 4, 
    binarize = FALSE,
    tstart = NULL,
    subThreads = 1,
    force = FALSE,
    verbose = TRUE,
    logFile = NULL
){
  
  ArrowFile <- ArrowFiles[i]
  sampleName <- ArchR:::.sampleName(ArrowFile)
  
  #Check
  o <- h5closeAll()
  o <- ArchR:::.createArrowGroup(ArrowFile = ArrowFile, group = matrixName, force = force, logFile = logFile)
  
  if(!is.null(tstart)){
    tstart <- Sys.time()
  }
  
  #Get all cell ids before constructing matrix
  if(is.null(cellNames)){
    cellNames <- ArchR:::.availableCells(ArrowFile)
  }
  
  if(!is.null(allCells)){
    cellNames <- cellNames[cellNames %in% allCells]
  }
  
  dfParams <- data.frame(
    ceiling = ceiling, 
    binarize = binarize * 1,
    stringsAsFactors = FALSE)
  
  if("name" %in% colnames(mcols(features))){
    featureDF <- data.frame(
      seqnames = paste0(seqnames(features)), 
      idx = mcols(features)$idx, 
      start = start(features), 
      end = end(features), 
      name = mcols(features)$name,
      stringsAsFactors = FALSE)
  }else{
    featureDF <- data.frame(
      seqnames = paste0(seqnames(features)), 
      idx = mcols(features)$idx, 
      start = start(features), 
      end = end(features), 
      stringsAsFactors = FALSE)
  }
  ArchR:::.logThis(featureDF, "featureDF", logFile = logFile)
  
  ######################################
  # Initialize SP Mat Group
  ######################################
  if(binarize){
    Class <- "binary"
    Units <- "BinarizedCounts"
  }else{
    Class <- "integer"
    Units <- "Counts"
  }
  o <- ArchR:::.initializeMat(
    ArrowFile = ArrowFile,
    Group = matrixName,
    Class = Class,
    Units = Units,
    cellNames = cellNames,
    params = dfParams,
    featureDF = featureDF,
    force = force
  )
  
  ######################################
  # Add To SP Mat Group
  ######################################
  strand(features) <- "*" #lets make sure nothing strand related occurs
  uniqueChr <- as.character(unique(seqnames(features)@values))
  insertionsInPeaks <- rep(0, length(cellNames))
  names(insertionsInPeaks) <- cellNames
  totalInsertions <- insertionsInPeaks
  
  for(z in seq_along(uniqueChr)){
    
    prefix <- sprintf("Chr %s (%s of %s)!", uniqueChr[z], z, length(uniqueChr))
    
    o <- tryCatch({
      
      o <- h5closeAll()
      chr <- uniqueChr[z]
      featurez <- features[BiocGenerics::which(seqnames(features)==chr)]
      ArchR:::.logDiffTime(sprintf("Adding %s to %s for Chr (%s of %s)!", sampleName, matrixName, z, length(uniqueChr)), tstart, verbose = verbose, logFile = logFile)
      
      #Read in Fragments
      fragments <- ArchR:::.getFragsFromArrow(ArrowFile, chr = chr, out = "IRanges", cellNames = cellNames)
      tabFrags <- table(mcols(fragments)$RG)
      
      #Count Left Insertion
      temp <- IRanges(start = start(fragments), width = 1)
      stopifnot(length(temp) == length(fragments))
      oleft <- findOverlaps(ranges(featurez), temp)
      oleft <- DataFrame(queryHits=Rle(queryHits(oleft)), subjectHits = subjectHits(oleft))
      
      #Count Right Insertion
      temp <- IRanges(start = end(fragments), width = 1)
      stopifnot(length(temp) == length(fragments))
      oright <- findOverlaps(ranges(featurez), temp)
      oright <- DataFrame(queryHits=Rle(queryHits(oright)), subjectHits = subjectHits(oright))
      remove(temp)
      ################ remove same fragment 
      z=rbind(oleft,oright)
      remove(oleft,oright)
      z=unique(z)
      
      
      if(F){
        #Feature Idx
        oleft$queryHits@values <- mcols(featurez)$idx[oleft$queryHits@values] 
        oright$queryHits@values <- mcols(featurez)$idx[oright$queryHits@values]
        #Correct to RG ID
        oleft$subjectHits <- as.integer(BiocGenerics::match(mcols(fragments)$RG[oleft$subjectHits], cellNames))
        oright$subjectHits <- as.integer(BiocGenerics::match(mcols(fragments)$RG[oright$subjectHits], cellNames))
        #Create Sparse Matrix
        mat <- Matrix::sparseMatrix(
          i = c( oleft$queryHits, oright$queryHits ),
          j = c( oleft$subjectHits, oright$subjectHits ),
          x = rep(1, nrow(oleft) + nrow(oright)),
          dims = c(max(mcols(featurez)$idx), length(cellNames))
        )
      }
      #
      #
      
      
      #
      #
      z$queryHits@values<- mcols(featurez)$idx[z$queryHits@values]
      z$subjectHits <- as.integer(BiocGenerics::match(mcols(fragments)$RG[z$subjectHits], cellNames))
      remove(fragments)

      mat <- Matrix::sparseMatrix(
        i = z$queryHits,
        j = z$subjectHits,
        x = rep(1, nrow(z)),
        dims = c(max(mcols(featurez)$idx), length(cellNames))
      )
      colnames(mat) <- cellNames
      
      #Compute total reads in Peak
      totalInsertions[names(tabFrags)] <- totalInsertions[names(tabFrags)] + 2 * tabFrags
      insertionsInPeaks <- insertionsInPeaks + Matrix::colSums(mat)
      
      #Ceiling
      if(!is.null(ceiling)){
        mat@x[mat@x > ceiling] <- ceiling
      }
      if(binarize){
        mat@x[mat@x > 0] <- 1
      }
      
      #Write sparseMatrix to Arrow File!
      o <- ArchR:::.addMatToArrow(
        mat = mat, 
        ArrowFile = ArrowFile, 
        Group = paste0(matrixName,"/", chr), 
        binarize = binarize,
        addColSums = TRUE,
        addRowSums = TRUE
      ) 
      gc()
      
      o
      
    }, error = function(e){
      
      errorList <- list(
        ArrowFile = ArrowFile,
        z = z,
        chr = uniqueChr[z],
        oleft = if(exists("oleft", inherits = FALSE)) oleft else "oleft",
        oright = if(exists("oright", inherits = FALSE)) oright else "oright",
        mat = if(exists("mat", inherits = FALSE)) mat else "mat",
        totalInsertions = if(exists("totalInsertions", inherits = FALSE)) totalInsertions else "totalInsertions",
        insertionsInPeaks = if(exists("insertionsInPeaks", inherits = FALSE)) insertionsInPeaks else "insertionsInPeaks",
        fragments = tryCatch({
          ArchR:::.getFragsFromArrow(ArrowFile, chr = chr, out = "IRanges", cellNames = cellNames)
        }, error = function(f){
          "fragments"
        })
      )
      
      ArchR:::.logError(e, fn = "myaddFeatureMatrix", info = prefix, errorList = errorList, logFile = logFile)
      
    })
    
  }
  
  out <- list(ArrowFile = ArrowFile, RIP = insertionsInPeaks, FRIP = insertionsInPeaks / totalInsertions)
  
  return(out)
  
}
