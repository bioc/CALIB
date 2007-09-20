setClass("RGList_CALIB",contains = "list")
setClass("SpikeList",contains = "list")
setClass("ParameterList",contains = "list")
   
read.rg<-function (files=NULL, source="generic", path=NULL, ext=NULL, 
    names=NULL, columns=NULL, other.columns=NULL, annotation=NULL, 
    wt.fun=NULL, verbose=TRUE, sep="\t", quote=NULL, 
    DEBUG=FALSE, ...) 
{
    debugVars <- function(DEBUG) {
        if (DEBUG) {
            print(paste("FILE:", fullname))
            print(paste("SKIP:", skip))
            print("SELECTED COLUMNS")
            print(colClasses)
        }
    }
    if (missing(files)) {
        if (is.null(ext)) 
            stop("Must specify input files")
        else {
            extregex <- paste("\\.", ext, "$", sep = "")
            files <- dir(path = ifelse(is.null(path), ".", path), 
                pattern = extregex)
            files <- sub(extregex, "", files)
        }
    }
    source <- match.arg(source, c("generic", "agilent", "arrayvision", 
        "bluefuse", "genepix", "genepix.median", "genepix.custom", 
        "imagene", "quantarray", "smd.old", "smd", "spot", "spot.close.open"))
    source2 <- strsplit(source, split = ".", fixed = TRUE)[[1]][1]
    if (is.null(quote)) 
        if (source == "agilent") 
            quote <- ""
        else quote <- "\""
    if (source2 == "imagene") 
        return(read.imagene(files = files, path = path, ext = ext, 
            names = names, columns = columns, wt.fun = wt.fun, 
            verbose = verbose, sep = sep, quote = quote, ...))
    slides <- as.vector(as.character(files))
    if (!is.null(ext)) 
        slides <- paste(slides, ext, sep = ".")
    nslides <- length(slides)
    if (is.null(names)) 
        names <- removeExt(files)
    if (is.null(columns)) {
        if (source2 == "generic") 
            stop("must specify columns for generic input")
        columns <- switch(source, agilent = list(G = "gMeanSignal", 
            Gb = "gBGMedianSignal", R = "rMeanSignal", Rb = "rBGMedianSignal"), 
            bluefuse = list(G = "AMPCH1", R = "AMPCH2"), genepix = list(R = "F635 Mean", 
                G = "F532 Mean", Rb = "B635 Median", Gb = "B532 Median",RArea = "Dia",GArea = "Dia"), 
            genepix.median = list(R = "F635 Median", G = "F532 Median", 
                Rb = "B635 Median", Gb = "B532 Median",RArea = "Dia",GArea = "Dia"), genepix.custom = list(R = "F635 Mean", 
                G = "F532 Mean", Rb = "B635", Gb = "B532"), quantarray = list(R = "ch2 Intensity", 
                G = "ch1 Intensity", Rb = "ch2 Background", Gb = "ch1 Background",RArea = "ch2 Diameter", GArea = "ch1 Diameter"), 
            smd.old = list(G = "CH1I_MEAN", Gb = "CH1B_MEDIAN", 
                R = "CH2I_MEAN", Rb = "CH2B_MEDIAN"), smd = list(G = "Ch1 Intensity (Mean)", 
                Gb = "Ch1 Background (Median)", R = "Ch2 Intensity (Mean)", 
                Rb = "Ch2 Background (Median)"), spot = list(R = "Rmean", 
                G = "Gmean", Rb = "morphR", Gb = "morphG"), spot.close.open = list(R = "Rmean", 
                G = "Gmean", Rb = "morphR.close.open", Gb = "morphG.close.open"), 
            NULL)
    }else {
        if (!is.list(columns)) 
            stop("columns must be a list")
        names(columns)[names(columns) == "Gf"] <- "G"
        names(columns)[names(columns) == "Rf"] <- "R"
        if (is.null(columns$G) || is.null(columns$R)) 
            stop("columns must specify foreground G and R")
        if (!all(names(columns) %in% c("G", "R", "Gb", "Rb","GArea","RArea"))) 
            warning("non-standard columns specified")
    }
    cnames <- names(columns)
    if (is.null(annotation)) 
        annotation <- switch(source, agilent = c("Row", "Col", 
            "Start", "Sequence", "SwissProt", "GenBank", "Primate", 
            "GenPept", "ProbeUID", "ControlType", "ProbeName", 
            "GeneName", "SystematicName", "Description"), bluefuse = c("ROW", 
            "COL", "SUBGRIDROW", "SUBGRIDCOL", "BLOCK", "NAME", 
            "ID"), genepix = , genepix.median = , genepix.custom = c("Block", 
            "Row", "Column", "ID", "Name"), smd = c("Spot", "Clone ID", 
            "Gene Symbol", "Gene Name", "Cluster ID", "Accession", 
            "Preferred name", "Locuslink ID", "Name", "Sequence Type", 
            "X Grid Coordinate (within sector)", "Y Grid Coordinate (within sector)", 
            "Sector", "Failed", "Plate Number", "Plate Row", 
            "Plate Column", "Clone Source", "Is Verified", "Is Contaminated", 
            "Luid"), smd.old = c("SPOT", "NAME", "Clone ID", 
            "Gene Symbol", "Gene Name", "Cluster ID", "Accession", 
            "Preferred name", "SUID"), NULL)
    fullname <- slides[1]
    if (!is.null(path)) 
        fullname <- file.path(path, fullname)
    switch(source2, quantarray = {
        firstfield <- scan(fullname, what = "", sep = "\t", flush = TRUE, 
            quiet = TRUE, blank.lines.skip = FALSE, multi.line = FALSE, 
            allowEscapes = FALSE)
        skip <- grep("Begin Data", firstfield)
        if (length(skip) == 0) 
            stop("Cannot find \"Begin Data\" in image output file")
        nspots <- grep("End Data", firstfield) - skip - 2
        allcnames <- scan(fullname, what = character(1), sep = sep, 
            quote = quote, skip = skip, nlines = 1, quiet = TRUE, 
            allowEscapes = FALSE)
        colClasses <- getColClasses(allcnames, annotation, columns, 
            other.columns, wt.fun)
        debugVars(DEBUG)
        obj <- read.table(fullname, skip = skip, header = TRUE, 
            sep = sep, quote = quote, as.is = TRUE, check.names = FALSE, 
            comment.char = "", nrows = nspots, colClasses = colClasses, 
            ...)
    }, arrayvision = {
        skip <- 1
        cn <- scan(fullname, what = "", sep = sep, quote = quote, 
            skip = 1, nlines = 1, quiet = TRUE, allowEscape = FALSE)
        fg <- grep(" Dens - ", cn)
        if (length(fg) != 2) 
            stop(paste("Cannot find foreground columns in", fullname))
        bg <- grep("Bkgd", cn)
        if (length(fg) != 2) 
            stop(paste("Cannot find background columns in", fullname))
        columns <- list(R = fg[1], Rb = bg[1], G = fg[2], Gb = bg[2])
        allcnames <- scan(fullname, what = character(1), sep = sep, 
            quote = quote, skip = skip, nlines = 1, quiet = TRUE, 
            allowEscapes = FALSE)
        colClasses <- getColClasses(allcnames, annotation, columns, 
            other.columns, wt.fun)
        debugVars(DEBUG)
        obj <- read.table(fullname, skip = skip, header = TRUE, 
            sep = sep, quote = quote, as.is = TRUE, check.names = FALSE, 
            comment.char = "", colClasses = colClasses, ...)
        nspots <- nrow(obj)
    }, bluefuse = {
        skip <- readBlueFuseHeader(fullname)$NHeaderRecords
        allcnames <- scan(fullname, what = character(1), sep = sep, 
            quote = quote, skip = skip, nlines = 1, quiet = TRUE, 
            allowEscapes = FALSE)
        colClasses <- getColClasses(allcnames, annotation, columns, 
            other.columns, wt.fun)
        debugVars(DEBUG)
        obj <- read.table(fullname, skip = skip, header = TRUE, 
            sep = sep, quote = quote, as.is = TRUE, check.names = FALSE, 
            comment.char = "", fill = TRUE, colClasses = colClasses, 
            ...)
        nspots <- nrow(obj)
    }, genepix = {
        h <- readGPRHeader(fullname)
        if (verbose && source == "genepix.custom") 
            cat("Custom background:", h$Background, "\n")
        skip <- h$NHeaderRecords
        allcnames <- scan(fullname, what = character(1), sep = sep, 
            quote = quote, skip = skip, nlines = 1, quiet = TRUE, 
            allowEscapes = FALSE)
        colClasses <- getColClasses(allcnames, annotation, columns, 
            other.columns, wt.fun, "X")
        debugVars(DEBUG)
        obj <- read.table(fullname, skip = skip, header = TRUE, 
            sep = sep, quote = quote, as.is = TRUE, check.names = FALSE, 
            comment.char = "", fill = TRUE, colClasses = colClasses, 
            ...)
        nspots <- nrow(obj)
    }, smd = {
        skip <- readSMDHeader(fullname)$NHeaderRecords
        allcnames <- scan(fullname, what = character(1), sep = sep, 
            quote = quote, skip = skip, nlines = 1, quiet = TRUE, 
            allowEscapes = FALSE)
        colClasses <- getColClasses(allcnames, annotation, columns, 
            other.columns, wt.fun)
        debugVars(DEBUG)
        obj <- read.table(fullname, skip = skip, header = TRUE, 
            sep = sep, quote = quote, as.is = TRUE, check.names = FALSE, 
            comment.char = "", fill = TRUE, colClasses = colClasses, 
            ...)
        nspots <- nrow(obj)
    }, {
        skip <- grep(protectMetachar(columns$R), readLines(fullname, 
            n = 80)) - 1
        if (length(skip) == 0) 
            stop("Cannot find column heading in image output file")
        else skip <- skip[1]
        allcnames <- scan(fullname, what = character(1), sep = sep, 
            quote = quote, skip = skip, nlines = 1, quiet = TRUE, 
            allowEscapes = FALSE)
        colClasses <- getColClasses(allcnames, annotation, columns, 
            other.columns, wt.fun)
        debugVars(DEBUG)
        obj <- read.table(fullname, skip = skip, header = TRUE, 
            sep = sep, quote = quote, as.is = TRUE, check.names = FALSE, 
            comment.char = "", fill = TRUE, colClasses = colClasses,
            ...)
        nspots <- nrow(obj)
    })
    Y <- matrix(NA, nspots, nslides)
    colnames(Y) <- names
    RG <- columns
    for (a in cnames) RG[[a]] <- Y
    if (!is.null(wt.fun)) 
        RG$weights <- Y
    RG$targets <- data.frame(FileName = I(files), row.names = names)
    if (!is.null(annotation)) {
        j <- match(annotation, colnames(obj), 0)
        if (any(j > 0)) 
            RG$genes <- data.frame(obj[, j, drop = FALSE], check.names = FALSE)
    }
    RG$source <- source
    if (source2 == "agilent") {
        if (!is.null(RG$genes$Row) && !is.null(RG$genes$Col)) {
            nr <- length(unique(RG$genes$Row))
            nc <- length(unique(RG$genes$Col))
            if (nspots == nr * nc) 
                RG$printer <- list(ngrid.r = 1, ngrid.c = 1, 
                  nspot.r = nr, nspot.c = nc)
        }
    }
    if (source2 == "genepix") {
        if (!is.null(RG$genes$Block) && !is.null(RG$genes$Row) && 
            !is.null(RG$genes$Column)) {
            RG$printer <- getLayout(RG$genes, guessdups = TRUE)
            nblocks <- RG$printer$ngrid.r * RG$printer$ngrid.c
            if (!is.na(nblocks) && (nblocks > 1) && !is.null(obj$X)) {
                blocksize <- RG$printer$nspot.r * RG$printer$nspot.c
                i <- (1:(nblocks - 1)) * blocksize
                ngrid.r <- sum(obj$X[i] > obj$X[i + 1]) + 1
                if (!is.na(ngrid.r) && nblocks%%ngrid.r == 0) {
                  RG$printer$ngrid.r <- ngrid.r
                  RG$printer$ngrid.c <- nblocks/ngrid.r
                }
                else {
                  warning("Can't determine number of grid rows")
                  RG$printer$ngrid.r <- RG$ngrid.c <- NA
                }
            }
        }
    }
    if (!is.null(other.columns)) {
        other.columns <- as.character(other.columns)
        j <- match(other.columns, colnames(obj), 0)
        if (any(j > 0)) {
            other.columns <- colnames(obj)[j]
            RG$other <- list()
            for (j in other.columns) RG$other[[j]] <- Y
        }
        else {
            other.columns <- NULL
        }
    }
    for (i in 1:nslides) {
        if (i > 1) {
            fullname <- slides[i]
            if (!is.null(path)) 
                fullname <- file.path(path, fullname)
            switch(source2, quantarray = {
                firstfield <- scan(fullname, what = "", sep = "\t", 
                  flush = TRUE, quiet = TRUE, blank.lines.skip = FALSE, 
                  multi.line = FALSE, allowEscapes = FALSE, nlines = 100)
                skip <- grep("Begin Data", firstfield)
            }, arrayvision = {
                skip <- 1
            }, bluefuse = {
                skip <- readBlueFuseHeader(fullname)$NHeaderRecords
            }, genepix = {
                skip <- readGPRHeader(fullname)$NHeaderRecords
            }, smd = {
                skip <- readSMDHeader(fullname)$NHeaderRecords
            }, {
                skip <- grep(protectMetachar(columns$R), readLines(fullname, 
                  n = 80)) - 1
            })
            if (verbose && source == "genepix.custom") 
                cat("Custom background:", h$Background, "\n")
            allcnames <- scan(fullname, what = character(1), 
                sep = sep, quote = quote, skip = skip, nlines = 1, 
                quiet = TRUE, allowEscapes = FALSE)
            colClasses <- getColClasses(allcnames, annotation, 
                columns, other.columns, wt.fun)
            debugVars(DEBUG)
            obj <- read.table(fullname, skip = skip, header = TRUE, 
                sep = sep, as.is = TRUE, quote = quote, check.names = FALSE, 
                comment.char = "", fill = TRUE, nrows = nspots, 
                colClasses = colClasses,...)
        }
        for (a in cnames) RG[[a]][, i] <- obj[, columns[[a]]]
        if (!is.null(wt.fun)) 
            RG$weights[, i] <- wt.fun(obj)
        if (!is.null(other.columns)) 
            for (j in other.columns) {
                RG$other[[j]][, i] <- obj[, j]
            }
        if (verbose) 
            cat(paste("Read", fullname, "\n"))
    }
    new("RGList_CALIB", RG)
}

read.spike<-function(RG,file=NULL,path=NULL,ext=NULL,sep="\t",conccol,regexpcol,different=FALSE,...)
{
   if (missing(RG))
   {
      stop("Must specify RG data")
   }
   if (is.null(RG$genes))
   {
      stop("Must specify gene information.")
   }
   if (is.null(RG$genes$Status))
   {
      stop("Must specify gene status.")
   }
   stname<-intersect(RG$genes$Status,c("Negative","Calibration","Ratio"))
   if (length(stname)!= 3)
   {
      stop("non-standard status value specified")
   }
   if (missing(file)) 
   {
      if (is.null(ext)) 
          stop("Must specify concentration files")
      else 
      {
          extregex <- paste("\\.", ext, "$", sep = "")
          file <- dir(path = ifelse(is.null(path), ".", path), pattern = extregex)
          file <- sub(extregex, "", files)
      }
   }
   if (missing(conccol))
   {
      conccol<-list(RConc="REDCONC",GConc="GREENCONC")
   }else{
      if (!is.list(conccol))
      {
         stop("Concentration column must be a list.")
      }
   } 
   filename<-as.vector(as.character(file))
   if (!is.null(ext))
   {
      filename<-paste(filename,ext,sep=".")
   }
   filenum<-length(filename)
   arraynum<-nrow(RG$targets)
   if (arraynum < 2)
   {
      RG$R <- as.matrix(RG$R)
      RG$G <- as.matrix(RG$G)
      RG$Rb <- as.matrix(RG$Rb)
      RG$Gb <- as.matrix(RG$Gb)
      RG$RArea <- as.matrix(RG$RArea)
      RG$GArea <- as.matrix(RG$GArea)
   }
   spotnum<-nrow(RG$R)
   if(!different)
   {
      if (filenum!=1)
      {
         stop("Only need one input file.")     
      }else{
         fullname<-file.path(path,filename)         
      }
   }else{
      if (filenum==1)
      {
         stop("Must specify more than one input files.")
      }
      if (filenum!=arraynum)
      {
         stop("The number of concentration files and the number of arrays don't match.")
      }
      fullname<-file.path(path,filename[1])
   }
   concinfo<-read.table(fullname,header=TRUE,as.is=TRUE,sep=sep,quote="\"",...)
   names<-colnames(concinfo)
   genes<-RG$genes
   genenames<-colnames(genes)
   concnames<-names(conccol)
   if (length(intersect(concnames,c("RConc","GConc")))!= 2)
   {
      stop("non-standard concentration column names.")
   }
   if (is.null(names) || is.null(genenames)) 
   {
      stop("concentration files and genes must have column names")
   }
   if (!missing(regexpcol) && is.numeric(regexpcol))
   {
      regexpcol<-names[regexpcol]
   }
   if (missing(regexpcol))
   {
      commoncol<-(names %in% genenames)
      if (!length(commoncol))
      {
         stop("concentration files and genes should have common columns")
      }
       regexpcol<-names[commoncol]
   }else{
      if (!all(is.element(regexpcol,names)& is.element(regexpcol,genenames)))
      {
         stop("concentration files and genes must both contain regexpcol columns") 
      }
   }
   for (i in regexpcol)
   {
      concinfo[,i]<-paste("^", gsub("\\*", "\\.*", concinfo[,i]), "$", sep = "")
      genes[,i]<-as.character(genes[,i])
   }
   ntypes<-nrow(concinfo)
   nregexp<-length(regexpcol)
   spikenum<-0
   allsel<-NULL
   conc<-conccol
   for (a in concnames )
   {
      conc[[a]]<-rep(NA,spotnum)
   }
   for(i in 1:ntypes)
   {
      sel<-grep(concinfo[i, regexpcol[1]],genes[, regexpcol[1]])
      if (nregexp>1)
      {
         for (j in regexpcol[-1])
         {
             sel<-intersect(sel,grep(concinfo[i, j],genes[, j]))
         }
       }
       for (a in concnames)
       {
         conc[[a]][sel]<-concinfo[i,conccol[[a]]]
       }
       spikenum<-spikenum + length(sel)
       allsel<-union(allsel,sel)
   } 
   fields<-list(R=0,G=0,Rb=0,Gb=0,RArea=0,GArea=0,RConc=0,GConc=0)
   fieldnames<-names(fields)
   spike<-fields
   for (a in setdiff(fieldnames,concnames))
   {
      if (arraynum < 2)
      {
         spike[[a]]<-RG[[a]][allsel] 
      }else{
         spike[[a]]<-RG[[a]][allsel,]
      }
   }
   allconc<-conccol
   for (a in concnames)
   {
      allconc[[a]]<-matrix(NA,spotnum,arraynum)
   }      
   if (!different)
   {
      for (a in concnames)
      {
         for (i in 1:arraynum)
         {
            allconc[[a]][,i]<-conc[[a]]
         }
      }
   }else{
      for (a in concnames)
      {          
         allconc[[a]][,1]<-conc[[a]]
      }
      for (l in 2:arraynum)
      {
         fullname<-file.path(path,filename[l])
         concinfo<-read.table(fullname,header=TRUE,as.is=TRUE,sep=sep,quote="\"",...)
         names<-colnames(concinfo)
         if (is.null(names) || is.null(genenames)) 
         {
            stop(paste("concentration files",i," and genes must have column names",sep=""))
         }
         for (i in regexpcol)
         {
            concinfo[,i]<-paste("^", gsub("\\*", "\\.*", concinfo[,i]), "$", sep = "")
         }
         ntypes<-nrow(concinfo)
         conc<-conccol
         for (a in concnames )
         {
            conc[[a]]<-rep(NA,spotnum)
         }
         for(i in 1:ntypes)
         {
             sel<-grep(concinfo[i, regexpcol[1]],genes[, regexpcol[1]])
             if (nregexp>1)
             {
                for (j in regexpcol[-1])
                {
                   sel<-intersect(sel,grep(concinfo[i, j],genes[, j]))
                }
             }
             for (a in concnames)
             {
                conc[[a]][sel]<-concinfo[i,conccol[[a]]]
             }
          }
          for (a in concnames)
          {
              allconc[[a]][,l]<-conc[[a]]
          }
        }
   }
   for (a in concnames)
   {
      colnames(allconc[[a]])<-colnames(RG$R)
      spike[[a]]<-allconc[[a]][allsel,]
   }        
   spike$genes<-data.frame(RG$genes[allsel,])
   new("SpikeList",spike)  
}  

estimateParameter<-function(spike,RG,bc=FALSE,area=TRUE,errormodel="M")
{
   if(!is(spike,"SpikeList"))
   {
      stop("Object of SpikeList class is required.")
   }
   if (!is(RG,"RGList_CALIB"))
   {
      stop("Object of RGList_CALIB class is required.")
   }
   if (!is.element(errormodel,c("A","M")))
   {
      stop("errormodel is expected to be A or M")
   }
   arraynum<-nrow(RG$targets)
   if (arraynum < 2)
   {
      RG$R <- as.matrix(RG$R)
      RG$G <- as.matrix(RG$G)
      RG$Rb <- as.matrix(RG$Rb)
      RG$Gb <- as.matrix(RG$Gb)
      RG$RArea <- as.matrix(RG$RArea)
      RG$GArea <- as.matrix(RG$GArea)
      spike$R <- as.matrix(spike$R)
      spike$G <- as.matrix(spike$G)
      spike$Rb <- as.matrix(spike$Rb)
      spike$Gb <- as.matrix(spike$Gb)
      spike$RArea <- as.matrix(spike$RArea)
      spike$GArea <- as.matrix(spike$GArea)
      spike$RConc <- as.matrix(spike$RConc)
      spike$GConc <- as.matrix(spike$GConc)   
   }
   if (bc)
   {
      if (area)
      {
         SRinten <- (spike$R - spike$Rb) * spike$RArea
         SGinten <- (spike$G - spike$Gb) * spike$GArea
         Rinten <- (RG$R - RG$Rb) * RG$RArea
         Ginten <- (RG$G - RG$Gb) * RG$GArea
      }else{
         SRinten <- spike$R - spike$Rb
         SGinten <- spike$G - spike$Gb
         Rinten <- RG$R - RG$Rb
         Ginten <- RG$G - RG$Gb
      }
   }else{
      if (area)
      {
         SRinten <- spike$R * spike$RArea
         SGinten <- spike$G * spike$GArea
         Rinten <- RG$R * RG$RArea
         Ginten <- RG$G * RG$GArea

      }else{
         SRinten <- spike$R
         SGinten <- spike$G
         Rinten <- RG$R
         Ginten <- RG$G

      }
   }
   SRinten <- pmax(SRinten, 0.5, na.rm=TRUE)
   SGinten <- pmax(SGinten, 0.5, na.rm=TRUE)
   SRconc <- spike$RConc
   SGconc <- spike$GConc
   meaninten <- max(log10(mean(Rinten)),log10(mean(Ginten)))
   meanconc <- max(log10(mean(SRconc)),log10(mean(SGconc)))
   mus <- 10^(ceiling(mean(meaninten,meanconc)))

   type <- spike$genes$Status
   storage.mode(SRinten)<-storage.mode(SGinten)<-"double"
   storage.mode(SRconc)<-storage.mode(SGconc)<-"double"
   storage.mode(mus)<- "double"
   storage.mode(type)<-storage.mode(errormodel)<-"character"
 
   maxinten <- max(max(Rinten[,1]),max(Ginten[,1]))
   storage.mode(maxinten)<-"double"

   par<-.Call("estimation",SRinten[,1],SGinten[,1],SRconc[,1],SGconc[,1],type,maxinten,errormodel,mus,PACKAGE="CALIB")
   parameter<-par
   fieldnames<-names(parameter)
   for (a in fieldnames)
   {
       Y<-matrix(NA,length(par[[a]]),arraynum)
       parameter[[a]]<-Y
       colnames(parameter[[a]])<-colnames(spike$R)
   }
   for (a in fieldnames)
   {
       parameter[[a]][,1] <- par[[a]]
   }
   if (arraynum > 1)
   {
      for (i in 2:arraynum)
      {
         maxinten <- max(max(Rinten[,i]),max(Ginten[,i]))
         storage.mode(maxinten)<-"double"
         par<-.Call("estimation",SRinten[,i],SGinten[,i],SRconc[,i],SGconc[,i],type,maxinten,errormodel,mus,PACKAGE="CALIB")
         for (a in fieldnames)
         {
             parameter[[a]][,i] <- par[[a]]
         }
      }
   }  
   parameter$Method$BC<-bc
   parameter$Method$Area<-area
   parameter$ErrorModel<-errormodel
   
   negindex <- spike$genes$Status %in% "Negative"
   index <- !negindex
   parameter$genes<-spike$genes[index,]

   new("ParameterList",parameter)
}

adjustP2<-function(RG,parameter,arrayindex=arrayindex,colorindex=colorindex)
{
   if (missing(RG))
   {
      stop("Must specify RG data.")
   }
   if (missing(parameter))
   {
      stop("Must specify parameter.")
   }
   if (is.null(RG))
   {
      stop("Must specify RG data.")
   }
   if (is.null(parameter))
   {
      stop("Must specify parameter.")
   }
   if (!is(RG,"RGList_CALIB"))
   {
      stop("Object of RGList_CALIB class is required.")
   }
   if (!is(parameter,"ParameterList"))
   {
      stop("Object of ParameterList class is required.")
   }
   if (!is.numeric(arrayindex))
   {
      stop("arrayindex should be numeric.")
   }
   if (!is.numeric(colorindex))
   {
      stop("colorindex should be numeric.")
   }
   lengtharray<-length(arrayindex)
   lengthcolor<-length(colorindex)
   if (lengtharray != lengthcolor)
   {
      stop("arrayindex and colorindex should have the same length")
   }
   for (i in 1:lengthcolor)
   {
      if (!is.element(colorindex[i],c(1,2)))
      {
         stop("colorindex should be either 1 or 2.")
      }
   }
   arraynum<-nrow(RG$targets)
   if (arraynum < 2)
   {
      RG$R <- as.matrix(RG$R)
      RG$G <- as.matrix(RG$G)
      RG$Rb <- as.matrix(RG$Rb)
      RG$Gb <- as.matrix(RG$Gb)
      RG$RArea <- as.matrix(RG$RArea)
      RG$GArea <- as.matrix(RG$GArea)      
   }
   bc <- parameter$Method$BC
   area <- parameter$Method$Area
   if (bc)
   {
      if (area)
      {
         Rinten <- (RG$R - RG$Rb) * RG$RArea
         Ginten <- (RG$G - RG$Gb) * RG$GArea
      }else{
         Rinten <- RG$R - RG$Rb
         Ginten <- RG$G - RG$Gb
      }
   }else{
      if (area)
      {
         Rinten <- RG$R * RG$RArea
         Ginten <- RG$G * RG$GArea
      }else{
         Rinten <- RG$R
         Ginten <- RG$G
      }
   }
   
   index <- 1
   adjust <- matrix(1,2,ncol(parameter$P2))
   parameter$AdjustFactor <- adjust

   for (i in arrayindex)
   {
      which <- colorindex[index]
      
      ka<-parameter$Ka[,i]
      mus<-parameter$MuS[,i]
      p1<-parameter$P1[,i]
      p2<-parameter$P2[,i]
      sigmaa<-parameter$SigmaA[,i]
      sigmam<-parameter$SigmaM[,i]
      sigmas<-parameter$SigmaS[,i]
      parlist<-list(Ka=ka,MuS=mus,P1=p1,P2=p2,SigmaA=sigmaa,SigmaM=sigmam,SigmaS=sigmas)
      par<-new("ParameterList",parlist)

      storage.mode(Rinten[,i]) <- storage.mode(Ginten[,i])<-"double"
      storage.mode(par)<-"list"
      storage.mode(which) <- "integer"
      newP2 <- .Call("adjustment",Rinten[,i],Ginten[,i],par,which,PACKAGE="CALIB")
      adjustfactor<- newP2 / parameter$P2[which,i]
      parameter$AdjustFactor[which,i] <- adjustfactor  
      index <- index + 1
   }

   return(parameter)
}

normalizeData<-function(RG,parameter,array=array,condition=condition,dye=dye,cloneid=cloneid,idcol=idcol)
{
   if (missing(RG))
   {
      stop("Must specify RG data.")
   }
   if (missing(parameter))
   {
      stop("Must specify parameter.")
   }
   if (missing(array))
   {
      stop("Must specify array design.")
   }
   if (missing(condition))
   {
      stop("Must specify condition design.")
   }
   if (missing(dye))
   {
      stop("Must specify dye design.")
   }
   if (is.null(RG))
   {
      stop("Must specify RG data.")
   }
   if (is.null(parameter))
   {
      stop("Must specify parameter.")
   }
   if (is.null(array))
   {
      stop("Must specify array design.")
   }
   if (is.null(condition))
   {
      stop("Must specify condition design.")
   }
   if (is.null(dye))
   {
      stop("Must specify dye design.")
   }
   if (!is(RG,"RGList_CALIB"))
   {
      stop("Object of RGList_CALIB class is required.")
   }
   if (!is(parameter,"ParameterList"))
   {
      stop("Object of ParameterList class is required.")
   }
   if (!is(array,"numeric"))
   {
      stop("array must be a numeric vector.")
   }
   if (!is(condition,"numeric"))
   {
      stop("condition must be a numeric vector.")
   }
   if (!is(dye,"numeric"))
   {
      stop("dye must be a numeric vector.")
   }

   arraynum<-nrow(RG$targets)
   if (arraynum < 2)
   {
      RG$R <- as.matrix(RG$R)
      RG$G <- as.matrix(RG$G)
      RG$Rb <- as.matrix(RG$Rb)
      RG$Gb <- as.matrix(RG$Gb)
      RG$RArea <- as.matrix(RG$RArea)
      RG$GArea <- as.matrix(RG$GArea)
   }
      
   ### define the clone id to be normalized. 
     ### if the parameter "cloneid" is missing, it means that all clones are to be normalized.
     ### if "cloneid" is not missing, normalize the specified ones.
   ### parameter "idcol" is needed.
     ### if missing (cloneid), use "idcol" to find out every gene
     ### if not missing(cloneid), use "idcol" to check whether those ids are valid.
 
   if (missing(idcol))
   {
      stop("Must specify clone id column.")
   }
   if (is.null(idcol))
   {
      stop("Must specify clone id column.")
   }
   if (! is.element(idcol,colnames(RG$genes)))
   { 
      stop("user specified clone id column is not in RG$genes columns.")
   }
   if (missing(cloneid))
   {      
      cloneid<-RG$genes[,idcol]
   }else
   {
      valid<-is.element(cloneid,RG$genes[,idcol])
      if (is.element(FALSE,valid))
      {
         stop("One or more specified cloneid is not a member of cloneids on the array.")
      }
      allcloneid<-RG[RG$genes[,idcol]%in%cloneid,]$genes[,idcol]
      cloneid<-allcloneid
   }     

   ### check whether cloneid contains NA

   if (NA %in% cloneid)
   {
      stop("Clone ID should not contain NA.")
   }
   if ("" %in% cloneid)
   {
      stop("Clone ID should not be empty.")
   } 
   larray <- length(array)
   lcond <- length(condition)
   ldye <- length(dye)
   if ((larray != lcond) || (larray != ldye) || (lcond != ldye))
   {
      stop("array, condition and dye should have the same length")      
   }
   narray <- ncol(RG$R)
   if (larray != narray * 2)
   {
      stop("design vector and array number don't match")
   }   
   
   bc <- parameter$Method$BC
   area <- parameter$Method$Area
   errormodel <- parameter$ErrorModel

   ### pick up Rinten and Ginten of the specified cloneid
   
   matchclone<-RG$genes[,idcol] %in% cloneid   

   R<-RG$R[matchclone,]
   G<-RG$G[matchclone,]
   Rb<-RG$Rb[matchclone,]
   Gb<-RG$Gb[matchclone,]
   RArea<-RG$RArea[matchclone,]
   GArea<-RG$GArea[matchclone,]

   ### calculate the R/G intensity according to the parameter "bc" and "area"

   if (bc)
   {
      if (area)
      {
         Rinten <- (R - Rb) * RArea
         Ginten <- (G - Gb) * GArea
      }else{
         Rinten <- R - Rb
         Ginten <- G - Gb 
      }
   }else{
      if (area)
      {
         Rinten <- R * RArea
         Ginten <- G * GArea
      }else{
         Rinten <- R
         Ginten <- G
      }
   }

   if (arraynum < 2)
   {
      Rinten <- as.matrix(Rinten)
      Ginten <- as.matrix(Ginten)
   }
   Rinten <- pmax(Rinten, 0.5, rm.na = TRUE)
   Ginten <- pmax(Ginten, 0.5, rm.na = TRUE)

   if (length(cloneid) < 2)
   {
      dimension <- c(1,arraynum)
   }else{
      dimension<-dim(Rinten)
   }

   storage.mode(Rinten)<-storage.mode(Ginten)<-"double"
   storage.mode(dimension)<-"integer"
   storage.mode(array)<-storage.mode(condition)<-storage.mode(dye)<-"integer"
   storage.mode(cloneid)<-storage.mode(errormodel)<-"character"
   if ("AdjustFactor" %in% names(parameter))
   {
      p2 <- parameter$P2 * parameter$AdjustFactor
   }else{
      p2 <- parameter$P2
   }
   tmppar<-list(MuS=parameter$MuS,Ka=parameter$Ka,P1=parameter$P1,P2=p2,SigmaA=parameter$SigmaA,SigmaM=parameter$SigmaM,SigmaS=parameter$SigmaS)
   par<-new("ParameterList",tmppar)
   storage.mode(par)<-"list"

   normalized<-.Call("normalization",Rinten,Ginten,dimension,cloneid,array,condition,dye,par,errormodel,PACKAGE="CALIB")
   
}

plotSpikeCI<-function(spike,parameter,array=1,bc=FALSE,area=TRUE,meanpoint=TRUE,xlab="log(Concentration)",ylab="log(Intensity)",main=colnames(spike$R)[array],
             onlycalib=TRUE,xlim=NULL,ylim=NULL,pch=19,cex=0.2,meanpch=21,meancex=1,lwd=1.5,cy5col="red",cy3col="green",...)
{
   if(missing(spike))
   {
      stop("Must specify spike.")
   }
   if (is.null(spike))
   {
      stop("Must specify spike.")
   }
   if (!is(spike,"SpikeList"))
   {
      stop("Object of SpikeList class is required.")
   }
   if (!missing(parameter))
   {
      if (!is(parameter,"ParameterList"))
      {
         stop("parameter must be object of ParameterList.")
      }
   }  
   ### pick out the right type of spike.   
   if (onlycalib)
   {
      index <- spike$genes$Status %in% "Calibration"
   }else{
      ### get rid of negative spike
      negindex <- spike$genes$Status %in% "Negative"
      index <- !negindex
   }

   if (!missing(parameter))
   {
      bcofpar = parameter$Method$BC
      areaofpar = parameter$Method$Area
      if (bcofpar != bc)
      {
         stop("bc must be the same as parameter estimation.")
      }
      if (areaofpar != area)
      {
         stop("area must be the same as parameter estimation.")
      }
   }
   spike$R <- as.matrix(spike$R)
   spike$G <- as.matrix(spike$G)
   spike$Rb <- as.matrix(spike$Rb)
   spike$Gb <- as.matrix(spike$Gb)
   spike$RArea <- as.matrix(spike$RArea)
   spike$GArea <- as.matrix(spike$GArea)
   spike$RConc <- as.matrix(spike$RConc)
   spike$GConc <- as.matrix(spike$GConc)   
  
   if (bc)
   {
      if (area)
      {
         SRinten <- (spike$R[index,array] - spike$Rb[index,array]) * spike$RArea[otherindex,array]
         SGinten <- (spike$G[index,array] - spike$Gb[index,array]) * spike$GArea[otherindex,array]
      }else{
         SRinten <- spike$R[index,array] - spike$Rb[index,array]
         SGinten <- spike$G[index,array] - spike$Gb[index,array]
      }
   }else{
      if (area)
      {
         SRinten <- spike$R[index,array] * spike$RArea[index,array]
         SGinten <- spike$G[index,array] * spike$GArea[index,array]
      }else{
         SRinten <- spike$R[index,array]
         SGinten <- spike$G[index,array]
      }
   }
   SRconc <- spike$RConc[index,array]
   SGconc <- spike$GConc[index,array]
   
   xred <- log10(SRconc)
   yred <- log10(SRinten)
   xgreen <- log10(SGconc)
   ygreen <- log10(SGinten)
         
   ### Sinten and Sconc are for the range of the plot
   x <- c(xred,xgreen)
   y <- c(yred,ygreen)

   if (is.null(xlim))
   {
      max_x <- max(x)
      min_x <- min(x)
      
      if (is.infinite(max_x))
      {
         index <- x %in% max_x
         index <- !index
         max_x <- max(x[index])
      }
      if (is.infinite(min_x))
      {
         index <- x %in% min_x
         index <- !index
         min_x <- min(x[index])
      }
      xlim <- c(min_x,max_x)
   }
   if (is.null(ylim))
   {
      max_y <- max(y)
      min_y <- min(y)

      if (is.infinite(max_y))
      {
         index <- y %in% max_y
         index <- !index
         max_y <- max(y[index])
      }
      if (is.infinite(min_y))
      {
         index <- y %in% min_y
         index <- !index
         min_y <- min(y[index])
      }
      ylim <- c(min_y,max_y)
   }

   plot(xred,yred,xlab=xlab,ylab=ylab,main=main,xlim=xlim,ylim=ylim,col=cy5col,pch=pch,cex=cex,...)
   points(xgreen,ygreen,col=cy3col,cex=cex,pch=pch)

   ### add mean points if meanpoint=TRUE
   if (meanpoint)
   {
      unique_conc <- unique(spike$RConc[index,array])
      len_unique <- length(unique_conc)
      for (i in 1:len_unique)
      {
         redindex <- spike$RConc[index,array] %in% unique_conc[i]
         redmean <- median(log10(SRinten[redindex]))
         greenindex <- spike$GConc[index,array] %in% unique_conc[i]
         greenmean <- median(log10(SGinten[greenindex]))
         points(log10(unique_conc[i]),redmean,cex=meancex,pch=meanpch,col=cy5col)
         points(log10(unique_conc[i]),greenmean,cex=meancex,pch=meanpch,col=cy3col)
      }       
   }
   if (!missing(parameter))
   {
      ka <- parameter$Ka[array]
      mus <- parameter$MuS[array]
      p1R <- parameter$P1[1,array]
      p2R <- parameter$P2[1,array]
      p1G <- parameter$P1[2,array]
      p2G <- parameter$P2[2,array]
 
      x <- xred
      curve(log10((p1R * 10^x * mus)/(1/ka + 10^x + 10^x ) + p2R),xlim[1],xlim[2],add=TRUE,col=cy5col,lwd=lwd)

      x<- xgreen
      curve(log10((p1G * 10^x * mus)/(1/ka + 10^x + 10^x ) + p2G),xlim[1],xlim[2],add=TRUE,col=cy3col,lwd=lwd)
   }   
}

plotSpikeRG<-function(spike,parameter,RG,array=1,bc=FALSE,area=TRUE,xlab="log(Rintensity)",ylab="log(Gintensity)",main=colnames(spike$R)[array],
              onlycalib=FALSE,xlim=NULL,ylim=NULL,pch=19,cex=0.3,col="black",allpch=19,allcex=0.05,allcol="lightgrey",diag=TRUE,diagcol="grey",diaglwd=1,
              curvecol="blue",curvelwd=1.5,calibtype=1,adjusttype=4,...)
{
   if (missing(spike))
   {
      stop("Must specify spike.")
   }
   if (is.null(spike))
   {
      stop("Must specify spike.")
   }
   if (!is(spike,"SpikeList"))
   {
      stop("Object of SpikeList class is required.")
   }
   if (!missing(parameter))
   {
      if (!is(parameter,"ParameterList"))
      {
         stop("parameter must be object of ParameterList.")
      }
      if ("AdjustFactor" %in% names(parameter))
      {
         if (missing(RG))
         {
            stop("Must specify RG data.")
         }
         if (is.null(RG))
         {
            stop("Must specify RG data.")
         }
         if (!is(RG,"RGList_CALIB"))
         {
            stop("RG must be object of RGList_CALIB.")
         }
      }   
   }
   
   ### pick out the right type of spike.   
   if (onlycalib)
   {
      index <- spike$genes$Status %in% "Calibration"
   }else{
      ### get rid of negative spike
      negindex <- spike$genes$Status %in% "Negative"
      index <- !negindex
   }
   
   ### pick out all the data except spike in RG
   if (!missing(RG))
   {
      index1 <- RG$genes$Status %in% "Calibration"
      index2 <- RG$genes$Status %in% "Ratio"
      index3 <- RG$genes$Status %in% "Negative"
      indexRG <- !(index1 | index2 | index3)

      RG$R <- as.matrix(RG$R)
      RG$G <- as.matrix(RG$G)
      RG$Rb <- as.matrix(RG$Rb)
      RG$Gb <- as.matrix(RG$Gb)
      RG$RArea <- as.matrix(RG$RArea)
      RG$GArea <- as.matrix(RG$GArea)
   } 

   if (!missing(parameter))
   {
      bcofpar = parameter$Method$BC
      areaofpar = parameter$Method$Area
      if (bcofpar != bc)
      {
         stop("bc must be the same as parameter estimation.")
      }
      if (areaofpar != area)
      {
         stop("area must be the same as parameter estimation.")
      }
   }

   spike$R <- as.matrix(spike$R)
   spike$G <- as.matrix(spike$G)
   spike$Rb <- as.matrix(spike$Rb)
   spike$Gb <- as.matrix(spike$Gb)
   spike$RArea <- as.matrix(spike$RArea)
   spike$GArea <- as.matrix(spike$GArea)
   spike$RConc <- as.matrix(spike$RConc)
   spike$GConc <- as.matrix(spike$GConc)   
      
   if (bc)
   {
      if (area)
      {
         SRinten <- (spike$R[index,array] - spike$Rb[index,array]) * spike$RArea[index,array]
         SGinten <- (spike$G[index,array] - spike$Gb[index,array]) * spike$GArea[index,array]
         if (!missing(parameter))
         {
            Rinten <- (RG$R[indexRG,array] - RG$Rb[indexRG,array]) * RG$RArea[indexRG,array]
            Ginten <- (RG$G[indexRG,array] - RG$Gb[indexRG,array]) * RG$GArea[indexRG,array]
         }
      }else{
         SRinten <- spike$R[index,array] - spike$Rb[index,array]
         SGinten <- spike$G[index,array] - spike$Gb[index,array]
         if (!missing(parameter))
         {
            Rinten <- RG$R[indexRG,array] - RG$Rb[indexRG,array]
            Ginten <- RG$G[indexRG,array] - RG$Gb[indexRG,array]
         }
      }
   }else{
      if (area)
      {
         SRinten <- spike$R[index,array] * spike$RArea[index,array]
         SGinten <- spike$G[index,array] * spike$GArea[index,array]
         if (!missing(parameter))
         {
            Rinten <- RG$R[indexRG,array] * RG$RArea[indexRG,array]
            Ginten <- RG$G[indexRG,array] * RG$GArea[indexRG,array]
         }
      }else{
         SRinten <- spike$R[index,array]
         SGinten <- spike$G[index,array]
         if (!missing(parameter))
         {
            Rinten <- RG$R[indexRG,array]
            Ginten <- RG$G[indexRG,array]
         }
      }
   }
 
   x <- log10(SRinten)
   y <- log10(SGinten)

   xy<-c(x,y)
         
   if (is.null(xlim))
   {
      max_xy <- max(xy)
      min_xy <- min(xy)
      
      if (is.infinite(max_xy))
      {
         index <- xy %in% max_xy
         index <- !index
         max_xy <- max(xy[index])
      }
      if (is.infinite(min_xy))
      {
         index <- xy %in% min_xy
         index <- !index
         min_xy <- min(xy[index])
      }
      xlim <- c(min_xy,max_xy)
   }
   if (is.null(ylim))
   {
      max_xy <- max(xy)
      min_xy <- min(xy)

      if (is.infinite(max_xy))
      {
         index <- xy %in% max_xy
         index <- !index
         max_xy <- max(xy[index])
      }
      if (is.infinite(min_xy))
      {
         index <- xy %in% min_xy
         index <- !index
         min_xy <- min(xy[index])
      }
      ylim <- c(min_xy,max_xy)
   }

   if (!missing(parameter))
   {
      plot(log10(Rinten),log10(Ginten),xlim=xlim,ylim=ylim,xlab=xlab,ylab=ylab,main=main,pch=allpch,cex=allcex,col=allcol,...)
      points(x,y,pch=pch,cex=cex,col=col,...)
   }else{
      plot(x,y,xlim=xlim,ylim=ylim,xlab=xlab,ylab=ylab,main=main,pch=pch,cex=cex,col=col,...)
   }
   
   if (diag)
   {
      abline(0,1,col=diagcol,lwd=diaglwd)
   }
   
   if (!missing(parameter))
   {
      spikenoneg <- spike[spike$genes$Status!="Negative",array] 
      ratio <- spikenoneg$RConc[,array]/spikenoneg$GConc[,array]
      uniqueratio <- unique(ratio)

      p1R <- parameter$P1[1,array]
      p1G <- parameter$P1[2,array]
      p2R <- parameter$P2[1,array]
      p2G <- parameter$P2[2,array]

      if ("AdjustFactor" %in% names(parameter))
      {
         p2Ra <- parameter$P2[1,array] * parameter$AdjustFactor[1,array]
         p2Ga <- parameter$P2[2,array] * parameter$AdjustFactor[2,array]
      }
      
      ### calculate the lower boundary of x
      if (p2R > p2G)
      {
         xlower <- log10(p2G)
      }else{
         xlower <- log10(p2R)
      }
      xupper <- xlim[2] 
      
      for (i in 1:length(uniqueratio))
      {
         index <- ratio %in% uniqueratio[i]
         x <- log10(SRinten[index])
         if (uniqueratio[i] == 1)
         {
            curve(log10((10^x-p2R)* p1G * uniqueratio[i]/p1R + p2G),xlower,xupper,add=TRUE,col=curvecol,lwd=curvelwd,lty=calibtype) 
            if ("AdjustFactor" %in% names(parameter))   
            {
               curve(log10((10^x-p2Ra)* p1G * uniqueratio[i]/p1R + p2Ga),xlower,xupper,add=TRUE,col=curvecol,lwd=curvelwd,lty=adjusttype)
            }    
         }
      }      
   }
}

plotSpikeHI<-function(spike,parameter,array=1,xlab="log(Hybridized)",ylab="log(Intensity)",main=colnames(spike$R)[array],
             xlim=NULL,ylim=NULL,pch=19,cex=0.2,cy5col="black",cy3col="black",noerror=TRUE,noepch=19,noecex=0.1,
             noecy5col="lightpink",noecy3col="lightblue",curve=TRUE,lwd=1.5,curvecy5col="red",curvecy3col="green",...)
{
   ### first check three parameters: spike, parameter and idcol
 
   if(missing(spike))
   {
      stop("Must specify spike.")
   }
   if (is.null(spike))
   {
      stop("Must specify spike.")
   }
   if (!is(spike,"SpikeList"))
   {
      stop("Object of SpikeList class is required.")
   }
   if(missing(parameter))
   {
      stop("Must specify parameter.")
   }
   if (is.null(parameter))
   {
      stop("Must specify parameter.")
   }
   if (!is(parameter,"ParameterList"))
   {
      stop("Object of ParameterList class is required.")
   }
   
   ### calculate Xs for each spot.
   bc <- parameter$Method$BC
   area <- parameter$Method$Area
   errormodel <- parameter$ErrorModel
      
   index <- spike$genes$Status %in% "Calibration"

   spike$R <- as.matrix(spike$R)
   spike$G <- as.matrix(spike$G)
   spike$Rb <- as.matrix(spike$Rb)
   spike$Gb <- as.matrix(spike$Gb)
   spike$RArea <- as.matrix(spike$RArea)
   spike$GArea <- as.matrix(spike$GArea)
   spike$RConc <- as.matrix(spike$RConc)
   spike$GConc <- as.matrix(spike$GConc)   

   if (bc)
   {
      if (area)
      {
         SRinten <- (spike$R[index,array] - spike$Rb[index,array]) * spike$RArea[index,array]
         SGinten <- (spike$G[index,array] - spike$Gb[index,array]) * spike$GArea[index,array]
      }else{
         SRinten <- spike$R[index,array] - spike$Rb[index,array]
         SGinten <- spike$G[index,array] - spike$Gb[index,array]
      }
   }else{
      if (area)
      {
         SRinten <- spike$R[index,array] * spike$RArea[index,array]
         SGinten <- spike$G[index,array] * spike$GArea[index,array]
      }else{
         SRinten <- spike$R[index,array]
         SGinten <- spike$G[index,array]
      }
   }

   SRconc <- spike$RConc[index,array]
   SGconc <- spike$GConc[index,array]

   #cloneid_spike <- spike$genes[index,idcol]
   
   index_par <- parameter$genes$Status %in% "Calibration"
  
   spoterror <- parameter$SpotError[index_par,array]
  
   mus <- parameter$MuS[,array]
   ka <- parameter$Ka[,array]
     
   if (errormodel == "A")
   {
      SRXs <- SRconc * (mus + spoterror) / ((1/ka) + SRconc + SGconc)
      SGXs <- SGconc * (mus + spoterror) / ((1/ka) + SRconc + SGconc)
   }else{
      if (errormodel == "M")
      {
         SRXs <- SRconc * mus * exp(spoterror) / ((1/ka) + SRconc + SGconc)
         SGXs <- SGconc * mus * exp(spoterror) / ((1/ka) + SRconc + SGconc)
      }else{
         stop("errormodel is not A or M.")
      }
   }
   
   xred <- log10(SRXs)
   yred <- log10(SRinten)
   xgreen <- log10(SGXs)
   ygreen <- log10(SGinten)

   x <- c(xred,xgreen)
   y <- c(yred,ygreen)

   if (is.null(xlim))
   {
      max_x <- max(x)
      min_x <- min(x)
      
      if (is.infinite(max_x))
      {
         index <- x %in% max_x
         index <- !index
         max_x <- max(x[index])
      }
      if (is.infinite(min_x))
      {
         index <- x %in% min_x
         index <- !index
         min_x <- min(x[index])
      }
	xlim <- c(0,max_x)
   }
   if (is.null(ylim))
   {
      max_y <- max(y)
      min_y <- min(y)

      if (is.infinite(max_y))
      {
         index <- y %in% max_y
         index <- !index
         max_y <- max(y[index])
      }
      if (is.infinite(min_y))
      {
         index <- y %in% min_y
         index <- !index
         min_y <- min(y[index])
      }
      ylim <- c(min_y,max_y)
   }

   plot(xred,yred,xlab=xlab,ylab=ylab,main=main,pch=pch,cex=cex,col=cy5col,xlim=xlim,ylim=ylim,...)
   points(xgreen,ygreen,col=cy3col,cex=cex,pch=pch)

   ### add points without spoterror
   if (noerror)
   {
      SRXs_noerror <- (SRconc * mus)  / ((1/ka) + SRconc + SGconc)
      SGXs_noerror <- (SGconc * mus) / ((1/ka) + SRconc + SGconc)
      
      xred_noerror <- log10(SRXs_noerror)
      xgreen_noerror <- log10(SGXs_noerror)

      points(xred_noerror,yred,col=noecy5col,cex=noecex,pch=noepch)
      points(xgreen_noerror,ygreen,col=noecy3col,cex=noecex,pch=noepch)
   }

   ### add parameter curve
   if (curve)
   {
      x <- xred
      p1R <- parameter$P1[1,array]
      p2R <- parameter$P2[1,array]
      curve(log10(p1R * 10^x + p2R),xlim[1],xlim[2],add=TRUE,col=curvecy5col,lwd=lwd)

      x <- xgreen
      p1G <- parameter$P1[2,array]
      p2G <- parameter$P2[2,array]
      curve(log10(p1G * 10^x + p2G),xlim[1],xlim[2],add=TRUE,col=curvecy3col,lwd=lwd)
   }   
}

plotNormalizedData<-function(data,condition=c(1,2),xlab=NULL,ylab=NULL,main=NULL,xlim=NULL,ylim=NULL,
                    pch=19,cex=0.2,col="black",diag=TRUE,diagcol="blue",diaglwd=1.5,...)
{
   if (missing(data))
   {
      stop("Must specify normalized data.")
   }
   condition <- as.character(condition)
   if (length(condition)!= 2)
   {
      stop("Only two conditions can be plotted at a time.")
   }
   if (!is.element(condition[1],colnames(data)))
   {
      stop("User specified conditions are not in data columns.")
   }
   if (!is.element(condition[2],colnames(data)))
   {
      stop("User specified conditions are not in data columns.")
   }
   if (is.null(xlab))
   {
      xlab <- paste("condition",condition[1])
   }
   if (is.null(ylab))
   {
      ylab <- paste("condition",condition[2])
   }
   if (is.null(main))
   {
      main <- paste("condition",condition[1],"vs","condition",condition[2])
   }
   
   x <- log10(data[,condition[1]])
   y <- log10(data[,condition[2]])
  
   ### calculate xlim and ylim taking into account there exists "infinite"

   if (is.null(xlim))
   {
      max_x <- max(x)
      min_x <- min(x)
      
      if (is.infinite(max_x))
      {
         index <- x %in% max_x
         index <- !index
         max_x <- max(x[index])
      }
      if (is.infinite(min_x))
      {
         index <- x %in% min_x
         index <- !index
         min_x <- min(x[index])
      }
      #xlim <- c(min_x,max_x)
      xlim <- c(0,max_x)
   }
   if (is.null(ylim))
   {
      max_y <- max(y)
      min_y <- min(y)

      if (is.infinite(max_y))
      {
         index <- y %in% max_y
         index <- !index
         max_y <- max(y[index])
      }
      if (is.infinite(min_y))
      {
         index <- y %in% min_y
         index <- !index
         min_y <- min(y[index])
      }
      #ylim <- c(min_y,max_y)
      ylim <- c(0,max_y)
   }

   plot(x,y,xlab=xlab,ylab=ylab,main=main,xlim=xlim,ylim=ylim,pch=pch,cex=cex,col=col,...)

   if (diag)
   {
      abline(0,1,col=diagcol,lwd=diaglwd)
   }
}

plotSpikeSpotError<-function(parameter,array=1,plottype="hist",width=1,plotnames=NULL,main=NULL,...)
{
   if(missing(parameter))
   {
      stop("Must specify parameter.")
   }
   if (is.null(parameter))
   {
      stop("Must specify parameter.")
   }
   if (!is(parameter,"ParameterList"))
   {
      stop("Object of ParameterList class is required.")
   }
   allplottype <- c("hist","boxplot","dens")
   if (!is.element(plottype,allplottype))
   {
       stop("User specified plot type is not variable.")
   }
   switch(plottype,
      hist = {
                spoterror <- parameter$SpotError[,array]
                if (is.null(main))
                {
                   main <- paste("Spike spot error Histogram of",colnames(parameter$SpotError)[array])
                }
                hist(spoterror,main=main,...)
             },
      boxplot = {
               spoterror <- parameter$SpotError
               if (is.null(main))
               {
                  main <- "Spike spot error Boxplot"
               } 
               if (is.null(plotnames))
               {
                  plotnames <- colnames(parameter$SpotError)
               }
               boxplot(data.frame(spoterror),names=plotnames,main=main,...)
            },
      dens = {
                spoterror <- parameter$SpotError[,array]
                if (is.null(main))
                {
                   main <- paste("Spike spot error Density plot of",colnames(parameter$SpotError)[array])
                }
                if (is.null(width))
                {
                   width <- 1
                }
                plot(density(spoterror,width=width),main=main,...)
             }
)
}

setMethod("show","RGList_CALIB",function(object){
   cnames <- names(object)
   for (a in cnames) 
   {
      cat("$",a,"\n",sep="")
      printHead(object[[a]])
      cat("\n")
   }
})

setMethod("show","SpikeList",function(object){
   cnames <- names(object)
   for (a in cnames) 
   {
      cat("$",a,"\n",sep="")
      printHead(object[[a]])
      cat("\n")
   }
})

setMethod("show","ParameterList",function(object){
   cnames <- names(object)
   for (a in cnames) 
   {
      cat("$",a,"\n",sep="")
      printHead(object[[a]])
      cat("\n")
   }
})

assign("[.RGList_CALIB",function(object,i,j,...)
{
   if (nargs()!= 3)
   {
      stop("Two subscripts required", call.= FALSE)
   }
   other <- names(object$other)
   if (missing(i))
   {
      if (missing(j))
      {
         return(object)
      }else{
         object$R <- object$R[,j,drop=FALSE]
	   object$G <- object$G[,j,drop=FALSE]
	   object$Rb <- object$Rb[,j,drop=FALSE]
	   object$Gb <- object$Gb[,j,drop=FALSE]
         object$RArea <- object$RArea[,j,drop=FALSE]
         object$GArea <- object$GArea[,j,drop=FALSE]
	   object$weights <- object$weights[,j,drop=FALSE]
	   object$targets <- object$targets[j,,drop=FALSE]
	   for(k in other) 
         {
            object$other[[k]] <- object$other[[k]][,j,drop=FALSE]
         }
	}
   }else{
      if(is.character(i))
      {
	   i <- match(i,rownames(object))
	   i <- i[!is.na(i)]
	}
	if(missing(j)) {
	   object$R <- object$R[i,,drop=FALSE]
	   object$G <- object$G[i,,drop=FALSE]
	   object$Rb <- object$Rb[i,,drop=FALSE]
	   object$Gb <- object$Gb[i,,drop=FALSE]
         object$RArea <- object$RArea[i,,drop=FALSE]
         object$GArea <- object$GArea[i,,drop=FALSE]
	   object$weights <- object$weights[i,,drop=FALSE]
	   object$genes <- object$genes[i,,drop=FALSE]
	   for(k in other)
         {
            object$other[[k]] <- object$other[[k]][i,,drop=FALSE]
         }
	} else {
	   object$R <- object$R[i,j,drop=FALSE]
	   object$G <- object$G[i,j,drop=FALSE]
	   object$Rb <- object$Rb[i,j,drop=FALSE]
	   object$Gb <- object$Gb[i,j,drop=FALSE]
         object$RArea <- object$RArea[i,j,drop=FALSE]
         object$GArea <- object$GArea[i,j,drop=FALSE]
	   object$weights <- object$weights[i,j,drop=FALSE]
	   object$genes <- object$genes[i,,drop=FALSE]
	   object$targets <- object$targets[j,,drop=FALSE]
	   for(k in other)
         {
            object$other[[k]] <- object$other[[k]][i,j,drop=FALSE]
         }
      }
   }
   return(object)        
})

assign("[.SpikeList",function(object,i,j,...)
{
   if (nargs()!= 3)
   {
      stop("Two subscripts required", call.= FALSE)
   }
   other <- names(object$other)
   if (missing(i))
   {
      if (missing(j))
      {
         return(object)
      }else{
         object$R <- object$R[,j,drop=FALSE]
	   object$G <- object$G[,j,drop=FALSE]
	   object$Rb <- object$Rb[,j,drop=FALSE]
	   object$Gb <- object$Gb[,j,drop=FALSE]
         object$RArea <- object$RArea[,j,drop=FALSE]
         object$GArea <- object$GArea[,j,drop=FALSE]
         object$RConc <- object$RConc[,j,drop=FALSE]
         object$GConc <- object$GConc[,j,drop=FALSE]
	   object$weights <- object$weights[,j,drop=FALSE]
	   object$targets <- object$targets[j,,drop=FALSE]
	   for(k in other) 
         {
            object$other[[k]] <- object$other[[k]][,j,drop=FALSE]
         }
	}
   }else{
      if(is.character(i))
      {
	   i <- match(i,rownames(object))
	   i <- i[!is.na(i)]
	}
	if(missing(j)) {
	   object$R <- object$R[i,,drop=FALSE]
	   object$G <- object$G[i,,drop=FALSE]
	   object$Rb <- object$Rb[i,,drop=FALSE]
	   object$Gb <- object$Gb[i,,drop=FALSE]
         object$RArea <- object$RArea[i,,drop=FALSE]
         object$GArea <- object$GArea[i,,drop=FALSE]
         object$RConc <- object$RConc[i,,drop=FALSE]
         object$GConc <- object$GConc[i,,drop=FALSE]
	   object$weights <- object$weights[i,,drop=FALSE]
	   object$genes <- object$genes[i,,drop=FALSE]
	   for(k in other)
         {
            object$other[[k]] <- object$other[[k]][i,,drop=FALSE]
         }
	} else {
	   object$R <- object$R[i,j,drop=FALSE]
	   object$G <- object$G[i,j,drop=FALSE]
	   object$Rb <- object$Rb[i,j,drop=FALSE]
	   object$Gb <- object$Gb[i,j,drop=FALSE]
         object$RArea <- object$RArea[i,j,drop=FALSE]
         object$GArea <- object$GArea[i,j,drop=FALSE]
         object$RConc <- object$RConc[i,j,drop=FALSE]
         object$GConc <- object$GConc[i,j,drop=FALSE]
	   object$weights <- object$weights[i,j,drop=FALSE]
	   object$genes <- object$genes[i,,drop=FALSE]
	   object$targets <- object$targets[j,,drop=FALSE]
	   for(k in other)
         {
            object$other[[k]] <- object$other[[k]][i,j,drop=FALSE]
         }
      }
   }
   return(object)   
})

assign("[.ParameterList",function(object,i,j,...)
{
   if (nargs()!= 3)
   {
      stop("Two subscripts required", call.= FALSE)
   }
   if (missing(i))
   {
      if (missing(j))
      {
         return(object)
      }else{
         object$MuS <- object$MuS[,j,drop=FALSE]
	   object$Ka <- object$Ka[,j,drop=FALSE]
	   object$P1 <- object$P1[,j,drop=FALSE]
	   object$P2 <- object$P2[,j,drop=FALSE]
         object$SigmaA <- object$SigmaA[,j,drop=FALSE]
         object$SigmaM <- object$SigmaM[,j,drop=FALSE]
         object$SigmaS <- object$SigmaS[,j,drop=FALSE]
         object$SpotError <- object$SpotError[,j,drop=FALSE]
	}
   }else{
      stop("It makes no sense to exact rows from ParameterList object.")
   }
   return(object)  
})

cbind.RGList_CALIB<-function (..., deparse.level = 1) 
{
   objects <- list(...)
   nobjects <- length(objects)
   out <- objects[[1]]
   other <- names(objects[[1]]$other)
   if (nobjects > 1) 
       for (i in 2:nobjects) {
           out$R <- cbind(out$R, objects[[i]]$R)
           out$G <- cbind(out$G, objects[[i]]$G)
           out$Rb <- cbind(out$Rb, objects[[i]]$Rb)
           out$Gb <- cbind(out$Gb, objects[[i]]$Gb)
           out$RArea <- cbind(out$RArea,objects[[i]]$RArea)
           out$GArea <- cbind(out$GArea,objects[[i]]$GArea)
           out$weights <- cbind(out$weights, objects[[i]]$weights)
           out$targets <- rbind(out$targets, objects[[i]]$targets)
           for (a in other) out$other[[a]] <- cbind(out$other[[a]], 
               objects[[i]]$other[[a]])
       }
   return(out)
}

cbind.SpikeList<-function(..., deparse.level = 1)
{
   objects <- list(...)
   nobjects <- length(objects)
   out <- objects[[1]]
   if (nobjects > 1) 
       for (i in 2:nobjects) {
           out$R <- cbind(out$R, objects[[i]]$R)
           out$G <- cbind(out$G, objects[[i]]$G)
           out$Rb <- cbind(out$Rb, objects[[i]]$Rb)
           out$Gb <- cbind(out$Gb, objects[[i]]$Gb)
           out$RArea <- cbind(out$RArea,objects[[i]]$RArea)
           out$GArea <- cbind(out$GArea,objects[[i]]$GArea)
           out$RConc <- cbind(out$RConc,objects[[i]]$RConc)
           out$GConc <- cbind(out$GConc,objects[[i]]$GConc)
       }
   return(out)
}

cbind.ParameterList<-function(..., deparse.level = 1)
{
   objects <- list(...)
   nobjects <- length(objects)
   out <- objects[[1]]
   if (nobjects > 1) 
       for (i in 2:nobjects) {
           out$MuS <- cbind(out$MuS, objects[[i]]$MuS)
           out$Ka <- cbind(out$Ka, objects[[i]]$Ka)
           out$P1 <- cbind(out$P1, objects[[i]]$P1)
           out$P2 <- cbind(out$P2, objects[[i]]$P2)
           out$SigmaA <- cbind(out$SigmaA,objects[[i]]$SigmaA)
           out$SigmaM <- cbind(out$SigmaM,objects[[i]]$SigmaM)
           out$SigmaS <- cbind(out$SigmaS,objects[[i]]$SigmaS)
           out$SpotError <- cbind(out$SpotError,objects[[i]]$SpotError)
           out$Adjust <- cbind(out$Adjust,objects[[i]]$Adjust)
       }
   return(out)
}

rbind.RGList_CALIB<-function (..., deparse.level = 1) 
{
   objects <- list(...)
   nobjects <- length(objects)
   out <- objects[[1]]
   other <- names(objects[[1]]$other)
   if (nobjects > 1) 
       for (i in 2:nobjects) {
           out$R <- rbind(out$R, objects[[i]]$R)
           out$G <- rbind(out$G, objects[[i]]$G)
           out$Rb <- rbind(out$Rb, objects[[i]]$Rb)
           out$Gb <- rbind(out$Gb, objects[[i]]$Gb)
           out$RArea <- rbind(out$RArea,objects[[i]]$RArea)
           out$GArea <- rbind(out$GArea,objects[[i]]$GArea)
           out$weights <- rbind(out$weights, objects[[i]]$weights)
           out$genes <- rbind(out$genes, objects[[i]]$genes)
           for (a in other) out$other[[a]] <- rbind(out$other[[a]], 
               objects[[i]]$other[[a]])
       }
   return(out)
}

rbind.SpikeList<-function (..., deparse.level = 1) 
{
   objects <- list(...)
   nobjects <- length(objects)
   out <- objects[[1]]
   if (nobjects > 1) 
       for (i in 2:nobjects) {
           out$R <- rbind(out$R, objects[[i]]$R)
           out$G <- rbind(out$G, objects[[i]]$G)
           out$Rb <- rbind(out$Rb, objects[[i]]$Rb)
           out$Gb <- rbind(out$Gb, objects[[i]]$Gb)
           out$RArea <- rbind(out$RArea,objects[[i]]$RArea)
           out$GArea <- rbind(out$GArea,objects[[i]]$GArea)
           out$RConc <- rbind(out$RArea,objects[[i]]$RConc)
           out$GConc <- rbind(out$GArea,objects[[i]]$GConc)
           out$genes <- rbind(out$genes, objects[[i]]$genes)
       }
   return(out)
}

dim.RGList_CALIB<-function(x)
{
   if(is.null(x$R))
   {
      c(0,0)
   }else{
      dim(as.matrix(x$R))
   }
}

dim.SpikeList<-function(x)
{
   if (is.null(x$R))
   {
      c(0,0)
   }else{
      dim(as.matrix(x$R))
   }
}

length.RGList_CALIB<-function(x)
{
  len<-prod(dim(x))
  return(len)
}

length.SpikeList<-function(x)
{
   len<-prod(dim(x))
   return(len)
}

dimnames.RGList_CALIB<-function(x)
{
   dimnames(x$R)
}

dimnames.SpikeList<-function(x)
{
   dimnames(x$R)
}

merge.RGList_CALIB<-function(x,y,...)
{
    if (!is(y, "RGList_CALIB")) 
        stop("both x and y must be RGList_CALIB objects")
    genes1 <- rownames(x$R)
    if (is.null(genes1)) 
        genes1 <- rownames(x$G)
    if (is.null(genes1)) 
        genes1 <- x$genes$ID
    genes2 <- rownames(y$R)
    if (is.null(genes2)) 
        genes2 <- rownames(y$G)
    if (is.null(genes2)) 
        genes2 <- y$genes$ID
    if (is.null(genes1) || is.null(genes2)) 
        stop("Need row names to align on")
    fields1 <- names(x)
    fields2 <- names(y)
    if (!identical(fields1, fields2)) 
        stop("The two RGLists_CALIB have different components")
    ord2 <- match(makeUnique(genes1), makeUnique(genes2))
    cbind(x, y[ord2, ])
}

merge.SpikeList<-function (x, y, ...) 
{
   if (!(is(y, "SpikeList") || is(x,"SpikeList"))) 
       stop("both x and y must be SpikeList objects")
   genes1 <- rownames(x$R)
   if (is.null(genes1)) 
       genes1 <- rownames(x$G)
   if (is.null(genes1)) 
       genes1 <- x$genes$ID
   genes2 <- rownames(y$R)
   if (is.null(genes2)) 
       genes2 <- rownames(y$G)
   if (is.null(genes2)) 
       genes2 <- y$genes$ID
   if (is.null(genes1) || is.null(genes2)) 
       stop("Need row names to align on")
   fields1 <- names(x)
   fields2 <- names(y)
   if (!identical(fields1, fields2)) 
       stop("The two SpikeLists have different components")
   ord2 <- match(makeUnique(genes1), makeUnique(genes2))
   cbind(x, y[ord2, ])
}

calibReadMe<-function(view=TRUE)
{
    f <- system.file("doc","readme.pdf",package="CALIB")
    if (view)
    {
        if (.Platform$OS.type == "windows")
        {
            shell.exec(f)
        }else{
            system(paste(Sys.getenv("R_PDFVIEWER"),f,"&")) 
        }
    }
    return(f)
}

getColClasses<-function (cols, ...) 
{
    cols <- as.character(cols)
    x <- rep("NULL", length(cols))
    names(x) <- cols
    wanted <- list(...)
    for (i in 1:length(wanted)) {
        if (is.null(wanted[[i]])) 
            next
        if (is.function(wanted[[i]])) 
            include <- namesInFun(cols, wanted[[i]])
        if (is.list(wanted[[i]])) 
            wanted[[i]] <- unlist(wanted[[i]])
        if (is.character(wanted[[i]])) 
            include <- wanted[[i]]
        ind <- match(include, cols, nomatch = 0)
        x[ind] <- NA
    }
    x
}