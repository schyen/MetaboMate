#' Importing 1D NMR spectra from Bruker
#' @export
#' @description  Reading-in 1D NMR spectra in Bruker file format. The function recursively descends the directory tree of a given location and is searching for corresponding processing status (\emph{procs}), acquisition status (\emph{acqus}) and spectrum (\emph{r1}) files.
#' @param path A character string, specifying the path of the experiment folder location.
#' @param filter Logical - if set to TRUE, then experiments with matching \emph{acqus}, \emph{procs} and \emph{1r} files are read-in, ignored are corrupt file systems and 2D experiments. If set to FALSE, the algorithm will abort if a file system is corrupt or if 2D experiments are present. Generally, the filter=TRUE option should be preferred except for some file system testing operations.
#' @return The following variables are automatically returned as globally defined workspace variables:
#' \describe{
#' \item{X}{two-dimensional NMR matrix with spectra in rows and chemical shift variables in columns.}
#' \item{ppm}{chemical shift vector (ppm).}
#' \item{meta}{spectrometer metadata extracted from acqus and procs files (eg time of spectral acquisition).}
#' }
#' @examples
#' path=system.file("extdata/", package = "MetaboMate")
#' readBruker(path)
#' @importFrom stats approxfun complete.cases
#' @author Torben Kimhofer \email{tkimhofer@@gmail.com}
#' @details The columns in meta represent spectrometer meta-variables and prefix \code{a} and \code{p} indicate if a parameter was extracted from the acquisition status (\emph{acqus}) or processing status (\emph{procs}) file, respectively. Parameter names without prefix were calculated on the fly (eg run-order). Spectrometer variables that often are of interest to s specific study include a_NS (number of scans), a_RG (receiver gain), a_Date (acquisition date), p_SF (spectrometer frequency) and a_AUNM (AU program).



readBruker=function (path, filter = T)
{
  warnDef <- options("warn")$warn
  warnRead <- options(warn = -1)
  if (grepl("/$", path)) {
    path <- gsub("/$", "", path)
  }
  afile <- list.files(path = path, pattern = "^acqus$", all.files = FALSE,
                      full.names = TRUE, recursive = TRUE, ignore.case = TRUE)
  pfile <- list.files(path = path, pattern = "^procs$", all.files = FALSE,
                      full.names = TRUE, recursive = TRUE, ignore.case = TRUE)
  pfile <- grep("pdata/1/", pfile, value = T)
  rfile <- list.files(path = path, pattern = "^1r$", all.files = FALSE,
                      full.names = TRUE, recursive = TRUE, ignore.case = TRUE)
  rfile <- grep("pdata/1/", rfile, value = T)
  Lp <- length(pfile)
  if (Lp == 0 || length(rfile) == 0 || length(afile) == 0) {
    stop("No existing Bruker files in specified path.")
  }
  a.exp <- gsub(paste("^", path, "/|/acqus$", sep = ""), "",
                afile)
  p.exp <- gsub(paste("^", path, "/|/pdata/1/procs$", sep = ""),
                "", pfile)
  r.exp <- gsub(paste("^", path, "/|/pdata/1/1r$", sep = ""),
                "", rfile)
  idx.a = a.exp %in% r.exp & a.exp %in% p.exp
  idx.p = p.exp %in% r.exp & p.exp %in% a.exp
  idx.r = r.exp %in% a.exp & r.exp %in% p.exp
  if (filter == T) {
    afile = afile[idx.a]
    pfile = pfile[idx.p]
    rfile = rfile[idx.r]
  }else {
    not.present <- list()
    not.present[[1]] = afile[!idx.a]
    not.present[[2]] <- pfile[!idx.p]
    not.present[[3]] <- rfile[!idx.r]
    if (length(not.present[[1]]) > 0) {
      cat("File structure depreciated: missing Bruker acquisition file(s). See output for file names.\n")
    }
    if (length(not.present[[2]]) > 0) {
      cat("File structure depreciated: missing Bruker processing file(s). See output for experiment folder names.\n")
    }
    if (length(not.present[[3]]) > 0) {
      cat("File structure depreciated: missing aquisition file(s). See output for experiment folder names.\n")
    }
    if (sum(sapply(not.present, length)) > 0) {
      return(sort(unlist(not.present)))
    }
  }
  Lp <- length(pfile)
  if (Lp == 0) {
    stop("No matching Bruker files in specified path.")
  }
  if (Lp == 1) {
    cat("Reading ", Lp, " spectrum.\n", sep = "")
  }else {
    cat("Reading ", Lp, " spectra.\n", sep = "")
  }
  out <- apply(rbind(1:Lp), 2, function(i, pf = pfile, rf = rfile,
                                        af = afile) {
    con <- file(af[i], open = "r")
    aLine <- readLines(con, n = -1, warn = FALSE)
    close(con)
    aLine <- aLine[-c(1:9)]
    aLine <- aLine[!grepl("##END", aLine)]
    idx <- grep("#", aLine)
    aLineC <- array()
    for (j in 1:length(idx)) {
      if (grepl("#", aLine[idx[j] + 1])) {
        aLineC[j] <- aLine[idx[j]]
        next
      }
      end_idx <- grep("#", aLine[(idx[j] + 1):length(aLine)])[1]
      aLineC[j] <- paste(aLine[idx[j]:(idx[j] + end_idx -
                                         1)], collapse = " ")
    }
    aLineC <- gsub(" \\(.*\\)", "", aLineC)
    myV <- strsplit(aLineC, "=")
    df.acqu <- unique(as.data.frame(do.call(rbind, myV),
                                    stringsAsFactors = F))
    df.acqu$V1 <- paste("a", gsub("##\\$", "", df.acqu$V1),
                        sep = "_")
    df.acqu$V2 <- gsub("^ | $", "", df.acqu$V2)
    con <- file(pf[i], open = "r")
    aLine <- readLines(con, n = -1, warn = FALSE)
    myV <- strsplit(aLine, "=")
    close(con)
    idx <- grepl("##", aLine) & sapply(myV, length) == 2
    myV <- myV[idx]
    df.proc <- unique(as.data.frame(do.call(rbind, myV[1:(length(myV) -
                                                            1)]), stringsAsFactors = F))
    df.proc$V1 <- paste("p", gsub("^##[\\$]?", "", df.proc$V1),
                        sep = "_")
    df.proc$V2 <- gsub("^ | $", "", df.proc$V2)
    meta <- rbind(df.acqu, df.proc)
    if (as.numeric(meta$V2[grep("p_BYTORDP", meta$V1)]) !=
        0) {
      endianness <- "big"
    }
    else {
      endianness <- "little"
    }
    spec <- readBin(rf[i], what = "int", n = as.numeric(meta$V2[grep("p_FTSIZE",
                                                                     meta$V1)]), size = 4, signed = T, endian = endianness)
    spec <- ((2^as.numeric(meta$V2[grep("p_NC_proc", meta$V1)])) *
               spec)
    swp <- as.numeric(meta$V2[grep("p_SW_p", meta$V1)])/as.numeric(meta$V2[grep("^p_SF$",
                                                                                meta$V1)])
    dppm <- swp/(length(spec) - 1)
    offset <- as.numeric(meta$V2[grep("OFFSET", meta$V1)])
    ppm <- seq(from = offset, to = (offset - swp), by = -dppm)
    return(list(meta, spec, ppm))
  })

  ids = unique(as.vector(unlist(sapply(out, function(x) x[[1]]$V1))))
  meta <- data.frame(t(sapply(out, function(x) {
    x[[1]] = unique(x[[1]])
    x[[1]]$V2[match(ids, x[[1]]$V1)]
  })), stringsAsFactors = F)
  colnames(meta) = ids

  nums = sapply(meta, function(x) {
    suppressWarnings(as.numeric(x))
  })
  idx = complete.cases(t(nums))
  meta[, idx] = nums[, idx]

  meta$a_Date =
    strptime("1970-01-01 00:00:00", format = "%Y-%m-%d %H:%M:%M") +
    as.numeric(meta$a_DATE)
  meta$a_RunOrder = rank(meta$a_Date, na.last = "min")
  rownames(meta) = r.exp
  ppm <- seq(meta$p_OFFSET[1], meta$p_OFFSET[1]-meta$a_SW[1], -(meta$a_SW[1]/meta$a_TD[1]))

  if(length(unique(meta$a_TD))>1){warning('Some experiments in this folder differ in digital resolution, you might want to filter these out.')}
  if(length(unique(round(meta$a_SW),4))>1){warning('The sweep width is different for some experiments in this folder, please make sure you read-in the correct experiments.')}

  X <- t(sapply(out, function(x, ppmt = ppm) {
    sInter = approxfun(x[[3]], x[[2]])
    sInter(ppmt)
  }))
  colnames(X) <- ppm

  # rownames should include basenames if these differ
  fnames=strsplit(dirname(afile), '/')
  rnames=do.call(rbind, fnames)
  idx=apply(rnames, 2, function(x) length(unique(x))==1)
  idx=which(idx==F)
  if(length(idx)>1){
    rownames(X) <- sapply(fnames, function(x){
      paste(x[idx[1]:length(x)], collapse='/')
    })
  }else{
    rownames(X)<-rnames[,idx]
  }


  # The offset is slightly different in each exp, leading to NA values for few ppm variables at the ends of some spectra. Remove these if not too many!
  idx=unique(which(is.na(X), arr.ind = T)[,2])
  if(length(idx)<300){
    X=X[,-idx]
    ppm=ppm[-idx] }else{
      warning('The ppm range is way to high for some experiments, please double check if you read-in the correct spectra!')

  }


  assign("X", X, envir = .GlobalEnv)
  assign("ppm", ppm, envir = .GlobalEnv)
  assign("meta", meta, envir = .GlobalEnv)
}
