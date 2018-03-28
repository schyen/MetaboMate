#' Reading 1D NMR spectra from Bruker
#' @export
#' @description  Reading 1D NMR spectra from Bruker. This function is recurisvely searching for matching processing, acquisition and spectrum data (r1) files in a given experiment folder.
#' @param path Path to experiment folder.
#' @param filter Logical, if set TRUE then experiments with matching \emph{acqus}, \emph{procs} and \emph{1r} files are read-in. If set FALSE, the algorihm will stop if file systems are corrupt or 2D experiments are presenr. The latter option is for more for running file system tests.
#' @return The following variables are automatically returned as globally defined variables:
#' \item{X}{Imported NMR matrix where NMR spectra in rows and ppm variables in columns.}
#' \item{ppm}{Ppm vector.}
#' \item{meta}{Spectrometer metadata extracted from acqus and procs files (eg time of spectral acquisition).}
#' @importFrom stats approxfun
#' @author Torben Kimhofer, Imperial College London (2017)
#' @details Paremeters with prefix \code{a} or \code{p} were extracted from the \emph{acqus} or \emph{procs} file, respectively. Paremeters without any prefix were calculated on the fly.


readBruker=function(path, filter=F)
{
  warnDef <- options("warn")$warn
  warnRead <- options(warn = -1)

  if(grepl('/$', path)){path <- gsub('/$', '', path)}

  # searches for acquisition files
  afile <- list.files(path = path, pattern = "^acqus$",
                        all.files = FALSE, full.names = TRUE, recursive = TRUE,
                        ignore.case = TRUE)

  # searches for parameter files
  pfile <- list.files(path = path, pattern = "^procs$",
                      all.files = FALSE, full.names = TRUE, recursive = TRUE,
                      ignore.case = TRUE)
  pfile <- grep('pdata/1/', pfile, value=T)

  # searches for 1r files (real part of spectrum)
  rfile <- list.files(path = path, pattern = "^1r$", all.files = FALSE,
                      full.names = TRUE, recursive = TRUE, ignore.case = TRUE)
  rfile <- grep('pdata/1/', rfile, value=T)


  Lp<-length(pfile)
  if (Lp == 0 || length(rfile) == 0 || length(afile)==0) {
    stop("No existing Bruker files in specified path.")
  }

  a.exp<-gsub(paste('^', path, '/|/acqus$', sep=''),'', afile)
  p.exp<-gsub(paste('^', path, '/|/pdata/1/procs$', sep=''),'', pfile)
  r.exp<-gsub(paste('^', path, '/|/pdata/1/1r$', sep=''),'', rfile)

  # check files match
  idx.a=a.exp %in% r.exp & a.exp %in% p.exp
  idx.p=p.exp %in% r.exp & p.exp %in% a.exp
  idx.r=r.exp %in% a.exp & r.exp %in% p.exp

  # filter acqus and procs for 1r files (exclude 2D experiments and corrupt file systems)
  if(filter==T){
    afile=afile[idx.a]
    pfile=pfile[idx.p]
    rfile=rfile[idx.r]
  }else{
    not.present<-list()
    not.present[[1]]=afile[!idx.a]
    not.present[[2]]<-pfile[!idx.p]
    not.present[[3]]<-rfile[!idx.r]

    if(length(not.present[[1]])>0){
      cat('File structure depreciated: missing Bruker acquisition file(s). See output for file names.\n')
    }
    if(length(not.present[[2]])>0){
      cat('File structure depreciated: missing Bruker processing file(s). See output for experiment folder names.\n')
    }
    if(length(not.present[[3]])>0){
      cat('File structure depreciated: missing aquisition file(s). See output for experiment folder names.\n')
    }

    if(sum(sapply(not.present, length))>0) {
      return(sort(unlist(not.present)))}

  }

  Lp<-length(pfile)
  if (Lp == 0) {
    stop("No matching Bruker files in specified path.")
  }


  if(Lp==1){cat('Reading ', Lp, ' spectrum.\n', sep='')}else{cat('Reading ', Lp, ' spectra.\n', sep='')}

    out<-apply(rbind(1:Lp), 2, function(i, pf=pfile, rf=rfile, af=afile){

      con <- file(af[i], open = "r")
      aLine <- readLines(con, n = -1, warn = FALSE)
      close(con)
      aLine<-aLine[-c(1:9)]
      aLine<-aLine[!grepl('##END', aLine)]

      idx<-grep('#', aLine)
      aLineC<-array()
      for(j in 1:length(idx)){
        if(grepl('#', aLine[idx[j]+1])){
          aLineC[j]<-aLine[idx[j]]
          next}
        end_idx<-grep('#', aLine[(idx[j]+1):length(aLine)])[1]
        aLineC[j]<-paste(aLine[idx[j]:(idx[j]+end_idx-1)], collapse=' ')
      }

      aLineC<-gsub(' \\(.*\\)', '', aLineC)
      myV <- strsplit(aLineC, "=")

      df.acqu<-unique(as.data.frame(do.call(rbind, myV), stringsAsFactors = F))
      df.acqu$V1<-paste('a', gsub('##\\$', '', df.acqu$V1), sep='_')
      df.acqu$V2<-gsub('^ | $', '', df.acqu$V2)


      # procs files (processing status parameters)
      con <- file(pf[i], open = "r")
      aLine <- readLines(con, n = -1, warn = FALSE)
      myV <- strsplit(aLine, "=")
      close(con)
      idx<-grepl('##', aLine) & sapply(myV, length)==2
      myV<-myV[idx]

      df.proc<-unique(as.data.frame(do.call(rbind, myV[1:(length(myV)-1)]), stringsAsFactors = F))
      df.proc$V1<-paste('p', gsub('^##[\\$]?', '', df.proc$V1), sep='_')
      df.proc$V2<-gsub('^ | $', '', df.proc$V2)

      meta<-rbind(df.acqu, df.proc)

      # read-in spec
      if(as.numeric(meta$V2[grep('p_BYTORDP', meta$V1)])!=0){endianness<-'big'}else{endianness<-'little'}

      # read binary real part of spec
      spec <- readBin(rf[i],
                   what = "int",
                   n = as.numeric(meta$V2[grep('p_FTSIZE', meta$V1)]),
                   size = 4,
                   signed = T,
                   endian = endianness)

      spec <- ((2^as.numeric(meta$V2[grep('p_NC_proc', meta$V1)])) * spec)

      swp <- as.numeric(meta$V2[grep('p_SW_p', meta$V1)])/as.numeric(meta$V2[grep('^p_SF$', meta$V1)]) # Spectral acquisition width in ppm
      dppm <- swp/(length(spec) - 1) # how much difference is there between each ppm variable (swp/number of spectral points -1)

      offset <- as.numeric(meta$V2[grep('OFFSET', meta$V1)])
      ppm <- seq(from=offset, to=(offset - swp), by = -dppm)

      return(list(meta, spec, ppm))
    })

    ids=unique(unlist(sapply(out, function(x) x[[1]]$V1)))
    meta<<-data.frame(t(sapply(out, function(x) {
      x=unique(x)
      x[[1]]$V2[match(ids, x[[1]]$V1)]
  })), stringsAsFactors = F)
    colnames(meta)=ids

    meta$a_Date=strptime('1970-01-01 00:00:00', format='%Y-%m-%d %H:%M:%M')+as.numeric(meta$a_DATE)
    meta$a_RunOrder=rank(meta$a_Date, na.last = 'min')

    rownames(meta)=r.exp

  ppm<-seq(14.5, -5, -0.0001526)
  X<-t(sapply(out, function(x, ppmt=ppm){
    sInter=approxfun(x[[3]],x[[2]])
    sInter(ppmt)
  }))

  colnames(X)<-ppm
  rownames(X)<-r.exp

  # return data
  assign('X', X, envir=.GlobalEnv)
  assign('ppm', ppm, envir=.GlobalEnv)
  assign('meta', meta, envir=.GlobalEnv)

}
