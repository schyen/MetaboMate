#' @name bariatric
#' @docType data
#' @title Bariatric surgery data
#' @author Torben Kimhofer \email{tkimhofer@@gmail.com}
#' @usage data(bariatric)
#' @description Metabolic profiling NMR dataset of murine urine samples before and after undergoing bariatric surgery (Roux-en-Y gastric bypass).
#' @details Proton NMR spectra of urine samples were obtained using a Bruker 600 MHz spectrometer (Bruker, Rheinstetten, Germany) operating a frequency of 600.13 MHz at a temperature of 300 K. A standard NMR pulse sequence (\emph{recycle delay-90°-t1-90°-tm-90°-acquisition}) was applied to acquire one-dimensional Proton NMR spectral data, where \emph{t1} was set to 3 μs and \emph{tm} (mixing time) was set to 100 ms. Water suppression was achieved using selective irradiation during a \emph{recycle delay} of 2 s and \emph{tm}. A 90° pulse was adjusted to 10 μs. A total of 128 scans were collected into 64 k data points with a spectral width of 20 ppm.
#'
#' For further information on sample collection and processing please see reference below. Pre-processing of raw NMR data is described in the \emph{MetaboMate} vignette \strong{Data Import and Proprocessing}.
#' @format A list with three elements:
#' \describe{
#'   \item{X.pqn}{Matrix of pre-processed Proton NMR spectra with rows representing spectra}
#'   \item{ppm}{Chemical shift vector in ppm with its length equals the number of columns in \code{X.pqn}}
#'   \item{meta}{Spectrometer metadata}
#' }
#' @source Dr. Jia V. Li, Imperial College London, Department of Surgery and Cancer, Centre for Computational and Systems Medicine. \url{https://www.imperial.ac.uk/people/jia.li}
#' @references Li, Jia V., \emph{et al.} (2011) Metabolic surgery profoundly influences gut microbial-host metabolic cross-talk. \emph{Gut}. 60.9, 1214-1223.
#' @references \url{https://www.ncbi.nlm.nih.gov/pubmed/21572120}
#' @keywords datasets
"bariatric"
