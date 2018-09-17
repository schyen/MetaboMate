# MetaboMate
The R toolbox for NMR-based metabolic profiling.

## MetaboMate
This package is all you need to perform data analysis in the filed of NMR-based metabolic profiling. It covers all basic processing steps, which includes importing 1D spectra from Bruker file format, data pre-processing (referencing, baseline correction, excision of spectral areas, etc.), spectral normalisation methods (eg., probablistic quotient normalisation), unsupervised and supervised multivariate statistis combined with aestethically pleasing visualisations using ggplot2. Many of the provided functions can also be used for mass spectrometry data analysis. All functions are well documented.

## Installation 
Installing the MetaboMate package from GitHub requires the developer tools R package:
```r
install.packages('devtools')
library(devtools)
```

The following command installs the MetaboMate package, including tutorials. The installation can take up to about 3 minutes as the installation builds the tutorials, which come with high-resolution NMR example data.

`devtools::install_github('kimsche/MetaboMate', build_vignettes=TRUE)`


### Installation Notes for R beginners

There might come up **warning messages** (eg. 'unknown time zone'). These are not a problem! :pray:
However, **error messages** are a problem as these stop the installation process. 

#### Solutions to common error messages: 

>`Installation failed: 'exdir' does not exist`

Workspace directory not writable, change location of the workspace

>`Dependency package missing`

Install the required package(s) from Bioconductor.org

You can check if the installation was successful by loading the package into your R session: `library(MetaboMate)` - all is fine if no error message is returned.


## Tutorials / Vignettes
Tutorials are also called vignettes in R. MetaboMate comes currently with two tutorials describing

1. Data Import and Preprocessing,
2. Multivariate Statistical Analysis.

To open these on your computer, type the following commands into your R console:
```r
vignette('Data_Import_and_Preprocessing')
vignette('Multivariate_Statistical_Analysis')
```

## Feedback
Please let me know what you have in mind: tkimhofer@gmail.com
