# MetaboMate
The R toolbox for NMR-based metabolic profiling

## R Installation 
Installing MetaboMate package from GitHub requires the developer tools R package:
```r
install.packages('devtools')
library(devtools)
```

The following command installs the MetaboMate package, including tutorials. The installation can take up to about 3 minutes as the installation builds the tutorials, which come with high-resolution NMR example data.

`devtools::install_github('kimsche/MetaboMate', build_vignettes=TRUE)`

> For R beginners: There might come up **warning messages** (eg. 'unknown time zone'). These are not a problem. However, **error messages** are a problem as these stop the installation process. 
>
>Solutions to common error messages: 
>
> **Installation failed: 'exdir' does not exist** - Workspace directory not writable, change location of the workspace
>
> ***Dependency package missing*** - Install the required package from CRAN or Bioconductor.org


## Tutorials / Vignettes
Tutorials are also called vignettes in R. MetaboMate comes currently with two tutorials describing

1. Data Import and Preprocessing
2. Multivariate Statistical Analysis

To open these on your computer, type the following commands:

`vignette('Data_Import_and_Preprocessing')`
`vignette('Multivariate_Statistical_Analysis')`


`vignette(package='MetaboMate')`

You can open these by simply clicking on them. For 

