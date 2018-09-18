
# Installation 
Installing the MetaboMate package from GitHub requires the developer tools R package:
```r
install.packages('devtools')
library(devtools)
```

The following command installs the MetaboMate package, including two tutorials.

**`devtools::install_github('kimsche/MetaboMate', build_vignettes=TRUE)`**

The installation can take up to about 3 minutes as the installation builds the tutorials, which come with high-resolution NMR example data.

#### Installation Notes for R Beginners

There might come up **warning messages** (eg. *time zone unknown*). These are not a problem! :pray:
However, **error messages** are a problem as these stop the installation process. 

#### Solutions to Common Installation Errors:

>`Installation failed: 'exdir' does not exist`

Workspace directory not writable, change the location of the current workspace

>`Dependency package missing`

Install the required package(s) from Bioconductor.org

You can check if the installation was successful by loading the package into your R session: `library(MetaboMate)` - all is fine if no error message is returned.

