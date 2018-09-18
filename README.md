# Installation 

The MetaboMate package is currently hosted on GitHub. Installing it requires the developer tools R package:
```r
install.packages('devtools')
library(devtools)
```

Now, the following command installs the MetaboMate package:

**`devtools::install_github('kimsche/MetaboMate', build_vignettes=TRUE)`**

The installation process can take up to about 3 minutes as the tutorials are build locally with high-resolution NMR example data.


&nbsp;

#### Installation Notes for R Starters

There might come up **warning messages** (eg. *unknown time zone*). These are not a problem! :pray:
However, **error messages** are a problem as these stop the installation process. 

&nbsp;

#### Solutions to Common Installation Errors:

>`Installation failed: 'exdir' does not exist`

Workspace directory not writable, change the location of the current workspace

>`Dependency package missing`

Install the required package(s) from Bioconductor.org

&nbsp;

You can check if the installation was successful by loading the package into your R session: `library(MetaboMate)` - all is fine if no error message is returned.

