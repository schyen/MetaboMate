
# [MetaboMate](https://kimsche.github.io/MetaboMate/)
MetaboMate is a R toolbox for performing data analysis in the field of NMR-based metabolic profiling.

## Getting started

### Instsallation
The MetaboMate package is currently hosted on GitHub. Installing it requires the developer tools R package:
```r
install.packages('devtools')
library(devtools)
```

The following command installs the MetaboMate package:

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



## Tutorials
Tutorials (also called vignettes in R) are distributed with the package (run `vignette(package='MetaboMate')`) and can also be accessed via the [MetaboMate](https://kimsche.github.io/MetaboMate/) website. Currently there are two tutorials availabe:

* Data Import and Pre-processing
* Mutlivariate Statistics




