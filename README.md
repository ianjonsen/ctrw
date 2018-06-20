# ctrw

[![Project Status: WIP – Initial development is in progress, but there has not yet been a stable, usable release suitable for the public.](http://www.repostatus.org/badges/latest/wip.svg)](http://www.repostatus.org/#wip)
[![Travis-CI Build Status](https://travis-ci.org/ianjonsen/ctrw.svg?branch=master)](https://travis-ci.org/ianjonsen/ctrw)

**ctrw** - a Continuous-Time Random Walk state-space model to filter Argos Least Squares or Kalman Filter location data

`ctrw` is an R package that fits a continuous-time random walk model in state-space form to filter Argos satellite location data. Template Model Builder (`TMB`) is used for fast estimation. The Argos data can be either (older) Least Squares-based locations or (newer) Kalman Filter-based locations with error ellipse information. Separate measurement models are used for these two data types. The model estimates two sets of location states: 1) corresponding to each observation, which are irregularly timed; and 2) corresponding to regular time intervals specified by the user. Locations are estimated and output returned on the Mercator projection. 

Read `?fit_ssm` for details and an example of how to use the package 

## Installation
First, ensure you have R version >= 3.3.0 installed:

```R
R.Version()
```

On PC's running Windows, ensure you have installed [Rtools](https://cran.r-project.org/bin/windows/Rtools/) 

On Mac's, ensure you have installed [Xcode](https://developer.apple.com/xcode/) and Xcode developer tools. If installation is needed, make sure you start Xcode after install to ensure final setup of developer tools is completed. Both Xcode and Xcode developer tools can be installed from the [Mac App Store](https://itunes.apple.com/au/app/xcode/id497799835?mt=12)

Next, you will need to install `TMB` and it's dependencies from within R:
```R
install.packages("TMB")
```

Then install `devtools` and it's dependencies and finally install `ctrw` from GitHub:

```R
install.packages("devtools")  
devtools::install_github("ianjonsen/ctrw")
```
