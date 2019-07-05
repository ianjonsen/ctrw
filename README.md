# ctrw

[![Project Status: Abandoned – Initial development has started, but there has not yet been a stable, usable release; the project has been abandoned and the author(s) do not intend on continuing development.](https://www.repostatus.org/badges/latest/abandoned.svg)](https://www.repostatus.org/#abandoned)

# This package has been merged into [`foieGras`](https://github.com/ianjonsen/foieGras) all future development will happen there

**ctrw** - a Continuous-Time Random Walk state-space model to filter Argos Least Squares or Kalman Filter location data

`ctrw` is an R package that fits a continuous-time random walk model in state-space form to filter Argos satellite location data. Template Model Builder (`TMB`) is used for fast estimation. The Argos data can be either (older) Least Squares-based locations or (newer) Kalman Filter-based locations with error ellipse information. Separate measurement models are used for these two data types. The model estimates two sets of location states: 1) corresponding to each observation, which are irregularly timed; and 2) corresponding to regular time intervals specified by the user. Locations are estimated and output returned on the Mercator projection. 
