Mutliple Aggregation Prediction Algorithm (MAPA)
=======
[![GitHub version](https://badge.fury.io/gh/trnnick%2Fmapa.svg)](https://badge.fury.io/gh/trnnick%2Fmapa)
[![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/MAPA?color=blue)](https://CRAN.R-project.org/package=MAPA)
[![Downloads](http://cranlogs.r-pkg.org/badges/MAPA?color=blue)](https://CRAN.R-project.org/package=MAPA)


Development repository for the MAPA package for R.
Stable version can be found at: https://cran.r-project.org/web/packages/MAPA/index.html

<img src="https://github.com/trnnick/mapa/blob/4f3ecc96c06064a91e6f26651493c72ac9120abd/mapa-hex.PNG" height="150"/>

Functions and wrappers for using the Multiple Aggregation Prediction Algorithm (MAPA) for time series forecasting. MAPA models and forecasts time series at multiple temporal aggregation levels, thus strengthening and attenuating the various time series components for better holistic estimation of its structure. 

More details can be found at: http://tinyurl.com/mapa-paper


To install use:

> if (!require("devtools")){install.packages("devtools")}

> devtools::install_github("trnnick/mapa")

_Happy forecasting!_
