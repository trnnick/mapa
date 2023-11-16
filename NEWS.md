# CHANGES IN MAPA VERSION 2.0.6 (12 NOV 2023)
- Fixed smooth related dependencies.

# CHANGES IN MAPA VERSION 2.0.5 (10 JUL 2021)
- Updated class ckecks.
- Updated links in help files.

# CHANGES IN MAPA VERSION 2.0.4 (05 JAN 2018)
- Updated statetranslate() to handle changes in smooth package (state name "seasonality" changed to "seasonal")
- By default all functions does not produce plots anymore.
- Default combination operator changed to w.mean.

# CHANGES IN MAPA VERSION 2.0.3 (04 Jul 2017)
- Minor improvement in plotting mapa.fit class objects. Now levels for which no seasonality is estimated are in grey.
- Accomondate changes in smooth::es output.

# CHANGES IN MAPA VERSION 2.0.2 (06 Jun 2017)
- Fixed bug in mapacalc when mapaest uses pr.comp=-1.
- Fixed bug in plotting th output of mapaest.
- Minor changes in input argument checking. 
- Added CITATION file.

# CHANGES IN MAPA VERSION 2.0.1 (01 Nov 2016)

## NEW FEATURES

- Introduced es from smooth package as core for MAPA. This can handle any seasonal length and allows some additional flexibility in model search.
- MAPAx implemented. Now you introduce exogenous variables using xreg input argument, which the user can pre-process using PCA.
- Introduced new combinations: "w.mean" and "w.median". These use mean or median for the combination of level and trend, but weight seasonality and xreg inversly proportional to the aggregation level. It is recommended to use these with high frequency time seres. 

## CHANGES

- mapaest, mapa & mapasimple: the user can select to use either ets (forecast) or es (smooth).
- Users can pass additional arguments to ets or es through ellipsis.
- plotmapa: became obsolete, simply use plot().
- mapacalc, mapasimple: the "output" input argument is removed.
- mapasimple: simpified function wrapper.
- Introduced S3methods.
- Improved visualisation. 
- Colour scheme from RColorBrewer package.
- smooth and RColorBrewer packages added to dependencies. 
- Switched ChangeLog to NEWS.md

# CHANGES IN MAPA VERSION 1.9.1 (24 May 2016)
- Fixed bug introduced by changes in the initialisation of ETS when n < 3*s and seasonality (s) = 2. The fix does not permit seasonal models in that case. If n >= 3*s seasonality is permitted. 

# CHANGES IN MAPA VERSION 1.9 (21 November 2014)
- Corrected bug in mapafor when ifh == 0 and prediction intervals were requested.
- Added more detailed package description.

# CHANGES IN MAPA VERSION 1.8 (03 October 2014)
- Corrected bug for very short time series.

# CHANGES IN MAPA VERSION 1.7 (31 July 2014)
- Added plotmapa function to visualise fitted MAPA without need to re-estimate
- Removed dependence on 'miscTools' package for calculation of medians by columns
- Corrections in .Rd files

# CHANGES IN MAPA VERSION 1.6 (03 July 2014)
- Added non-overlapping temporal aggregation function (tsaggr).
- Added URL info in description.
- Corrected bug in mapaest when actual was constant (ets returned different structure).

# CHANGES IN MAPA VERSION 1.5 (04 June 2014)
- Added empirical prediction intervals for mapa and mapafor.
- Added return to some functions that were missing it.
- Plots in mapaest and mapafor re-initialise the plot figure.
- Variable name insample replaced with y.
- Added more help details for mapasimple and mapacalc outplot input.
- Corrected bug in mapa that would not transfer model to mapaest.
- Corrected bug in mapafor when ppy was equal to 1.
- Corrected bug in mapasimple when fh was equal to 1. The bug was caused by the introduction of the model option in v1.4. 

# CHANGES IN MAPA VERSION 1.4 (13 May 2014)
- Added model option for mapa, mapasimple and mapaest.

# CHANGES IN MAPA VERSION 1.3 (15 April 2014)
- Corrected bug with minimumAL. Now minimumAL can be > 1.
- Corrected ylim when outplot=2 for the component plots.
- Added links to the help files. 

# CHANGES IN MAPA VERSION 1.2 (24 March 2014)	
- Stopped using forecast:::pegelsresid.C. New forecast package (v5.3) allows keeping initial values constant.

# CHANGES IN MAPA VERSION 1.1 (Not in CRAN)
- Switched from doSNOW to parallel.
- Introduced NAMESPACE.

# CHANGES IN MAPA VERSION 1.0 (Not in CRAN)
- Initial version.
