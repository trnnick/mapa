# CHANGES IN MAPA VERSION 2.1.0 (07 NOV 2026)
- new: New plot methods for mapa.calc and mapa.for. Plots accept user arguments, such as main, xlab, ylab, xlim, ylim, col and lwd, use the time index of ts inputs, and respect the user's graphical layout.
- new: plot.mapa.calc can plot all or selected components.
- new: New summary methods for mapa.fit (model and parameters at each aggregation level), mapa.calc (combined components for each horizon) and mapa.for (in-sample ME, MAE and RMSE for each horizon).
- new; Forecasts of ts inputs are returned as ts objects.
- new: mapaest, mapa, mapasimple, mapafor: new input cl, to provide an existing cluster for parallel processing.
- new: mapafor and mapa: parallel processing is also used for the in-sample forecasts. mapa initialises a single cluster for estimation and forecasting.
- new: tsaggr: additional arguments are passed to plot.
- new: New function mapafit.old, to transform the output of mapaest to the structure of earlier versions. This is provided for compatibility and will be removed in a later version.
- new: Added tests.
- mapaest: the fitted MAPA is now a list that keeps the complete fitted model at each aggregation level, instead of a list-matrix. This makes MAPA robust to changes in the forecast and smooth packages. Code that indexes the fitted MAPA directly needs to be updated, or can use mapafit.old in the meantime. Fitted MAPA objects from earlier versions are still accepted.
- The default forecast horizon is now fh=NULL in all functions, which gives a full season. In mapafor and mapa fh=0 still gives only in-sample forecasts. The earlier defaults (fh=-1 in mapafor and fh=0 in mapacalc) are still accepted.
- forecast, smooth and RColorBrewer are now imported instead of attached. Code that relied on library(MAPA) to attach the forecast package needs to load it separately. Requires R >= 4.1.0.
- Prediction intervals are printed in order, from the lowest to the highest bound.
- Invalid options give informative errors.
- Internal code simplified, with repeated calculations moved to helper functions.
- Documentation and NAMESPACE are now generated with roxygen2 from the R files. Added package level help page (?MAPA).
- Expanded documentation and examples, covering the use of xreg and principal components, the use of ets (forecast) and es (smooth), and the other arguments.
- tsaggr now uses colours from Rcolorbrewer as the rest of the package does. 
- The default prediction interval shading in plot.mapa.for now follows the forecast colour.
- buggfix: mapaest: MAPAx (type="es" with xreg) failed whenever an aggregation level could not be estimated due to sample size, which is common when xreg is used.
- buggfix: mapacalc, mapafor: failed when the lowest requested aggregation level could not be estimated. Only the estimated aggregation levels are now used, with the correct weights.
- buggfix: mapacalc: with pr.comp != 0 the principal components were calculated again at forecast time, giving components that did not match the ones used in estimation. The rotation from estimation is now stored and reused. This changes MAPAx forecasts that use pr.comp.
- buggfix: mapacalc: the aggregation of future xreg values did not include the start offset, so the last aggregated future period used fewer observations than needed. This changes MAPAx forecasts at longer horizons for some aggregation levels.
- buggfix: print.mapa.fit: printed the position of the aggregation level instead of its value when minimumAL > 1, and failed when some levels were not estimated.
- buggfix: plot.mapa.fit: the grey shading for levels without seasonality was misplaced when minimumAL > 1.
- buggfix: Component plot: the combined level and trend shown for comb="w.mean" and "w.median" was the weighted one, instead of the mean or median used for the forecasts.
- buggfix: print.mapa.for: garbled output when ifh > 1.
- buggfix: mapafor: failed when fh=0 and conf.lvl was given.
- buggfix: Plot limits were not symmetric in mapafor and tsaggr.

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
