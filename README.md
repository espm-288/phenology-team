# Forecasting greenness with NEON data

Materials for [EFI-RCN NEON phenology forecasting challenge](https://ecoforecast.org/efi-rcn-forecast-challenges/) for team "greenbears." The goal of this challenge is to predict gcc_90 values (the 90th quantile of the common greenness index) 35 days out. NEON measures gcc_90 values every day at 8 target sites.

This repository contains code for retrieving NEON phenology data, generating predictions using generalized additive models (GAMs) or seasonal trend decomposition using LOESS (STL), and submitting data to the NEON contest.

The main script, which is run daily using GitHub Actions, is `automated_par.R`. This script uses GAMs to estimate an annual cycle in the gcc greeness index and supplements this prediction with an effect of photosynthetically active radiation (PAR) which itself is also forecasted using a separate GAM. It also contains code to retrieve up-to-date data inputs and to publish predictions to NEON.

The script `automated_workflow.R` contains an equivalent method but without looking at the effect of the PAR covaraite.




