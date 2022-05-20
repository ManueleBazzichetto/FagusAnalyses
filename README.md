# FagusAnalyses
This repo includes data and the R script for replicating the analyses of the distribution of Fagus sylvatica across west Europe (Italy, France and Spain).  The R script of the analyses provides a step-by-step guide on how to use the uniform sampling of the environmental space to collect background points.

Data used: [WorldClim](https://www.worldclim.org/data/index.html) bioclimatic variables (through the function raster::getData), sPlotOpen ([Sabatini et al., 2021](https://doi.org/10.1111/geb.13346)), and EU-Forest ([Mauri et al., 2017](https://doi.org/10.1038/sdata.2016.123)). The sPlotOpen and EU-Forest data version used in our analyses can be found in the R folder (in this repo). sPlotOpen and EU-Forest are open datasets and can be downloaded at the data sources reported in the related publications.

The R script is roughly divided in three sections:
1) Data from the above-mentioned sources are downloaded (WorldClim) and cleaned for the analyses.
2) The uniform sampling of the environmental space is implemented (using the UEsampling package) to i) apply kind of a spatial thinning procedure for subsetting presence data of Fagus sylvatica, and ii) collecting background points within the environmental space.
3) The obtained training and testing datasets are used to calibrate and validate, respectively, species ditribution models for Fagus sylvatica fitted using binary generalised linear models and random forest.
