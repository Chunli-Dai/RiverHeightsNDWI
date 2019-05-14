Same as https://github.com/Chunli-Dai/RiverHeightsNDWI, except in parallel.

Things to do manually:
5\ run qsub jobp48.pbs
4\ copy the entire folder rivercore to your local source folder from https://github.com/Chunli-Dai/RiverHeightsNDWI.
3\ Change constant.m, poolsize=24; %24;
2\ Chang Tilemain.m, line 115 to the line number of target station : for i=32
1\ Change Tilemain.m, line 9 to the directory of code: addpath(genpath(['/home/dai.56/arcticdemapp/river/rivergithub2v2/']));


########Same as RiverHeightsNDWI
1\ In .bashrc, add line export PATH=$PATH:setsmdir
where setsmdir is your directory of setsm code, e.g. /home/dai.56/arcticdemapp/river/rivergithub2/SETSM_coreg/
2\ Change the code directory in constant.m, e.g. addpath(genpath([currentdir,'/../rivergithub2/']));
3\ Check the image directory in constant.m
4\ Download USGS gage time series from website (e.g. https://waterdata.usgs.gov/nwis/uv?site_no=15908000) save as file usgsgage.txt.
   get usgsgagewidth.dat
5\ Edit Tilemain, let i be selected station number in "for i=32"; Run matlab< Tilemain.m
6\ To do: automatically select the direction in rivercenterline.

Step 1: Shoreline detection using entropy and brightness: 
maskentropy.m: water classification using brightness and optimal thresholding method (Gonzalez and Woods, 1992); need a priori water mask.
multispecmono.m: water classification using NDWI.

Step 2: Elevation extraction.
riverprof.m: data processing for extracting river heights. 

Step 3: Filtering and Fitting.
ProcessTananaFairbanks.m: data processing and plotting of Fig.2.

Other codes and subroutines see https://github.com/mikedurand/SmoothRiverElevations

Step 4: Plot river height time series and discharge time series.
gageheights.m: plotting the time series of river heights at the USGS gage in Fairbanks, Alaska (Fig 3a).
Preparation: 
gagefo.txt, river height time series at the gage for all seasons.
gageft.txt, same as above but for winter season only.
usgsgage.txt, USGS gage height time series.

stagedischarge.m: Plot Fig.3b for the discharge time series.
Preparation:
uv10to16.txt: data files for stage-discharge rating curve.
legs.m, legsd.m: computation of legendre function for fitting rating curve.

The comparison of river height time series by two different methods.
comparemethods.m


%% Other versions

pxpymonobp1.m:reference images are the mono multispectral images in strip files
pxpymono.m: reference images are the panchormatic orthoimage strip files.For the mono image of strip, make sure to use its own strip orthoimage as reference.
pxpymonobp2.m:reference images not fixed to its own strip files. Not used.

multispecstripbp1.m: assume mono image aligned to strip DEM, which turns out can be off by 2 m.
multispecstrip.m: use the translational parameters to align mono image to the strip DEM.

multispecmonobp1.m: use the mean of land mean and ocean mean.


Note:
cp ~/chunliwork/river/riverpic/run2/gageflagyx/gageheights_test.m widthgageprof.m

Version:  rivergithub2v2
Main change: do not save datar, which use too large memories.
Difference with rivergithub2 on results: datarsv is water mask interpolated to 2m resolutio. whileas, in rivergithub2v2, water mask is the original resolution, which is finer than 2m.

