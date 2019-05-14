This is the update for the work in the paper by Dai et al. (Estimating river surface elevation from ArcticDEM, GRL, 2018). Here the main change is that we use mainly multispectral images rather than panchromatic images. The NDWI thresholding method is used for water classification.

Things to do manually:
1\ In .bashrc, add line export PATH=$PATH:setsmdir
where setsmdir is your directory of setsm code, e.g. /home/dai.56/arcticdemapp/river/rivergithub2/SETSM_coreg/
2\ Change the code directory in Tilemain.m, e.g. addpath(genpath([currentdir,'/../rivergithub2/']));
3\ Check the image directory in constant.m
4\ Download USGS gage time series from website (e.g. https://waterdata.usgs.gov/nwis/uv?site_no=15908000) save as file usgsgage.txt (Format:yyyy-mm-dd HH:MM height(feet)).
   get usgsgagewidth.dat (Format:yyyy-mm-dd HH:MM discharge(ft^3/s) channel_width(ft) gage_height_va(ft))
  See below for details;
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

####################################################

## Steps for getting usgs height and width time series:
a\ Getting usgs height
https://waterdata.usgs.gov/nwis/uv?site_no=15908000
https://nwis.waterdata.usgs.gov/nwis/inventory ; search lat lon
click " Current / Historical Observations "
Choose Discharge, Gage Height, Tab-separated, Begin date and End date; then click go
awk -F'\t' '{if ($7!="") print $3, $7}' uv10to16.txt > usgsgage.txt  #tab separated; ^_^; manually checked

b\ USGS gagewidth steps:
1\ download data from https://waterdata.usgs.gov/nwis/measurements?site_no=15908000&agency_cd=USGS&format=html_table_expanded (search field measurements from https://waterdata.usgs.gov/nwis/inventory/?site_no=15908000&agency_cd=USGS )
2\ download data with tab-separated data (measurementstab.txt)
   click "Tab-separated data with channel data"
3\ import to excel, and then edit, get usgsgagewidth.xlsx
   when import, choose Data-> from text->Delimited; started from row 15;
   delete columns except (measurement_dt  discharge_va    chan_width      gage_height_va)
   change date format to be (yyyy-mm-dd hh:mm)
4\ get txt format, usgsgagewidth.dat 
% measurement_dt  discharge_va    chan_width      gage_height_va




