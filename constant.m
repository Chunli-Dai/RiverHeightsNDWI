
%directory
macdir=[]; % Leave it blank for linux;Absolute main directory part.
% macdir='/Users/chunlidai/surge/';
currentdir=pwd;
%addpath(genpath(['/Users/chunlidai/Google Drive/OpticalImages/ArcticDEM/paper/ReCoordinates/']));


codedir='../rivergithub2/';
multidir='/fs/byo/howat-data2/river/2018dec03/ortho/';
multidir='/fs/byo/howat-data5/pgc_deliv/chunli/sag/deliv3/ortho_imagery_non-max_ona/';
% stripdir='/data3/ArcticDEM/region_*/strips/2m/';
stripdir='/Users/chunlidai/Google Drive/NASAHydro/manuscript2/data/';
stripdir='/home/dai.56/data2/river/2018dec03/setsm_results/strips/2m_filt001/';
stripdir='/fs/byo/howat-data2/ArcticDEM/region_34_alaska_north/strips/2m/';
tiledir=[macdir,'/fs/byo/howat-data3/ArcticDEMmosaics/'];
tiledirnew=tiledir; %to store new downloaded dems

%control parameters
method=1; %1 direct method; 2 imagery-altimetry method
flagcoreg=1; % 1 do coregistration % 0 do not apply coregistration
flagplot=0; % 1 plot; 0 do not plot
flagmono=1; %1 use mono images only; 2 use stereor
widthstat=0;  %buffer width of a priori river mask for calculating ndwi statistics, and the buffer of tile for calculation.
mons=5;mone=10; %mon>=5&mon<=10; %mons, start month of running/non-frozen rivers. %mone, end month of running water.
%mons=1;mone=12; %all season.
poolsize=1; %24;
flagsect=0;  %1 work on a section of river; 0 work on a gage
flagbuf=1; %apply the river buffer or not. 1 apply; 0 do not apply. check riverbuffStep3.png, then decide.
flaglowest=0; % use the designated dem as lowest stage DEM for retrieving height profiles. 0 don't use; 1 use; 2 use but compare stage.
idlowest=1; %give a number based the output ('Use these strip DEMs').

%control parameters for multispec.m
threshold=0.3; %Suggest value 0.5 or 0.3; % general threshold of NDWI for water classification.
probthre=50.;% threshold for water probability. %Best for migitating random coregistration offset: 50.
stdthres=0.5; % if NDWI STD > stdthres, discard the image. Suggest value 0.5
dmthres=0.4; % if mean_ocean - mean_land > dmthres, discard the image. Suggest value 0.6 for ocean; 0.4 for river;
%For tanana river, it seems lots of good images has dmthres < 0.6, try 0.5. Need to get a statistical value based on clouds.

comerrthres=0.25;%2;%0.2;% commission error threshold; if  comerr >= comerrthres, image is bad.
		%Based on 55 images test for sag river, suggest value 0.25.
resr=2;
%lakearea=1000*500*4/4; %remove clusters smaller than this.
lakearea=100*500*4/4; %remove clusters smaller than this.
cloudarea=1000*5*4/4; %1000*5 pixels (2m); refill holes smaller than this.
cntmin=25*25; %unit:pixels. The size of a priori land/ocean area should be big enough to ensure reliable statistical analysis of the histogram of the region. (Liu and Jezek, 2004)
novlmt=3; %if number of repeats <= novlmt, set the area as edges/void.

%tmp directory
picdir='./pic/';
odircoreg='./outcoreg/';

if ~exist(picdir,'dir')
  mkdir(picdir)
end
if ~exist(odircoreg,'dir')
  mkdir(odircoreg)
end


