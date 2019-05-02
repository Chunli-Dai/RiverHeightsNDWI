function test

% target DEM filename
target = '/data/aster/unzip_folder/AST14DMO_00307102005142948_20081010195344_16272/AST14DMO_00307102005142948_20081010195344_16272_DEM_ps.tif';
% reference DEM filename for loading Tinfo
reference_dem = '/data/GimpDEM/gimpdem1_1.tif';

TarImg              = readGeotiff(target);
Tinfo               = imfinfo(reference_dem);

% finding reference DEM overlapped with target DEM
tempboundary(1)     = min(TarImg.x);
tempboundary(2)     = min(TarImg.y);
tempboundary(3)     = max(TarImg.x);
tempboundary(4)     = max(TarImg.y);
tempRefImg          = subsetGimpDEM(tempboundary(1),tempboundary(3),tempboundary(2),tempboundary(4));
RefImg              = tempRefImg;

% target DEM gridspace
tG                  = TarImg.info.map_info.dx;  
% reference DEM gridspace
rG                  = 30;

%% coregistration function - case 1 : using default values (no input params)
% out                 = coregisterdems(TarImg, RefImg, [tG rG]);

%% coregistration function - case 2 : input params
%  params attributes
%  1. I : Initial 7 parameters for adjustment
%   params.I[s w p k tx ty tz]
%  If these values don't be known, use as following values.
%   params.I[1 0 0 0 0 0 0] (default)
%  2. C : Stop criteria (criteria for stopping the iterative least squares adjustment)
%   params.C[MI Mds Mdr Mds MdS]
%   where, (default values)
%   MI = 100 (default): the number of maximum iteration
%   Mds = 0.001 (default): the maximum correction value of scale
%   Mdr = 0.001 (default): the maximum correction value of rotations (unit = deg)
%   Mds = 0.05  (default): the maximum correction value of shifts (unit = m)
%   MdS = 0.0001 (default): the maximum changing-rate of sigma0
%  3. G : General information
%   params.G[MH AH Method]
%   where,
%   MH = 3000 (default) : Maximum Height value for removing null heights such as clouds or false matching value (unit = m)
%   AH = 20   (default) : DEM height accuracy (unit = m)  
%  4. M : Method of surface representation. Select between '2.5D' and '3D'
%   params.M = 3D (default)
%  5. V : minimum similarity value for selecting control surfaces
%   params.V[weight minimum_angle minimum_undulation]
%   where,
%   weight = 10 (default): maximum percentage of weights accumulated from right end on weight histogram
%   minumum_angle = 20 (default) (unit = deg) : minimum angle of surface slope for selecting control surface
%   minumum_undulation = 10 (default) (unit = m): minimum undulation of height differences within windows mask (3X3) for selecting control surface

params.I            = [1 0 0 0 0 0 0];
params.C            = [100 0.001 0.001 0.05 0.0001];
params.G            = [3000 20];
params.M            = '3D';
params.V            = [10 20 10];

out                 = coregisterdems(TarImg, RefImg, [tG rG], params);
%%

% results
Zc                  = out{1};
Zcs                 = out{2};
Zref                = out{3};
stats               = out{4};

% save tiff files
Boundary(1)     = min(Zc.x);
Boundary(2)     = min(Zc.y);
Boundary(3)     = max(Zc.x);
Boundary(4)     = max(Zc.y);

[minlon, minlat]                    = ps2wgs(Boundary(1), Boundary(2));   
[maxlon, maxlat]                    = ps2wgs(Boundary(3), Boundary(4));
bbox                                = [minlon minlat;maxlon maxlat];
%Radiometric resolution setting
BitDepth            = 32;
%Samples_per_pixel setting (band count)
SamplesPerPixel     = 1;
%Geotiff option setting
option              = SettingGeotiffOption(Tinfo, Boundary, rG, BitDepth, SamplesPerPixel);
geotiffwrite('/scratch/ngnmj_test/coregistered_DEM.tif', bbox, Zc.z,option.bit_depth,option);
geotiffwrite('/scratch/ngnmj_test/control_surfaces.tif', bbox, Zcs.z,option.bit_depth,option);
geotiffwrite('/scratch/ngnmj_test/reference_DEM.tif', bbox, Zref.z,option.bit_depth,option);
