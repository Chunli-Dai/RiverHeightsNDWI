function [option]=SettingGeotiffOption(varargin)

Tinfo = varargin{1};
boundary = varargin{2};
gridspace = varargin{3};
%Radiometric resolution setting
BitDepth = varargin{4};
%Samples_per_pixel setting (band count)
SamplesPerPixel = varargin{5};

%Pixel size setting
ModelPixelScaleTag(1) = gridspace;
ModelPixelScaleTag(2) = gridspace;
ModelPixelScaleTag(3) = 0;

%Image top-left geo-coordinate
ModelTiepointTag(1) = boundary(1);
ModelTiepointTag(2) = boundary(4);
ModelPixelScaleTag(3) = 0;

GeoKeyDirectoryTag = Tinfo.GeoKeyDirectoryTag;

%Image origin setting
GeoDoubleParamsTag = zeros(1,4);
%Latitude setting
GeoDoubleParamsTag(1) = Tinfo.GeoDoubleParamsTag(1);
%Longitude setting
GeoDoubleParamsTag(4) = Tinfo.GeoDoubleParamsTag(4);

option.ModelPixelScaleTag = ModelPixelScaleTag;
option.ModelTiepointTag   = [0 0 0 ModelTiepointTag(1) ModelTiepointTag(2) 0] ;
option.GeoKeyDirectoryTag=GeoKeyDirectoryTag;
option.GeoDoubleParamsTag=GeoDoubleParamsTag;
option.bit_depth = BitDepth/SamplesPerPixel;

option.GTModelTypeGeoKey  = 1;
option.GTRasterTypeGeoKey = 1;
option.GeographicTypeGeoKey = 4326;

option.ProjectedCSTypeGeoKey = 32767; % user define
option.ProjectionGeoKey      = 32767; % user define
option.ProjCoordTransGeoKey = 15;     % polar stereographic
option.ProjNatOriginLatGeoKey = 70.00;
option.ProjNatOriginLongGeoKey = -45.00;
option.ProjFalseEastingGeoKey = 0;
option.ProjFalseNorthingGeoKey = 0;
option.ProjScaleAtNatOriginGeoKey = 1;
option.ProjStraightVertPoleLongGeoKey = -45.00;
option.SampleFormat = BitDepth/8;

option.GeogLinearUnitsGeoKey=9001;
option.ProjLinearUnitsGeoKey=9001;
option.VerticalUnitsGeoKey = 1;
%option.GeogAngularUnitsGeoKey=9107;
%option.GeogAzimuthUnitsGeoKey=9107;