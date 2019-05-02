function [I]= readGeotiff(varargin)
% readGeotiff: read geotiff using imread and assign map info from infinfo.
%
% I = GEOTIFF_READ('filename');
% Reads whole images
% I = GEOTIFF_READ('filename','pixel_subset', [minrow maxrow mincol maxcol]);
% I = GEOTIFF_READ('filename','map_subset'  , [minx maxx miny maxy]);
% extract subset of the specified.
% output:
% I.z, image data
% I.x, x coordinate in map
% I.y, y coordinate in map
% I.info, misc. info
% imshow(I.z, 'xdata', I.x, 'ydata', I.y);
% shows image with map coordinate
% Version by Yushin Ahn, ahn.74@osu.edu
% Glacier Dynamics Laboratory, 
% Byrd Polar Resear Center, Ohio State University 
% Referenced enviread.m (Ian Howat)

name = varargin{1};

Tinfo       = imfinfo(name);
info.cols   = Tinfo.Width;
info.rows   = Tinfo.Height;
info.imsize = Tinfo.Offset;
info.bands  = Tinfo.SamplesPerPixel;

info.map_info.dx    = Tinfo.ModelPixelScaleTag(1);
info.map_info.dy    = Tinfo.ModelPixelScaleTag(2);
info.map_info.mapx  = Tinfo.ModelTiepointTag(4);
info.map_info.mapy  = Tinfo.ModelTiepointTag(5);

subrows = [1 info.rows];
subcols = [1 info.cols];

minx = info.map_info.mapx;
maxy = info.map_info.mapy;

x = minx + ((0:info.cols-1).*info.map_info.dx);
y = maxy - ((0:info.rows  -1).*info.map_info.dy);
 

if nargin == 3;
    
    if strcmp(varargin{2},'pixel_subset');
        subrows = varargin{3}(1:2);
        subcols = varargin{3}(3:4);
        
    elseif strcmp(varargin{2},'map_subset');
        
        map_subset  = varargin{3};
        subcols = (map_subset(1:2)-info.map_info.mapx)./info.map_info.dx+1;
        subrows = (info.map_info.mapy - map_subset([4,3]))./info.map_info.dy+1;
        subcols= round(subcols);
        subrows = round(subrows);
        
        subcols(subcols < 1) = 1;
        subrows(subrows < 1) = 1;
        subcols(subcols > info.cols) = info.cols;
        subrows(subrows > info.rows) = info.rows;
    end

end
       
I.x = x(subcols(1):subcols(2));
I.y = y(subrows(1):subrows(2));

I.z = imread(name,'PixelRegion',{[subrows],[subcols]});

I.info = info;





