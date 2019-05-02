function [option, bbox]=geotiff_option(filename,mapx,mapy)

Tinfo = imfinfo(filename);
[lon, lat]=ps2wgs(fbnd(:,1), fbnd(:,2));
minlon = min(lon);
maxlon = max(lon);
minlat = min(lat);
maxlat = max(lat);
bbox = [minlon minlat;maxlon maxlat];



    
