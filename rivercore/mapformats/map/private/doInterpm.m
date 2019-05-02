function [lat, lon] = doInterpm(lat, lon, minimumSeparation, method, angleUnit)
% Core computations needed by INTERPM, MAPPROFILE, and VEC2MTX.

% Copyright 2006-2015 The MathWorks, Inc.

% Ensure column vectors.
lat = lat(:);
lon = lon(:);

% Compute max angular separation between each pair of adjacent vertices.
separation = max([abs(diff(lat))'; abs(diff(lon))'])';

% Find separations that exceed the specified limit.
indx = find(separation > minimumSeparation);

if ~isempty(indx)
    % Interpolate additional vertices where needed.
    steps = ceil(separation(indx)/minimumSeparation);
    switch method
        case {'gc','rh'}
            % Interpolate along great circles or rhumb lines.
            [lat, lon] = toRadians(angleUnit, lat, lon);
            interpfun = @(lat, lon, loc, step) interpTrack(lat, lon, loc, step, method);
            [lat, lon] = interpAll(lat, lon, indx, steps, interpfun);
            lon = unwrapMultipart(lon);
            [lat, lon] = fromRadians(angleUnit, lat, lon);
            
        otherwise
            % Interpolate in a planar system.
            [lat, lon] = interpAll(lat, lon, indx, steps, @interpPlanar);
    end
end

%--------------------------------------------------------------------------

function [latinsert, loninsert] = interpTrack(lat, lon, e, step, method)

[latinsert, loninsert] = doTrack(method, ...
    lat(e + [0 1]), lon(e + [0 1]), [1 0], step + 1);

% Remove trailing NaNs inserted by doTrack.
latinsert(isnan(latinsert)) = [];
loninsert(isnan(loninsert)) = [];

% Remove starting and ending points to avoid duplication when inserted.
latinsert = latinsert(2:length(latinsert)-1);
loninsert = loninsert(2:length(loninsert)-1);

%--------------------------------------------------------------------------

function [latinsert, loninsert] = interpPlanar(lat, lon, e, step)

factors = (1:step-1)';
latinsert = ((lat(e+1)-lat(e))/step)*factors + lat(e);
loninsert = ((lon(e+1)-lon(e))/step)*factors + lon(e);

%--------------------------------------------------------------------------

function [ilat, ilon] = interpAll(lat, lon, indx, steps, interpfun)

% Length of output vectors
N = numel(lat) + sum(steps) - numel(steps);

% Pre-allocate outputs.
ilat = zeros(N,1);
ilon = zeros(N,1);

s = 1;
sNew = 1;
for k = 1:numel(indx)
    % Insert vertices between lat(indx(k)), lon(indx(k))
    % and lat(1 + indx(k)), lon(1 + indx(k)).
    
    % Copy the elements that precede the current insertion point.
    e = indx(k);
    eNew = sNew + e - s;
    ilat(sNew:eNew) = lat(s:e);
    ilon(sNew:eNew) = lon(s:e);
    
    % Advance s to be ready for the next iteration (or final copy),
    % and also for use in the insertion step below.
    s = 1 + e;
    
    % Compute the new vertices that are to be inserted
    [latinsert, loninsert] = interpfun(lat, lon, e, steps(k));
    
    % Copy the new vertices into the output arrays.
    sNew = eNew + 1;
    eNew = eNew + numel(latinsert);
    ilat(sNew:eNew) = latinsert;
    ilon(sNew:eNew) = loninsert;
    
    % Advance sNew to be ready for the next iteration (or final copy).
    sNew = eNew + 1;
end

% Append any vertices that come after the last insertion point. Note
% that if indx is empty, we'll have skipped the loop, s will equal 1, and
% we'll be making an exact copy of the inputs.
ilat(sNew:end) = lat(s:end);
ilon(sNew:end) = lon(s:end);
