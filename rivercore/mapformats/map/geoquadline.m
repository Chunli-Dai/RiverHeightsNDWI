function [latlim, lonlim] = geoquadline(lat, lon)
%GEOQUADLINE Geographic quadrangle bounding multi-part line
%
%   [LATLIM, LONLIM] = GEOQUADLINE(LAT, LON) returns the limits of the
%   tightest possible quadrangle that bounds a line connecting vertices
%   with geographic coordinates given by vectors LAT and LON.  The line may
%   be broken into multiple parts, delimited by values of NaN.  LATLIM is a
%   latitude limit vector of the form [southern_limit northern_limit].
%   LONLIM is a longitude limit vector of the form [western_limit
%   eastern_limit].  The elements of LONLIM are wrapped to the interval
%   [-180 180] and are not necessarily in ascending order.  All angles are
%   in degrees.
%
%   See also BUFGEOQUAD, GEOQUADPT, OUTLINEGEOQUAD

% Copyright 2012 The MathWorks, Inc.

% Work with row vectors.
lon = lon(:)';
lat = lat(:)';

% Keep only points for which both coordinates have finite values.
discard = isnan(lon) | isinf(lon) | isnan(lat) | isinf(lat);
lat(discard) = NaN;
lon(discard) = NaN;
[lat, lon] = removeExtraNanSeparators(lat, lon);

if all(isnan(lat))
    % No vertices were provided; meaningful limits cannot be computed.
    % (The empty input case ends up here, as well.)
    latlim = [];
    lonlim = [];
else
    % We expect -90 <= lat <= 90, but clamp latlim just in case.
    latlim = [max(min(lat),-90) min(max(lat),90)];
    
    % Unwrap each part and determine its western and eastern limits.
    lon = unwrapMultipart(lon,'degrees');
    [first, last] = internal.map.findFirstLastNonNan(lon);
    wlim = zeros(size(first));
    elim = wlim;
    for k = 1:numel(first)
        segment = lon(first(k):last(k));
        wlim(k) = min(segment);
        elim(k) = max(segment);
    end
    
    % There are 3 cases to consider ...
    if any(elim - wlim >= 360)
        % At least one part wraps all the way around the sphere.
        lonlim = [-180 180];
    elseif isscalar(wlim)
        % There's only one part, and it doesn't wrap the sphere.
        lonlim = wrapTo180([wlim elim]);
    else
        % There are at least two parts, and none wrap the sphere.
        lonlim = lonlimMultipart(wlim, elim);
    end
    
    if lonlim(1) ~= lonlim(2)
        if lonlim(1) == 180
            lonlim(1) = -180;
        end
        if lonlim(2) == -180
            lonlim(2) = 180;
        end
    end
    
end

%--------------------------------------------------------------------------

function lonlim = lonlimMultipart(wlim, elim)
% Longitude limits from western and eastern limit vectors, wlim and elim.
% The number of parts, numel(wlim) == numel(elim), is 2 or greater.

% Wrap part limits to 360, then sort in order of increasing wlim,
% before starting to identify gaps (intervals in longitude in which no
% lines are present).
wlim = mod(wlim,360);
elim = mod(elim,360);
[wlim,index] = sort(wlim);
elim = elim(index);
elim(elim == 0) = 360;

% Pre-allocate vectors to hold western and eastern gap limits.  These
% vectors are guaranteed to be large enough, but may have excess capacity.
wgap = zeros(size(wlim));
egap = wgap;

% Prepare to sweep through the sorted part limts by choosing an initial
% "sweep point" longitude. Being aware that the first part might or might
% not terminate a gap that wraps at lon == 0: If any parts wrap around at
% lon == 0, find the most eastern limit of such parts. Otherwise, take the
% most eastern limit overall and subtract 360.
if any(elim < wlim)
    sweep = max(elim(elim < wlim));
else
    sweep = max(elim) - 360;
end

% Sweep through the sorted part limits, looking for gaps in longitude.
n = 0;
k = 1;
numParts = numel(wlim);
while k <= numParts && sweep < 360
    if wlim(k) > sweep
        % The k-th part starts past the sweep point, which means that a
        % gap has been identified.
        n = n + 1;
        wgap(n) = sweep;
        egap(n) = wlim(k);
    end
    
    % Move the sweep point to the end of the k-th part, unless it
    % wraps around at lon == 0.
    if elim(k) >= wlim(k)
        sweep = elim(k);
    else
        % The k-th part wraps at lon == 0, so there are no more gaps to
        % be found. Set the sweep point to 360 to terminate the loop.
        sweep = 360;
    end
    k = k + 1;
end

% Remove unused elements.
wgap(n+1:end) = [];
egap(n+1:end) = [];

if isempty(wgap)
    % There are no gaps. Together, the parts span a full 360 degrees.
    lonlim = [-180 180];
else
    % Find the widest gap and take its complement: Assign the _eastern_
    % limit to lonlim(1) and the _western_ limit to lonlim(2).
    gapwidth = mod(egap - wgap, 360);
    k = find(gapwidth == max(gapwidth), 1);
    lonlim = wrapTo180([egap(k) wgap(k)]);
end
