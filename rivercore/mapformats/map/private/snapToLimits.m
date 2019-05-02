function [x,y] = snapToLimits(x, y, xLimit, tolSnap)
% For each open curve in the NaN-delimited arrays (x,y), if the first or
% last vertex is within distance tolSnap of one of the limits defined
% by xLimit, snap it to the limit.

% Copyright 2016 The MathWorks, Inc.

% Adapted from subfunction in private/closePolygonInRectangle.

% Indices of the first and last vertex in each curve.
[first, last] = internal.map.findFirstLastNonNan(x);

% Identify open curves.
isOpen = ~((x(first) == x(last)) & (y(first) == y(last)));

% Snap first vertices that belong to open
% curves and are close to the limits.
x(first(isOpen & (abs(x(first) - xLimit(1)) < tolSnap))) = xLimit(1);
x(first(isOpen & (abs(x(first) - xLimit(2)) < tolSnap))) = xLimit(2);

% Snap last vertices that belong to open
% curves and are close to the limits.
x(last(isOpen & (abs(x(last) - xLimit(1)) < tolSnap))) = xLimit(1);
x(last(isOpen & (abs(x(last) - xLimit(2)) < tolSnap))) = xLimit(2);
