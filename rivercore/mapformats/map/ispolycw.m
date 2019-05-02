function tf = ispolycw(x, y)
%ISPOLYCW True if polygon vertices are in clockwise order
%
%   TF = ISPOLYCW(X, Y) returns true if the polygonal contour vertices 
%   represented by X and Y are ordered in the clockwise direction.  X and Y
%   are numeric vectors with the same number of elements.
%
%   Alternatively, X and Y can contain multiple contours, either in
%   NaN-separated vector form or in cell array form.  In that case,
%   ISPOLYCW returns a logical array containing one true or false value
%   per contour.
%
%   ISPOLYCW always returns true for polygonal contours containing two or
%   fewer vertices.
%
%   Vertex ordering is not well defined for self-intersecting polygonal
%   contours.  For such contours, ISPOLYCW returns a result based on the
%   order or vertices immediately before and after the left-most of the 
%   lowest vertices.  In other words, of the vertices with the lowest Y
%   value, find the vertex with the lowest X value.  For a few special
%   cases of self-intersecting contours, the vertex ordering cannot be
%   determined using only the left-most of the lowest vertices; for these
%   cases, ISPOLYCW uses a signed area test to determine the ordering.
%
%   Class Support
%   -------------
%   X and Y may be any numeric class.
%
%   Example
%   -------
%   Orientation of a square:
%
%       x = [0 1 1 0 0];
%       y = [0 0 1 1 0];
%       ispolycw(x, y)                     % Returns 0
%       ispolycw(fliplr(x), fliplr(y))     % Returns 1
%
%   See also POLY2CW, POLY2CCW, POLYBOOL

% Copyright 2004-2015 The MathWorks, Inc.

if isempty(x)
   tf = true;
elseif iscell(x)
    tf = false(size(x));
    for k = 1:numel(x)
        tf(k) = isContourClockwise(x{k}, y{k});
    end
else
    checkxy(x, y, mfilename, 'X', 'Y', 1, 2)
    [first, last] = internal.map.findFirstLastNonNan(x);
    numParts = numel(first);
    if isrow(x)
        tf = false(1,numParts);
    else
        tf = false(numParts,1);
    end
    for k = 1:numParts
        s = first(k);
        e = last(k);
        tf(k) = isContourClockwise(x(s:e), y(s:e));
    end
end

%----------------------------------------------------------------------
function tf = isContourClockwise(x, y)

[x, y] = removeDuplicates(x, y);
if numel(x) > 2
    tf = (signedArea(x, y) <= 0);    
else
    tf = true;
end

%----------------------------------------------------------------------
function [x, y] = removeDuplicates(x, y)
% ... including duplicate start and end points.

is_closed = (x(1) == x(end)) && (y(1) == y(end));
if is_closed
    x(end) = [];
    y(end) = [];
end

dups = [false; (diff(x(:)) == 0) & (diff(y(:)) == 0)];
x(dups) = [];
y(dups) = [];

%----------------------------------------------------------------------
function a = signedArea(x, y)
% a = signedArea(x,y) returns twice the signed area of the polygonal
% contour represented by vectors x and y.  Assumes (x,y) is NOT closed.

% Reference: 
% http://geometryalgorithms.com/Archive/algorithm_0101/algorithm_0101.htm

x = x - mean(x);
n = numel(x);
if n <= 2
    a = 0;
else
    i = [2:n 1];
    j = [3:n 1 2];
    k = (1:n);
    a = sum(x(i) .* (y(j) - y(k)));
end
