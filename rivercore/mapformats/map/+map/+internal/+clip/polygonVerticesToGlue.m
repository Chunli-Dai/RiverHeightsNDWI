function [efirst, elast, eskip] = polygonVerticesToGlue(x, y, xlimits)
%[efirst, elast, eskip] = polygonVerticesToGlue(x, y, xlimits)
%
%   Prepare to glue polygons (x,y) on a cylinder such that xlimts(2) wraps
%   around to xlimits(1). The gluing should then proceed by (a) connecting
%   the vertices in pairs: x(elast(k)), y(elast(k)) to x(efirst(k)),
%   y(efirst(k)), keeping only one vertex per pair when y(elast(k)) ==
%   y(efirst(k)), and (b) skipping the vertices whose indices are listed in
%   eskip.

%   Copyright 2016 The MathWorks, Inc.

    % Find edge vertices.
    edgesWest = find(x == xlimits(1));
    edgesEast = find(x == xlimits(2));
    if isempty(edgesWest) || isempty(edgesEast)
        efirst = [];
        elast = [];
        eskip = [];
    else
        % Find indices and y-coordinates for vertical segments along the
        % western and eastern limits, sort in order of increasing y.
        [westseg1, westseg2, eastseg1, eastseg2, yWest1, yWest2, ...
            yEast1, yEast2] = findAndSortEdgeSegments(edgesWest, edgesEast, y);
        
        % Analyze overlap between segments along the western edge and segments
        % along the eastern edge.
        [overlapWest, overlapEast] = analyzeOverlap(yWest1, yWest2, yEast1, yEast2);
        
        % Get the endpoint indices for the segments which overlap (at more
        % than a single vertex).
        west1 = westseg1(overlapWest);
        west2 = westseg2(overlapWest);
        
        east1 = eastseg1(overlapEast);
        east2 = eastseg2(overlapEast);
        
        % When segments that meet end-to-end along either edge are being
        % removed due to overlap, skip over the vertices that connect them.
        eskip = [intersect(west1,west2); intersect(east1,east2)];
        
        % Organize connections from west-to-east and east-to-west.
        [elast, efirst] = organizeConnections(west1, west2, east1, east2, y);
                
        % Sort in order of elast
        [elast, ilast] = sort(elast);
        efirst = efirst(ilast);
    end
end


function [westseg1, westseg2, eastseg1, eastseg2, yWest1, yWest2,  ...
    yEast1, yEast2] = findAndSortEdgeSegments(edgesWest, edgesEast, y)

    % Find vertical segments between pairs of adjacent vertices where both
    % vertices fall at the same limit.
    n = length(y);
    [first, last] = internal.map.findFirstLastNonNan(y);
    [westseg1, westseg2] = findEdgeSegments(edgesWest, n, first, last);
    [eastseg1, eastseg2] = findEdgeSegments(edgesEast, n, first, last);
    
    yWest1 = y(westseg1);
    yEast1 = y(eastseg1);
    
    [yWest1, iWest] = sort(yWest1);
    [yEast1, iEast] = sort(yEast1);
    
    westseg1 = westseg1(iWest);
    westseg2 = westseg2(iWest);
    
    eastseg1 = eastseg1(iEast);
    eastseg2 = eastseg2(iEast);
    
    yWest2 = y(westseg2);
    yEast2 = y(eastseg2);
    
    % Sanity check.
    assert(all(yWest1 <= yWest2), ...
        'map:gluingPolygon:expectedDecreasingY', ...
        'Assert failed: Expected increasing y in segments along western edge.')
    
    assert(all(yEast1 >= yEast2), ...
        'map:gluingPolygon:expectedIncreasingY', ...
        'Assert failed: Expected increasing y in segments along eastern edge.')
end


function [seg1, seg2] = findEdgeSegments(edges, n, first, last)
% edge - List of edge vertices
% n - Total number of vertices
% first - Index to first element of each NaN-delimited part
% last - Index to last element of each NaN-delimited part
% seg1 - First index for each segment
% seg2 - Second (and last) index for each segment

    % Most edge segments correspond to adjacent pairs in the list of edge
    % indices.
    adjacent = find(diff(edges) == 1);
    seg1 = edges(adjacent);
    seg2 = edges(adjacent + 1);
     
    % Others may correspond to the last/first vertices of parts that begin
    % and end on an edge.
    if n > 2
        for k = 1:numel(first)
            if any(first(k) == edges) && any(last(k) == edges)
                seg1 = [seg1; last(k)];  %#ok<AGROW>
                seg2 = [seg2; first(k)]; %#ok<AGROW>
            end
        end
    end
end


function [overlapWest, overlapEast] = analyzeOverlap(yWest1, yWest2, yEast1, yEast2)
% Given the y valaues for two sets of intervals (yWest1(w), yWest2(w)) and
% (yEast1(e), yEast2(e)) for integer indices w and e, return logical
% indices indicating which intervals on the have a non-trivial overlap with
% intervals from the east, and vice-versa.

    overlapWest = false(size(yWest1));
    overlapEast = false(size(yEast1));
    
    w = 1;
    e = 1;
    while (w <= length(yWest1)) && (e <= length(yEast1))
         if yWest2(w) <= yEast2(e)
             % No overlap, or overlap only at vertex 2
             w = w + 1;
         elseif yEast1(e) <= yWest1(w)
             % No overlap, or overlap only at vertex 1
             e = e + 1;
         else
             % Overlap
             overlapWest(w) = true;
             overlapEast(e) = true;
             
             if yWest2(w) < yEast1(e)
                 w = w + 1;
             elseif yEast1(e) < yWest2(w)
                 e = e + 1;
             else
                 % yEast1(e) == yWest2(w)
                 % (includes perfect overlap as a special case)
                 w = w + 1;
                 e = e + 1;
             end
         end
    end
end


function [elast, efirst] = organizeConnections(west1, west2, east1, east2, y)
% Oganize the indices of vertices that remain (from the overlapping edge
% segments) into sets efirst and elast. When sorted by y, these should pair
% up such that the vertex with index elast(k) needs to be connected to the
% vertex with index efirst(k).
    elast  = [setdiff(west1,west2); setdiff(east1,east2)];
    efirst = [setdiff(west2,west1); setdiff(east2,east1)];
    [~, iLast]  = sort(y(elast));
    [~, iFirst] = sort(y(efirst));
    elast = elast(iLast);
    efirst = efirst(iFirst);
end
