function [xg, yg] = gluePolygonOnCylinder(x, y, xlimits)
%gluePolygonOnCylinder Glue polygons such that they warp around a cylinder
%
%   [xg,yg] = gluePolygonsOnCylinder(x,y,xlimits) glues polygons (x,y) on a
%   cylinder such that xlimits(2) wraps around to xlimits(1). x and y
%   vectors of matching size and can be double or single.  The output
%   vectors, xg and yg, have the same class as x and y. xlimits is
%   typically double, but that is not enforced.

% Copyright 2016 The MathWorks, Inc.

    ex = zeros([0 0],'like',x);
    ey = zeros([0 0],'like',y);
    
    if all(isnan(x))
        x = ex;
        y = ey;
    end
    
    [x, y] = removeExtraNanSeparators(x, y);
    if isempty(x)
        xg = ex;
        yg = ey;
        return
    end
    
    % Unless the polygons extend to both limits, there's no work to do.
    if (xlimits(1) < min(x)) || (max(x) < xlimits(2))
        xg = x;
        yg = y;
        return
    end
        
    x = x(:);
    y = y(:);
    
    [efirst, elast, eskip] = map.internal.clip.polygonVerticesToGlue(x, y, xlimits);
    [xg, yg] = linkPolygonParts(x, y, xlimits, efirst, elast, eskip);
end


function [xg, yg] = linkPolygonParts(x, y, xlimits, efirst, elast, eskip)
% Trace and link the polygon parts: (a) connecting the vertices in pairs:
% x(elast(k)), y(elast(k)) to x(efirst(k)), y(efirst(k)), keeping only one
% vertex per pair when y(elast(k)) == y(efirst(k)), and (b) skipping the
% vertices whose indices are listed in eskip.

    [x, y, efirst, elast] = removeVerticesToSkip(x, y, efirst, elast, eskip);
    
    % Allocate xg and yg to hold the "glued" polygons. As a result of
    % gluing, the number of vertices should decrease, so the original size
    % of x and y should be more than enough. There are two advantages to
    % starting with NaN-filled vectors: (1) We can drop in parts without
    % having to explicitly add NaN-separators -- instead we just skip a
    % slot between subsequent parts and (2) we can easily remove unused
    % slots when done.
    xg = NaN(size(x),'like',x);
    yg = NaN(size(y),'like',y);
    
    % Merge the efirst/elast indices into first and last for the existing
    % parts, breaking up the parts into smaller pieces. Sorts are required
    % to set up the correct associations between the new first-last pairs.
    [first, last] = internal.map.findFirstLastNonNan(x);
    
%     % The next 3 lines are equivalent to:
%     %     [first, ~, next] = unique([first; efirst]);
%     [first, iFirst] = sort([first; efirst]);
%     next = (1:length(iFirst))';
%     next(iFirst) = next;
%     
%     [last, iLast] = sort([last; elast]);
%     next = next(iLast);
    
    [first, ~, next] = unique([first; efirst]);    
    [last, iLast] = unique([last; elast]);
    next = next(iLast);
    
    % Keep track of the number of pieces traced in relation to the total
    % number, to ensure that the while loop terminates (see assert below).
    numPieces = numel(first);
    numTraced = 0;
    
    n = 1;
    traced = false(numPieces,1);
    for k = 1:numPieces
        if ~traced(k)
            kFirst = first(k);
            lastxg = NaN('like',x);
            lastyg = NaN('like',y);
            
            % Copy the rest of the k-th piece and the pieces to which
            % it connects.
            j = k;
            while ~traced(j)
                % Ensure that the loop cannot run forever. The following
                % assertion should never be triggered.
                numTraced = numTraced + 1;
                assert(numTraced <= numPieces, ...
                    'map:gluePolygonsOnVerticalEdges:tracingFailed', ...
                    'Failed to converge when tracing open curves.')
                
                jFirst = first(j);
                bothOnEdge = any(x(jFirst) == xlimits) && any(lastxg == xlimits);
                firstIsDuplicate = (y(jFirst) == lastyg) ...
                    && (x(jFirst) == lastxg || bothOnEdge);
                if firstIsDuplicate
                    % Skip over the first vertex in the j-th piece, copying
                    % vertices 2 through end.
                    s = 1 + jFirst;
                else
                    % Copy the entire j-th piece.
                    s = jFirst;
                end
                
                e = last(j);
                m = n + e - s;
                xg(n:m) = x(s:e);
                yg(n:m) = y(s:e);
                
                % Set up for next iteration
                traced(j) = true;
                j = next(j);
                lastxg = xg(m);
                lastyg = yg(m);
                n = m + 1;
            end
            
            % If necessary, duplicate a vertex to close the loop.
            bothOnEdge = any(x(kFirst) == xlimits) && any(lastxg == xlimits);
            if (lastyg ~= y(kFirst)) || (lastxg ~= x(kFirst) && ~bothOnEdge)
                xg(n) = x(kFirst);
                yg(n) = y(kFirst);                
                n = n + 1;
            end            
            
            % Leave a NaN-separator in xg and yg to separate this curve
            % from the next one.
            xg(n) = NaN;
            yg(n) = NaN;
            n = n + 1;
        end
    end
    
    % Remove trailing NaNs
    k = length(xg);
    while isnan(xg(k))
        xg(k) = [];
        yg(k) = [];
        k = k - 1;
    end
end


function [x, y, efirst, elast] = removeVerticesToSkip(x, y, efirst, elast, eskip)
% Skip the vertices whose indices are listed in eskip:  Reduce the value of
% each element in efirst and elast by one for each element of eskip that
% precedes it. Remove the elements listed in eskip from x and y.

    s = zeros(size(x));
    s(eskip) = 1;
    skipcount = cumsum(s);
    
    efirst = efirst - skipcount(efirst);
    elast = elast - skipcount(elast);
    
    x(eskip) = [];
    y(eskip) = [];
end
