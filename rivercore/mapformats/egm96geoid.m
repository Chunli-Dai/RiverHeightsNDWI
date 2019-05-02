function [N, refvec] = egm96geoid(varargin)
% EGM96GEOID Read 15-minute gridded geoid heights from EGM96
%
%   [N, REFVEC] = EGM96GEOID(SAMPLEFACTOR) imports global geoid height in
%   meters from the EGM96 geoid model. The data set is gridded at 15-minute
%   intervals, but may be down sampled as specified by the positive integer
%   SAMPLEFACTOR. The result is returned in the regular data grid N along
%   with referencing vector REFVEC. At full resolution (SAMPLEFACTOR = 1),
%   N will be 721-by-1441. The data grid has a raster interpretation of
%   'postings'.
%
%   The gridded EGM96 data set must be on your path in a file named
%   'WW15MGH.GRD'.
% 
%   [N, REFVEC] = EGM96GEOID(SAMPLEFACTOR, LATLIM, LONLIM) imports data for
%   the part of the world within the specified latitude and longitude
%   limits. The limits must be two-element vectors in units of degrees.
%   Longitude limits can be defined in the range [-180 180] or [0 360]. For
%   example, lonlim = [170 190] returns data centered on the dateline, while
%   lonlim = [-10 10] returns data centered on the prime meridian.
%
%   For details on locating map data for download over the Internet, see
%   the following documentation at the MathWorks web site:
%
%   <a href="matlab:
%   web('https://www.mathworks.com/help/map/finding-geospatial-data.html')
%   ">https://www.mathworks.com/help/map/finding-geospatial-data.html</a>

% Copyright 1996-2016 The MathWorks, Inc.

    narginchk(1,3)
    if nargin == 2
        error(message('map:validate:invalidArgCount'))
    end

    [samplefactor, latlim, lonlim] = parseInputs(varargin);

    % Import geoid height grid N and referencing object R.
    [N, R] = readEGM96(latlim, lonlim, samplefactor);

    if isempty(N)
        refvec = [];
    else
        % Convert referencing object to referencing vector.
        refvec = geopostings2refvec(R);
    end
end


function [samplefactor, latlim, lonlim] = parseInputs(inputs)

    numInputs = numel(inputs);
    try
        % Validate SAMPLEFACTOR.
        samplefactor = inputs{1};
        maxSampleFactor = 4*180;
        validateattributes(samplefactor, {'numeric'}, ...
            {'positive','scalar','integer','<=', maxSampleFactor}, ...
            mfilename, 'SAMPLEFACTOR', 1)

        % Set or validate LATLIM.
        if numInputs < 2
            latlim = [-90 90];
        else
            latlim = inputs{2};

            if isempty(latlim)
                latlim = [-90 90];
            end

            validateattributes(latlim, {'double'}, ...
                {'real','vector','finite','>=',-90,'<=',90}, mfilename, 'LATLIM')

            map.internal.assert(numel(latlim) == 2, ...
                'map:validate:expectedTwoElementVector', 'LATLIM');

            latlim = sort(latlim);

            map.internal.assert(latlim(1) < latlim(2), ...
                'map:maplimits:expectedAscendingLatlim')

            latlim = latlim(:)';
        end

        % Set or validate LONLIM.
        if numInputs < 3
            lonlim = [0 360];
        else
            lonlim = inputs{3};

            if isempty(lonlim)
                lonlim = [0 360];
            end

            validateattributes(lonlim, {'double'}, ...
                {'real','vector','finite'}, mfilename, 'LONLIM')

            map.internal.assert(numel(lonlim) == 2, ...
                'map:validate:expectedTwoElementVector', 'LONLIM');

            map.internal.assert(all(lonlim <= 360), ...
                'map:validate:expectedRange', 'LONLIM',  '0', 'lonlim', '360')

            if any(lonlim < 0) && (any(lonlim < -180) || any(lonlim > 180))
                error(message('map:validate:expectedRange', ...
                    'LONLIM', '-180', 'lonlim', '180'))
            end

            if lonlim(2) < lonlim(1)
                lonlim(2) = lonlim(2) + 360;
            end

            lonlim = lonlim(:)';
        end
    catch me
        throwAsCaller(me)
    end
end


function [N, R] = readEGM96(latlim, lonlim, samplefactor)
% Read, crop, and subsample the full geoid height grid, returning the
% output grid and a geographic postings raster reference object.

    fileID = openFile();
    if fileID == -1
        % User has canceled file dialog.
        N = [];
        R = [];
    else
        validateFileSize(fileID)
        clean = onCleanup(@() fclose(fileID));
        globalN = readGlobalGrid(fileID);
        
        glatlim = [-90 90];
        glonlim = [0 360];
        spacing = 0.25;
        globalR = georefpostings(glatlim, glonlim, spacing, spacing);
        
        [R, rows, cols] = setupCropAndSubsampleForGlobalPostings( ...
            globalR, latlim, lonlim, samplefactor);
        
        N = globalN(rows, cols);
    end
end


function validateFileSize(fileID)
% Compare to expected sizes for Unix and Windows.

    expectedSizeUnix = 9618935;  % One line-ending character
    expectedSizeWin  = 9756647;  % Two line-ending characters

    fseek(fileID,0,'eof');
    fileSize = ftell(fileID);
    if fileSize ~= expectedSizeUnix && fileSize ~= expectedSizeWin
        error(message('map:egm96geoid:invalidFileSize'))
    end
    fseek(fileID,0,'bof');
end


function fileID = openFile()
% Try to open the EGM96 geoid file.

    filename = 'WW15MGH.GRD';
    fileID = fopen(filename,'r');

    if fileID == -1
        filename = lower(filename);
        fileID = fopen(filename,'r');
        if fileID == -1
            [filename,filepath] = uigetfile(filename,['Where is ',filename,'?']);
            if filename == 0
                fileID = -1;
            else
                fileID = fopen([filepath,filename],'r');
                if fileID == -1
                    error(message('map:fileio:fileNotFound',[filepath,filename]))
                end
            end
        end
    end
end


function N = readGlobalGrid(fileID)
% Read the full data grid, after reading and validating the header line.

% File layout from NGS documentation:
%
% This file contains 1038961 point values in grid form.  The first row of the file
% is the "header" of the file and shows the south, north, west, and east limits of
% the file followed by the grid spacing in n-s and e-w. All values in the "header"
% are in DECIMAL DEGREES.
% 
% The geoid undulation grid is computed at 15 arc minute spacings in north/south
% and east/west with the new "EGM96" spherical harmonic potential coefficient set
% complete to degree and order 360 and a geoid height correction value computed
% from a set of spherical harmonic coefficients ("CORRCOEF"), also to degree and
% order 360.  The file is arranged from north to south, west to east (i.e., the
% data after the header is the north most latitude band and is ordered from west
% to east).
% 
% The coverage of this file is:
% 
%                90.00 N  +------------------+
%                         |                  |
%                         | 15' spacing N/S  |
%                         |                  |
%                         |                  |
%                         | 15' spacing E/W  |
%                         |                  |
%               -90.00 N  +------------------+
%                        0.00 E           360.00 E

    % Skip header line
    fgetl(fileID);
    
    nrows =  721;
    ncols = 1441;

    N = textscan(fileID,'%f','EndOfLine','\r\n','CollectOutput',true);
    N = N{1};
    N = reshape(N,[ncols nrows]);
    N = flipud(N');
end


function [R, rows, cols] = setupCropAndSubsampleForGlobalPostings( ...
    globalR, latlim, lonlim, samplefactor)

% This function can be very simple for two reasons: (1) the input globalR
% is assumed to correspond to a raster of global extent, and (2) both the
% input and ouput rasters are snapped to multiples of the input sample
% spacings.  These conditions guarantee that every output sample exactly
% coincides (in latitude and longitude) with an input sample.  The
% function also assumes that globalR.ColumnsStartFrom is 'south' and
% globalR.RowsStartFrom is 'west.'
%
% See also map.rasterref.GeographicRasterReference/setupCropAndSubsample

    % Pre-snap the limits because, in the special cases in which
    % diff(latlim) or diff(lonlim) is an exact multiple of latspace or
    % lonspace, georefpostings won't snap the limits for us.
    [latlim, lonlim] = snapLimitsToRaster(globalR, latlim, lonlim);

    latspace = samplefactor * globalR.SampleSpacingInLatitude;
    lonspace = samplefactor * globalR.SampleSpacingInLongitude;
    
    R = georefpostings(latlim, lonlim, latspace, lonspace);

    rows = 1:R.RasterSize(1);
    cols = 1:R.RasterSize(2);

    if ~isequal(R, globalR)
        lat = intrinsicYToLatitude(R, rows);
        lon = intrinsicXToLongitude(R, cols);

        rows = round(latitudeToIntrinsicY(globalR, lat));
        cols = round(longitudeToIntrinsicX(globalR, lon));
    end
end


function [latlim, lonlim] = snapLimitsToRaster(R, latlim, lonlim)
% Assume that R.ColumnsStartFrom is 'south', R.RowsStartFrom is 'west',
% R.LongitudeLimits(1) = 0, R.SampleSpacingInLongitude = 1/4:
    
    ylimi = latitudeToIntrinsicY(R, latlim);
    ylimi = [floor(ylimi(1)) ceil(ylimi(2))];
    latlim = intrinsicYToLatitude(R, ylimi);
    
    % Need to scale and shift without wrapping:
    % avoid the longitudeToIntrinsicX method.
    xlimi = 1 + 4 * lonlim;
    xlimi = [floor(xlimi(1)) ceil(xlimi(2))];
    lonlim = intrinsicXToLongitude(R, xlimi);
end


function refvec = geopostings2refvec(R)
%   This function assumes:
%
%       'postings'
%       columns start from west
%       rows start from south
%       sample density is the same in both dimensions
%
%   This is not a general purpose conversion. Instead, the limits are
%   extrapolated 1/2 sample north and west of their true locations, for
%   consistency with the previous behavior of egm96geoid.

    cellsPerDegree = sampleDensity(R);

    % Step north 1/2 cell/sample-spacing relative to northern limit.
    northernLimit = R.LatitudeLimits(2) + 0.5 / cellsPerDegree;

    % Step west 1/2 cell/sample-spacing relative to the western limit.
    westernLimit  = R.LongitudeLimits(1) - 0.5 / cellsPerDegree;

    refvec = [cellsPerDegree northernLimit westernLimit];
end
