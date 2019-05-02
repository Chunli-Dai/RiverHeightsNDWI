function [out] = Normalize_coord(varargin)

coord = varargin{1};
Mean = varargin{2};
Scale = varargin{3};
Z_only_check = varargin{4};

if strcmp(Z_only_check,'Z');
    out    = (coord{3} - Mean(3))/Scale(3);
elseif strcmp(Z_only_check,'XY');
    for i=1:2
        out{i}    = (coord{i} - Mean(i))/Scale(i);
    end
else
    for i=1:3
        out{i}    = (coord{i} - Mean(i))/Scale(i);
    end
end

