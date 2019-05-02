function [out] = GeoToImage(varargin)

geo = varargin{1};
boundary = varargin{2};
gridspace = varargin{3};

out{1} = (geo{1}   - boundary(1))/gridspace;
out{2} = (boundary(4) - geo{2})  /gridspace;