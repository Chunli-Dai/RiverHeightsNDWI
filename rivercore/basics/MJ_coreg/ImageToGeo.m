function [out] = ImageToGeo(varargin)

image = varargin{1};
boundary = varargin{2};
gridspace = varargin{3};

out{1} = image{1}*gridspace + boundary(1);
out{2} = boundary(4) - image{2}*gridspace;