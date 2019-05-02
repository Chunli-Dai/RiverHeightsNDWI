function [out] = FindImageFromCoord(varargin)

img = varargin{1};
boundary = varargin{2};
gridspace = varargin{3};
coord = varargin{4};

image_coord = GeoToImage(coord,boundary,gridspace);
image_coord{1} = floor(image_coord{1});
image_coord{2} = floor(image_coord{2});
[size_row size_col] = size(img);
temp_array = img(:);

null_row  = (image_coord{2} > size_row | image_coord{2} < 1);
null_col  = (image_coord{1} > size_col | image_coord{1} < 1);
index_null= null_row | null_col;
index     = index_null == 0;

t_size_x  = size(image_coord{1});
t_size_y  = size(image_coord{2});

out(index) = temp_array((image_coord{1}(index)-1)*size_row + image_coord{2}(index));
out(index_null) = 1; 

out        = out.';
% a = find(index_null);
% b = find(index);