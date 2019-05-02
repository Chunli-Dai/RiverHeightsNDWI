function [out]=savefromindex(varargin)

array_size = varargin{1};
index = varargin{2};
value = varargin{3};
null_value = varargin{4};

save_array_t = ones(array_size(1)*array_size(2),1)*null_value;
save_array_t(index) =  value;
out = reshape(save_array_t,array_size(1),array_size(2));
