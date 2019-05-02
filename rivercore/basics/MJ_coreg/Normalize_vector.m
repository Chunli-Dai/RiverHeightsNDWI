function [out]=Normalize_vector(varargin)

vector = varargin{1};
scale = varargin{2};
pt = varargin{3};

for i=1:3
   norvector{i} = vector{i}*scale(i);
end

out{4} = -norvector{1}.*pt{1} - norvector{2}.*pt{2} - norvector{3}.*pt{3};