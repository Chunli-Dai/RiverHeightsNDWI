function [Boundary Refboundary Tarboundary RefHeight_min_max TarHeight_min_max] = BoundarySetting(varargin)

RefImg              = varargin{1};
TarImg              = varargin{2};
Null_height         = varargin{3};
Buffer_grid         = varargin{4};
RefImg_gridspace    = varargin{5};

Refboundary(1)      = min(RefImg.x);
Refboundary(2)      = min(RefImg.y);
Refboundary(3)      = max(RefImg.x);
Refboundary(4)      = max(RefImg.y);
RefH                = RefImg.z(RefImg.z>0.0 & RefImg.z < Null_height);
t_size              = size(RefH);
if t_size(2) > 10
    RefHeight_min_max(1)= min(RefH);
    RefHeight_min_max(2)= max(RefH);
else
    RefHeight_min_max(1)= 0;
    RefHeight_min_max(2)= Null_height;
end


Tarboundary(1)      = min(TarImg.x);
Tarboundary(2)      = min(TarImg.y);
Tarboundary(3)      = max(TarImg.x);
Tarboundary(4)      = max(TarImg.y);
TarH                = TarImg.z(TarImg.z>0.0 & TarImg.z < Null_height);
t_size              = size(TarH);
if t_size(2) > 10
    TarHeight_min_max(1)= min(TarH);
    TarHeight_min_max(2)= max(TarH);
else
    TarHeight_min_max(1)= RefHeight_min_max(1);
    TarHeight_min_max(2)= RefHeight_min_max(2);
end

buffer              = double(Buffer_grid*RefImg_gridspace);
Boundary            = zeros(1,4);
Boundary(1)         = max(Refboundary(1), Tarboundary(1)) + buffer;
Boundary(2)         = max(Refboundary(2), Tarboundary(2)) + buffer;
Boundary(3)         = min(Refboundary(3), Tarboundary(3)) - buffer;
Boundary(4)         = min(Refboundary(4), Tarboundary(4)) - buffer;

for i=1:2
    if Tarboundary(i) > Refboundary(i)
        Boundary(i) = Refboundary(i) + (floor((Tarboundary(i) - Refboundary(i))/RefImg_gridspace)+1)*RefImg_gridspace;
    end
end

for i=3:4
    if Tarboundary(i) < Refboundary(i)
        Boundary(i) = Refboundary(i) - (floor((Refboundary(i) - Tarboundary(i))/RefImg_gridspace)+1)*RefImg_gridspace;
    end
end

clear RefH TarH RefImg TarImg ;