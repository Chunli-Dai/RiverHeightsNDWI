function [S A]=SlopeAspect(varargin)

a = varargin{1};
b = varargin{2};
c = varargin{3};
scale = varargin{4};

a = a./(scale(1));
b = b./(scale(2));
c = c./(scale(3));

denominator = sqrt(a.*a + b.*b + c.*c);
value = c./denominator;

S = acos(value)*(180/pi);

t_size = size(a);
aspect = zeros(t_size(1),1);
index_2 = find(a ~= 0);
value  = b(index_2)./a(index_2);
aspect(index_2,1) = 90 - atan2(b(index_2),a(index_2))*(180/pi);
index_1 = find(aspect < 0);
aspect(index_1) = aspect(index_1) + 360;

A = aspect;