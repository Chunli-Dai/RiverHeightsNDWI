function [out] = ConformalTransform(varargin)

coord = varargin{1};
X = varargin{2};

R = RotationMatrix([X(2,1) X(3,1) X(4,1)]);

out{1} = X(1,1)*(R(1,1)*coord{1} + R(2,1)*coord{2} + R(3,1)*coord{3}) + X(5,1);
out{2} = X(1,1)*(R(1,2)*coord{1} + R(2,2)*coord{2} + R(3,2)*coord{3}) + X(6,1);
out{3} = X(1,1)*(R(1,3)*coord{1} + R(2,3)*coord{2} + R(3,3)*coord{3}) + X(7,1);