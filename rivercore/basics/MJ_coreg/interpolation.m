function [Z flag]=interpolation(varargin)

Grid = varargin{1};
Pos_x = varargin{2};
Pos_y = varargin{3};
r = varargin{4};
c = varargin{5};
t_size = size(Pos_x);
X_diff = zeros(t_size(1),1);
Y_diff = zeros(t_size(1),1);
re_grid1 = ones(t_size(1),1)*(-9999);
re_grid2 = ones(t_size(1),1)*(-9999);
re_grid3 = ones(t_size(1),1)*(-9999);
re_grid4 = ones(t_size(1),1)*(-9999);
index11  = zeros(t_size(1),1);
index12  = zeros(t_size(1),1);
index21  = zeros(t_size(1),1);
index22  = zeros(t_size(1),1);

col = floor(Pos_x);
row = floor(Pos_y);

flag = col > 1 & col < c-1 & row > 1 & row < r-1;


X_diff(flag) = Pos_x(flag) - col(flag);
Y_diff(flag) = Pos_y(flag) - row(flag);

index11(flag) = row(flag)     + (col(flag)-1)*r;
index12(flag) = row(flag)     + (col(flag)  )*r;
index21(flag) = (row(flag)+1) + (col(flag)-1)*r;
index22(flag) = (row(flag)+1) + (col(flag)  )*r;

re_grid1(flag) = double(Grid(index11(flag)));
re_grid2(flag) = double(Grid(index12(flag)));
re_grid3(flag) = double(Grid(index21(flag)));
re_grid4(flag) = double(Grid(index22(flag)));

nei_flag = re_grid1 > 0 & re_grid2 > 0 & re_grid3 > 0 & re_grid4 > 0;
flag = flag & nei_flag;

Z = ones(t_size(1),1)*(-9999);
Z(flag) = (X_diff(flag).*Y_diff(flag)).*re_grid4(flag) + ((1.0-X_diff(flag)).*Y_diff(flag)).*re_grid3(flag) + (X_diff(flag).*(1-Y_diff(flag))).*re_grid2(flag) + ((1.0-X_diff(flag)).*(1-Y_diff(flag))).*re_grid1(flag);



clear Grid Pos col row X_diff Y_diff Patch;


