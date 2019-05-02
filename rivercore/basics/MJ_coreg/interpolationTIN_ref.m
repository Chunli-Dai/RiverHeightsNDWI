function [Z Pa Pb Pc Pd hroh NH index_select] = interpolationTIN_ref(varargin)

Grid = varargin{1};

Gridspace = varargin{2};
%Normalized Coord
Pos1 = varargin{3};
Pos2 = varargin{4};
%Normalized Mean
Mean = varargin{5};
%Normalized Scale
Scale = varargin{6};
%DSM boundary
Boundary = varargin{7};
Size = varargin{8};
cal = varargin{9};

total_size = size(Pos1);

geo_X = Pos1*Scale(1) + Mean(1);
geo_Y = Pos2*Scale(2) + Mean(2);

col_float = (geo_X- Boundary(1))/Gridspace;
row_float = (Boundary(4)-geo_Y)/Gridspace;

col = floor(col_float);
row = floor(row_float);

X_diff = col_float - col;
Y_diff = row_float - row;

index_select = false(total_size(1),1);

index_select = (col > 2 & col < Size(2)-2 & row > 2 & row < Size(1)-2);


col_float = col_float(index_select);
row_float = row_float(index_select);
col = col(index_select);
row = row(index_select);



clear geo_X geo_Y col_float row_float ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%normalized height setting
Patch11 = ones(total_size(1),1)*(-9999);
Patch12 = ones(total_size(1),1)*(-9999);
Patch21 = ones(total_size(1),1)*(-9999);
Patch22 = ones(total_size(1),1)*(-9999);
if cal == 1
    Patch31 = ones(total_size(1),1)*(-9999);
    Patch32 = ones(total_size(1),1)*(-9999);
    Patch41 = ones(total_size(1),1)*(-9999);
    Patch33 = ones(total_size(1),1)*(-9999);
    Patch51 = ones(total_size(1),1)*(-9999);
    hroh    = zeros(total_size(1),1);
end

Patch11(index_select) = (double(Grid(row  +(col-1)*Size(1)  )) - Mean(3))/Scale(3);
Patch12(index_select) = (double(Grid(row  +(col  )*Size(1)  )) - Mean(3))/Scale(3);
Patch21(index_select) = (double(Grid(row+1+(col-1)*Size(1)  )) - Mean(3))/Scale(3);
Patch22(index_select) = (double(Grid(row+1+(col  )*Size(1)  )) - Mean(3))/Scale(3);

if cal == 1
    Patch31(index_select) = (double(Grid(row-1+(col-2)*Size(1)  )) - Mean(3))/Scale(3);
    Patch32(index_select) = (double(Grid(row-1+(col-1)*Size(1)  )) - Mean(3))/Scale(3);
    Patch41(index_select) = (double(Grid(row  +(col-2)*Size(1)  )) - Mean(3))/Scale(3);

    Patch33(index_select) = (double(Grid(row-1+(col  )*Size(1)  )) - Mean(3))/Scale(3);
    Patch51(index_select) = (double(Grid(row+1+(col-2)*Size(1)  )) - Mean(3))/Scale(3);
    
    NH{1} = Patch11;
    NH{2} = Patch12;
    NH{3} = Patch21;
    NH{4} = Patch22;
    NH{5} = Patch31;
    NH{6} = Patch32;
    NH{7} = Patch41;
    NH{8} = Patch33;
    NH{9} = Patch51;
    al = [Patch11 Patch12 Patch21 Patch22 Patch31 Patch32 Patch41 Patch33 Patch51]';
    hroh(:,1) = std(al)'.*Scale(3);
else
    NH    = zeros(1,1);
    hroh  = zeros(1,1);
end

clear Grid al;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%coordinate setting
P1 = zeros(total_size(1),3);
P2 = zeros(total_size(1),3);
P3 = zeros(total_size(1),3);
P4 = zeros(total_size(1),3);

P1(index_select,1) =(    (col  )*Gridspace + Boundary(1)   - Mean(1))/Scale(1);
P1(index_select,2) =( - ((row  )*Gridspace - Boundary(4))  - Mean(2))/Scale(2);
P1(index_select,3) = Patch11(index_select);
 
P2(index_select,1) =(    (col+1)*Gridspace + Boundary(1)   - Mean(1))/Scale(1);
P2(index_select,2) =( - ((row  )*Gridspace - Boundary(4))  - Mean(2))/Scale(2);
P2(index_select,3) = Patch12(index_select);
 
P3(index_select,1) =(    (col+1)*Gridspace + Boundary(1)   - Mean(1))/Scale(1);
P3(index_select,2) =( - ((row+1)*Gridspace - Boundary(4))  - Mean(2))/Scale(2);
P3(index_select,3) = Patch22(index_select);
 
P4(index_select,1) =(    (col  )*Gridspace + Boundary(1)  - Mean(1))/Scale(1);
P4(index_select,2) =( - ((row+1)*Gridspace - Boundary(4)) - Mean(2))/Scale(2);
P4(index_select,3) = Patch21(index_select);

clear Patch11 Patch12 Patch22 Patch21 Patch31 Patch32 Patch41 Patch33 Patch51 col row;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%vector calculation
v1 = zeros(total_size(1),3);
v2 = zeros(total_size(1),3);
n1 = zeros(total_size(1),3);

v1 = P3-P1;
v2 = P2-P1;
n1(:,1) =    v1(:,2).*v2(:,3) - v2(:,2).*v1(:,3) ;
n1(:,2) = - (v1(:,1).*v2(:,3) - v2(:,1).*v1(:,3));
n1(:,3) =    v1(:,1).*v2(:,2) - v2(:,1).*v1(:,2);

clear v1 v2 ;

v1 = zeros(total_size(1),3);
v2 = zeros(total_size(1),3);
n2 = zeros(total_size(1),3);

v1 = P4-P1;
v2 = P3-P1;
n2(:,1) =    v1(:,2).*v2(:,3) - v2(:,2).*v1(:,3) ;
n2(:,2) = - (v1(:,1).*v2(:,3) - v2(:,1).*v1(:,3));
n2(:,3) =    v1(:,1).*v2(:,2) - v2(:,1).*v1(:,2);

clear v1 v2 ;

clear v1 v2 mag P2 P3 P4;

%plane vector decide
% X_diff == 0 && Y_diff == 0
index_1 = find(X_diff == 0 & Y_diff == 0);
% X_diff == 0 && Y_diff ~= 0
index_2 = find(X_diff == 0 & Y_diff ~= 0);
% X_diff ~= 0 && Y_diff == 0
index_3 = find(X_diff ~= 0 & Y_diff == 0);
% X_diff > 0 && Y_diff > 0 && (X_diff <= Y_diff)
index_4 = find(X_diff > 0 & Y_diff > 0 & (X_diff <= Y_diff));
% X_diff > 0 && Y_diff > 0 && (X_diff > Y_diff)
index_5 = find(X_diff > 0 & Y_diff > 0 & (X_diff > Y_diff));

clear X_diff Y_diff;

%pt_size = size(index_select);
n = zeros(total_size(1),3);

pt_size = size(index_1);
if pt_size(1) > 0
    n(index_1,:) = n1(index_1,:);
end
pt_size = size(index_2);
if pt_size(1) > 0
    n(index_2,:) = n2(index_2,:);
end
pt_size = size(index_3);
if pt_size(1) > 0
    n(index_3,:) = n1(index_3,:);
end
pt_size = size(index_4);
if pt_size(1) > 0
    n(index_4,:) = n2(index_4,:);
end
pt_size = size(index_5);
if pt_size(1) > 0
    n(index_5,:) = n1(index_5,:);
end

clear n1 n2 index_1 index_2 index_3 index_4 index_5;

mag = zeros(total_size(1),1);
mag(index_select,1) = sqrt(n(index_select,1).*n(index_select,1) + n(index_select,2).*n(index_select,2) + n(index_select,3).*n(index_select,3));
t_size = size(mag);
if t_size(1) == 0
    mag_flag = index_select;
else
    mag_flag = (mag > 0);
end
index_select = index_select & mag_flag;
n(index_select,1) = n(index_select,1)./mag(index_select,1);
n(index_select,2) = n(index_select,2)./mag(index_select,1);
n(index_select,3) = n(index_select,3)./mag(index_select,1);

d = zeros(total_size(1),1);
d(index_select,1) = -P1(index_select,1).*n(index_select,1) - P1(index_select,2).*n(index_select,2) - P1(index_select,3).*n(index_select,3);

Z = ones(total_size(1),1)*(-9999);

tttt    = find(index_select);
tt_size = size(tttt);
if tt_size(1) > 0
    Z(index_select,1) = -((Pos1(index_select).*n(index_select,1) + Pos2(index_select).*n(index_select,2) +d(index_select,1) )./n(index_select,3));
end

index_select = Z > -1 & Z < 1;

Pa = zeros(total_size(1),1);
Pb = zeros(total_size(1),1);
Pc = zeros(total_size(1),1);
Pd = zeros(total_size(1),1);

Pa(:,1) = n(:,1);
Pb(:,1) = n(:,2);
Pc(:,1) = n(:,3);
Pd(:,1) = d;

clear mag d n Pos1 Pos2 P1;


 