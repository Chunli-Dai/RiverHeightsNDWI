function [Check Z Coeff] = interpolationTIN(varargin)

Grid = varargin{1};

Gridspace = varargin{2};
%Normalized Coord
Pos = varargin{3};
%Normalized Mean
Mean = varargin{4};
%Normalized Scale
Scale = varargin{5};
%DSM boundary
Boundary = varargin{6};
%Z_shift
shift = varargin{7};
imgchk = varargin{8};
X = varargin{9};

Check = 'False';

Size = size(Grid);

geo_X = Pos(1)*Scale(1) + Mean(1);
geo_Y = Pos(2)*Scale(2) + Mean(2);

col_float = (geo_X- Boundary(1))/Gridspace +1;
row_float = (Boundary(4)-geo_Y)/Gridspace +1;

col = floor(col_float);
row = floor(row_float);

X_diff = col_float - col;
Y_diff = row_float - row;
P1 = zeros(3,1);
P2 = zeros(3,1);
P3 = zeros(3,1);
P4 = zeros(3,1);
P5 = zeros(3,1);
P6 = zeros(3,1);
P7 = zeros(3,1);

Coeff = zeros(4,1);

if row+1 < Size(1) && col+1 < Size(2) && row > 1 && col > 1
    Check = 'True';
    
    Patch(1,1) = (double(Grid(row  ,col  )+shift(3)) - Mean(3))/Scale(3);
    Patch(1,2) = (double(Grid(row  ,col+1)+shift(3)) - Mean(3))/Scale(3);
    Patch(2,1) = (double(Grid(row+1,col  )+shift(3)) - Mean(3))/Scale(3);
    Patch(2,2) = (double(Grid(row+1,col+1)+shift(3)) - Mean(3))/Scale(3);
    
    %upper TIN
    Patch(3,1) = (double(Grid(row-1,col-1)+shift(3)) - Mean(3))/Scale(3);
    Patch(3,2) = (double(Grid(row-1,col  )+shift(3)) - Mean(3))/Scale(3);
    Patch(4,1) = (double(Grid(row  ,col-1)+shift(3)) - Mean(3))/Scale(3);
    
    for i=1:4
        for j=1:2
            if Patch(i,j) < -1.0 || Patch(i,j) > 1.0
                Check = 'False';
            end
        end
    end
    
    if strcmp(Check,'True')
        P1(1) = (col-1)*Gridspace + Boundary(1);
        P1(2) = - ((row-1)*Gridspace - Boundary(4));
        P1(3) = Patch(1,1);

        P2(1) = (col)*Gridspace + Boundary(1);
        P2(2) = - ((row-1)*Gridspace - Boundary(4));
        P2(3) = Patch(1,2);

        P3(1) = (col)*Gridspace + Boundary(1);
        P3(2) = - ((row)*Gridspace - Boundary(4));
        P3(3) = Patch(2,2);

        P4(1) = (col-1)*Gridspace + Boundary(1);
        P4(2) = - ((row)*Gridspace - Boundary(4));
        P4(3) = Patch(2,1);

        %upper TIN
        P5(1) = (col-2)*Gridspace + Boundary(1);
        P5(2) = - ((row-2)*Gridspace - Boundary(4));
        P5(3) = Patch(3,1);

        P6(1) = (col-1)*Gridspace + Boundary(1);
        P6(2) = - ((row-2)*Gridspace - Boundary(4));
        P6(3) = Patch(3,2);

        P7(1) = (col-2)*Gridspace + Boundary(1);
        P7(2) = - ((row-1)*Gridspace - Boundary(4));
        P7(3) = Patch(4,1);

        for num=1:2
            P1(num) = (P1(num) - Mean(num))/Scale(num);
            P2(num) = (P2(num) - Mean(num))/Scale(num);
            P3(num) = (P3(num) - Mean(num))/Scale(num);
            P4(num) = (P4(num) - Mean(num))/Scale(num);

            P5(num) = (P5(num) - Mean(num))/Scale(num);
            P6(num) = (P6(num) - Mean(num))/Scale(num);
            P7(num) = (P7(num) - Mean(num))/Scale(num);
        end

        if strcmp(imgchk,'Tar')
            R = RotationMatrix([X(2,1) X(3,1) X(4,1)]);
            %R = RotationMatrix([0.0 0.0 0.0]);
            for num = 3:3;

            P1(num) = X(1,1)*( R(1,num)*P1(1) + R(2,num)*P1(2) + R(3,num)*P1(3) ) + X(num+4,1);
            P2(num) = X(1,1)*( R(1,num)*P2(1) + R(2,num)*P2(2) + R(3,num)*P2(3) ) + X(num+4,1);
            P3(num) = X(1,1)*( R(1,num)*P3(1) + R(2,num)*P3(2) + R(3,num)*P3(3) ) + X(num+4,1);
            P4(num) = X(1,1)*( R(1,num)*P4(1) + R(2,num)*P4(2) + R(3,num)*P4(3) ) + X(num+4,1);

            P5(num) = X(1,1)*( R(1,num)*P5(1) + R(2,num)*P5(2) + R(3,num)*P5(3) ) + X(num+4,1);
            P6(num) = X(1,1)*( R(1,num)*P6(1) + R(2,num)*P6(2) + R(3,num)*P6(3) ) + X(num+4,1);
            P7(num) = X(1,1)*( R(1,num)*P7(1) + R(2,num)*P7(2) + R(3,num)*P7(3) ) + X(num+4,1);
            end
        end

        v1 = P3-P1;
        v2 = P2-P1;
        n1 = cross(v1,v2);
        mag = norm(n1);
        n1 = n1/mag;

        v1 = P4-P1;
        v2 = P3-P1;
        n2 = cross(v1,v2);
        mag = norm(n2);
        n2 = n2/mag;

        v1 = P6-P1;
        v2 = P5-P1;
        n3 = cross(v1,v2);
        mag = norm(n3);
        n3 = n3/mag;

        v1 = P5-P1;
        v2 = P7-P1;
        n4 = cross(v1,v2);
        mag = norm(n4);
        n4 = n4/mag;

        v1 = P2-P1;
        v2 = P6-P1;
        n5 = cross(v1,v2);
        mag = norm(n5);
        n5 = n5/mag;

        v1 = P7-P1;
        v2 = P4-P1;
        n6 = cross(v1,v2);
        mag = norm(n6);
        n6 = n6/mag;

        if X_diff == 0 && Y_diff == 0
            n=(n1+n2+n3+n4+n5+n6)/6.0;
        elseif X_diff == 0 && Y_diff ~= 0
            n = (n2 + n6)/2.0;
        elseif X_diff ~= 0 && Y_diff == 0
            n = (n1 + n5)/2.0;
        elseif X_diff > 0 && Y_diff > 0 && (X_diff <= Y_diff) %First TIN
            n = n2;
        elseif  X_diff > 0 && Y_diff > 0 && (X_diff > Y_diff) %Second TIN
            n = n1;
        end

        d = dot([-P1(1) -P1(2) -P1(3)],n);

        Coeff(1) = n(1);
        Coeff(2) = n(2);
        Coeff(3) = n(3);
        Coeff(4) = d;

        Z = -(Coeff(1)*Pos(1) + Coeff(2)*Pos(2) + Coeff(4))/Coeff(3);
    else
        Z = 0.0;
    end
else
    Check = 'False';
    Z=0.0;
end