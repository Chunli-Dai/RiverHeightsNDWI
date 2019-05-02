%Point, Angle, Scale
function [C11 C12 C13 C14 C15 C16 C17 C21 C22 C23 C24 C25 C26 C27 C31 C32 C33 C34 C35 C36 C37]=ConformalCoeff3D(varargin)

P1 = varargin{1};
P2 = varargin{2};
P3 = varargin{3};
A = varargin{4};
S = varargin{5};

R = RotationMatrix(A);
total_size = size(P1,1);
%X equation
C11 = R(1,1)*P1 + R(2,1)*P2 + R(3,1)*P3; %d_scale 
C12 = zeros(total_size,1); %d_omega 
C13 = ( -sin(A(2))*cos(A(3))*P1 + sin(A(2))*sin(A(3))*P2 + cos(A(2))*P3 )*S; %d_phi
C14 = ( R(2,1)*P1 - R(1,1)*P2 )*S; %d_kappa
C15 = ones(total_size,1);
C16 = zeros(total_size,1);
C17 = zeros(total_size,1);

%Y equation
C21 = R(1,2)*P1 + R(2,2)*P2 + R(3,2)*P3;
C22 = (-R(1,3)*P1 - R(2,3)*P2 - R(3,3)*P3)*S;
C23 =( sin(A(1))*cos(A(2))*cos(A(3))*P1 - sin(A(1))*cos(A(2))*sin(A(3))*P2 + sin(A(1))*sin(A(2))*P3 )*S;
C24 = ( R(2,2)*P1 - R(1,2)*P2 )*S;
C25 = zeros(total_size,1);
C26 = ones(total_size,1);
C27 = zeros(total_size,1);

%Z equation
C31 = R(1,3)*P1 + R(2,3)*P2 + R(3,3)*P3;
C32 = ( R(1,2)*P1 + R(2,2)*P2 + R(3,2)*P3 )*S;
C33 = ( -cos(A(1))*cos(A(2))*cos(A(3))*P1 + cos(A(1))*cos(A(2))*sin(A(3))*P2 - cos(A(1))*sin(A(2))*P3 )*S;
C34 = ( R(2,3)*P1 - R(1,3)*P2 )*S;
C35 = zeros(total_size,1);
C36 = zeros(total_size,1);
C37 = ones(total_size,1);

clear P A S R;



% P = varargin{1};
% A = varargin{2};
% S = varargin{3};
% 
% R = RotationMatrix(A);
% 
% C = zeros(3,7);
% 
% %X equation
% C(1,1) = R(1,1)*P(1) + R(2,1)*P(2) + R(3,1)*P(3); %d_scale 
% C(1,2) = 0.0; %d_omega 
% C(1,3) = ( -sin(A(2))*cos(A(3))*P(1) + sin(A(2))*sin(A(3))*P(2) + cos(A(2))*P(3) )*S; %d_phi
% C(1,4) = ( R(2,1)*P(1) - R(1,1)*P(2) )*S; %d_kappa
% C(1,5) = 1.0;
% C(1,6) = 0.0;
% C(1,7) = 0.0;
% 
% %Y equation
% C(2,1) = R(1,2)*P(1) + R(2,2)*P(2) + R(3,2)*P(3);
% C(2,2) = (-R(1,3)*P(1) - R(2,3)*P(2) - R(3,3)*P(3))*S;
% C(2,3) =( sin(A(1))*cos(A(2))*cos(A(3))*P(1) - sin(A(1))*cos(A(2))*sin(A(3))*P(2) + sin(A(1))*sin(A(2))*P(3) )*S;
% C(2,4) = ( R(2,2)*P(1) - R(1,2)*P(2) )*S;
% C(2,5) = 0.0;
% C(2,6) = 1.0; 
% C(2,7) = 0.0;
% 
% %Z equation
% C(3,1) = R(1,3)*P(1) + R(2,3)*P(2) + R(3,3)*P(3);
% C(3,2) = ( R(1,2)*P(1) + R(2,2)*P(2) + R(3,2)*P(3) )*S;
% C(3,3) = ( -cos(A(1))*cos(A(2))*cos(A(3))*P(1) + cos(A(1))*cos(A(2))*sin(A(3))*P(2) - cos(A(1))*sin(A(2))*P(3) )*S;
% C(3,4) = ( R(2,3)*P(1) - R(1,3)*P(2) )*S;
% C(3,5) = 0.0;
% C(3,6) = 0.0;
% C(3,7) = 1.0;
% 
% clear P A S R;