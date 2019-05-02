%point, grident, angle, scale
function [A AW L]= CoeffMatrix_25D(varargin)

X = varargin{1}; %3D point
Y = varargin{2};
Z = varargin{3};
Gx = varargin{4}; %gradient
Gy = varargin{5}; %gradient
Gz = varargin{6}; %gradient
Angle = varargin{7}; %angle
S = varargin{8}; %scale
observation = varargin{9}; %difference distance
weight = varargin{10};

%Cal 3D Conformal Coeff
[C11 C12 C13 C14 C15 C16 C17 C21 C22 C23 C24 C25 C26 C27 C31 C32 C33 C34 C35 C36 C37] = ConformalCoeff3D(X,Y,Z,Angle,S);

A1 = (Gx.*C11 + Gy.*C21 + Gz.*C31); %Scale
A2 = (Gx.*C12 + Gy.*C22 + Gz.*C32); %omega
A3 = (Gx.*C13 + Gy.*C23 + Gz.*C33); %phi
A4 = (Gx.*C14 + Gy.*C24 + Gz.*C34); %kappa
A5 = (Gx.*C15 + Gy.*C25 + Gz.*C35); %tx
A6 = (Gx.*C16 + Gy.*C26 + Gz.*C36); %ty
A7 = (Gx.*C17 + Gy.*C27 + Gz.*C37); %tz

A = [A1 A2 A3 A4 A5 A6 A7];

AW = [A1.*weight A2.*weight A3.*weight A4.*weight A5.*weight A6.*weight A7.*weight];

L = observation;

clear P G Angle S ref_H tar_H Coeff;



% function [A AW L]= CoeffMatrix_25D(varargin)
% 
% P = varargin{1}; %3D point
% G = varargin{2}; %gradient
% Angle = varargin{3}; %angle
% S = varargin{4}; %scale
% observation = varargin{5}; %difference distance
% weight = varargin{6};
% 
% %Cal 3D Conformal Coeff
% Coeff = ConformalCoeff3D(P,Angle,S);
% 
% A(1) = (G(1)*Coeff(1,1) + G(2)*Coeff(2,1) + G(3)*Coeff(3,1)); %Scale
% A(2) = (G(1)*Coeff(1,2) + G(2)*Coeff(2,2) + G(3)*Coeff(3,2)); %omega
% A(3) = (G(1)*Coeff(1,3) + G(2)*Coeff(2,3) + G(3)*Coeff(3,3)); %phi
% A(4) = (G(1)*Coeff(1,4) + G(2)*Coeff(2,4) + G(3)*Coeff(3,4)); %kappa
% A(5) = (G(1)*Coeff(1,5) + G(2)*Coeff(2,5) + G(3)*Coeff(3,5)); %tx
% A(6) = (G(1)*Coeff(1,6) + G(2)*Coeff(2,6) + G(3)*Coeff(3,6)); %ty
% A(7) = (G(1)*Coeff(1,7) + G(2)*Coeff(2,7) + G(3)*Coeff(3,7)); %tz
% 
% AW = A*weight;
% 
% L(1) = observation;
% 
% clear P G Angle S ref_H tar_H Coeff;
