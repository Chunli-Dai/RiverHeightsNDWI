function [B A L] = General3Dconformal(varargin)

P = varargin{1}; %3D point
G = varargin{2}; %gradient
Angle = varargin{3}; %angle
S = varargin{4}; %scale
observation = varargin{5}; %difference distance


%Cal 3D Conformal Coeff
Coeff = ConformalCoeff3D(P,Angle,S);

A(1) = G(1)*Coeff(1,1) + G(2)*Coeff(2,1) + G(3)*Coeff(3,1); %Scale
A(2) = G(1)*Coeff(1,2) + G(2)*Coeff(2,2) + G(3)*Coeff(3,2); %omega
A(3) = G(1)*Coeff(1,3) + G(2)*Coeff(2,3) + G(3)*Coeff(3,3); %phi
A(4) = G(1)*Coeff(1,4) + G(2)*Coeff(2,4) + G(3)*Coeff(3,4); %kappa
A(5) = G(1)*Coeff(1,5) + G(2)*Coeff(2,5) + G(3)*Coeff(3,5); %tx
A(6) = G(1)*Coeff(1,6) + G(2)*Coeff(2,6) + G(3)*Coeff(3,6); %ty
A(7) = G(1)*Coeff(1,7) + G(2)*Coeff(2,7) + G(3)*Coeff(3,7); %tz

L(1) = observation;

R = RotationMatrix(Angle);

B(1) = G(1)*S*R(1,1) + G(2)*S*R(1,2) + G(3)*S*R(1,3);
B(2) = G(1)*S*R(2,1) + G(2)*S*R(2,2) + G(3)*S*R(2,3);
B(3) = G(1)*S*R(3,1) + G(2)*S*R(3,2) + G(3)*S*R(3,3);
B(4) = -1;
B(5) = -1;
B(6) = -1;

clear P G Angle S ref_H tar_H Coeff R;
