function [R]=RotationMatrix(varargin)

Angle = varargin{1};

omega = Angle(1);
phi = Angle(2);
kappa = Angle(3);

R_omega = eye(3);
R_omega(2,:) = [0.0  cos(omega) sin(omega)];
R_omega(3,:) = [0.0 -sin(omega) cos(omega)];

R_phi = eye(3);
R_phi(1,:) = [cos(phi) 0 -sin(phi)];
R_phi(3,:) = [sin(phi) 0  cos(phi)];

R_kappa = eye(3);
R_kappa(1,:) = [ cos(kappa) sin(kappa) 0];
R_kappa(2,:) = [-sin(kappa) cos(kappa) 0];

R = R_kappa*R_phi*R_omega;


%R(1,1) = cos(phi)*cos(kappa);
%R(1,2) = sin(omega)*sin(phi)*cos(kappa) + cos(omega)*sin(kappa);
%R(1,3) = -cos(omega)*sin(phi)*cos(kappa) + sin(omega)*sin(kappa);

%R(2,1) = -cos(phi)*sin(kappa);
%R(2,2) = -sin(omega)*sin(phi)*sin(kappa) + cos(omega)*cos(kappa);
%R(2,3) = cos(omega)*sin(phi)*sin(kappa) + sin(omega)*cos(kappa);

%R(3,1) = sin(phi);
%R(3,2) = -sin(omega)*cos(phi);
%R(3,3) = cos(omega)*cos(phi);


clear Angle omega phi kappa;
