function [es]=getes(data)
% Do imshow(data.z) to find the sides (North East South West) of image that intersects river.
% Input: data: data.z, a water mask matrix, 1 is water, 0 is nonwater.
%	       data.x, data.y polar stereographic coordinates.
% Output: es, two-element character string specifying the sides of the image
%               that the channel intersects (e.g. 'SW', 'NS', etc.).
%		e.g. 'WE': % START stream at West, End stream at East.
%         

resx=mean(data.x(2:end)-data.x(1:end-1));resy=mean(data.y(2:end)-data.y(1:end-1));
res=mean([abs(resx),abs(resy)]);
resr=200;%m
dsr=res/resr;
datar.x= imresize(data.x,dsr); % 0.01
datar.y= imresize(data.y,dsr);
datar.z= imresize(data.z,dsr);
[X,Y]=meshgrid(datar.x,datar.y);
x=X(datar.z==1);y=Y(datar.z==1);

%Fitting x y using y=kx+b
n=length(x);
A=[x(:), ones(n,1)];
yobs=y(:);
ksi=A\yobs;
yfit=A*ksi;
k=ksi(1);

%figure;hold on;plot(x*1e-3,y*1e-3,'r.',x*1e-3,yfit*1e-3,'r-','linewidth',3)

if (abs(k)<=1) %
   es='EW';
else
   es='NS';
end


return
end
