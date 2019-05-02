% if 0
% addpath(genpath(['./rivergithub2/']))
% 
% rang0=[-2235 -2228 546 552]*1e3;
% x0=[rang0(1) rang0(2) rang0(2) rang0(1) rang0(1) ];y0=[rang0(4) rang0(4) rang0(3) rang0(3) rang0(4) ];
% hold on;plot(x0,y0)
% 
% infile='46_18_2_1_2m_v3.0_reg_dem.tif';
% data=readGeotiff(infile,'map_subset',rang0);
% hills=hillshade(double(data.z),data.x,data.y,'plotit');
% %from /home/dai.56/arcticarcticdemappdemapp/river/rivergithub2/maskentropybp1.m 
% end

function [m]=demslopemask(data0r)
% Usage:[m]=demslopemask(data0r);
%input: data0r 

flagplot=0;
windowm=300; %
anglethres=10;

% load RefDEM.mat

data=data0r;clear data0r
z=data;
% res=or.info.map_info.dx;
resx=mean(z.x(2:end)-z.x(1:end-1));resy=mean(z.y(2:end)-z.y(1:end-1));res=mean([abs(resx),abs(resy)]);
P=ones(size(z.z));

z = z.z;
window=round(windowm/res);%200m smoothing window
z = smoothdata(z,'gaussian',window);

resr=2.;dsr=res/resr;
dsr2=1;resr2=res/dsr2; %for better resolution

z(z == -9999) = NaN;
zsv=z;

z = imresize(z,dsr2);
% x = imresize(x,dsr); only needed for hillshade
% y = imresize(y,dsr); only needed for hillshade

% construct hillshade for debugging
%hill = hillshade(z,x,y);

[sx,sy] = gradient(z,resr2); % calculate slopes

[~,rho] = cart2pol(sx,sy); % vector slope

k=11; % convolution kernel size
rhomn =conv2(rho,ones(k)/(k.^2),'same'); % mean slope
rhosd=sqrt( conv2(rho.^2,ones(k)/(k.^2),'same') - rhomn.^2 ); % std dev
rhosd = imresize(rhosd,size(P),'nearest');
z = imresize(zsv,size(P),'nearest');
% Mrho=(rhosd<0.5);   

angle=atan(rho)*180/pi;
M=angle>anglethres;
m.x=data.x;m.y=data.y;m.z=M;

% figure;imagesc(data.x*1e-3,data.y*1e-3,rhosd);colorbar;colormap jet;title('rhosd');view(0,-90)
% hold on;plot(x0*1e-3,y0*1e-3,'r-','linewidth',2)

if flagplot==1
figure;imagesc(data.x*1e-3,data.y*1e-3,angle);colorbar;colormap jet;title('Slope angle');view(0,-90)
hold on;plot(x0*1e-3,y0*1e-3,'r-','linewidth',2)

figure;imagesc(data.x*1e-3,data.y*1e-3,angle<anglethres);colorbar;colormap jet;title('Slope angle<thres');view(0,-90)
hold on;plot(x0*1e-3,y0*1e-3,'r-','linewidth',2)
set(gcf,'Color','white')

[X,Y]=meshgrid(data.x,data.y);
hold on;plot(X(M)*1e-3,Y(M)*1e-3,'m.')
title(['Slope angle > ',num2str(anglethres)])

hold on;plot(X(M),Y(M),'m.')


axis([-2242 -2228 546 559])
end
end

