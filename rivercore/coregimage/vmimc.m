function [dx,dy, sigma]=vmimc(infile1,infile2)
% [dx, dy] the displacement infile2 minus infile1.
% sigma: the estmated variance of unit weight
%script_demo
%param_IM2_demo
% clear ; clc; close all
%load and set-up the parameter

% macdir='';
% addpath(genpath([macdir,'/data/chunli/coastline/']));

param_IM2_demo
if ~exist(param.path.module,'dir')
    mkdir param.path.module
end
if ~exist(param.path.subset,'dir')
    mkdir param.path.subset
end
if ~exist(param.path.result,'dir')
    mkdir param.path.result
end
if ~exist(param.path.vmap,'dir')
    mkdir param.path.vmap
end
if ~exist(param.path.tmpdir,'dir')
    mkdir param.path.tmpdir
end

% infile2=['./work/WV03_20170414_104001002B2B4700_104001002C651800_seg6_2m_orthos.tif'];
% infile1=['./work/WV02_20140925_1030010036BD9500_103001003727E300_seg5_2m_orthos.tif'];
% 
% infile2=['./work/WV03_20170414_104001002B2B4700_104001002C651800_seg6_2m_orthos1.tif'];
% infile1=['./work/WV02_20140925_1030010036BD9500_103001003727E300_seg5_2m_orthos1.tif'];
% infile2=['./data/WV02_20150412_10300100400AF500_1030010040C40400_seg1_2m_orthos1668.tif'];
% infile1=['./data/WV02_20140925_1030010036BD9500_103001003727E300_seg6_2m_orthos1668.tif'];

[dir,name,ext] =fileparts(infile1);
param.path.subset=dir;
filename_i0=[name,ext];
[dir,name,ext] =fileparts(infile2);
filename_i1=[name,ext];

%Creat the prior velocity matrix
if 1 %
ic=readGeotiff([param.path.subset,'/',filename_i0]);
ndx=20;res=ic.x(2)-ic.x(1);
v=81:ndx:(length(ic.y)-68);u=81:ndx:(length(ic.x)-85);[V,U]=meshgrid(v,u);
xyuvav(:,3:4)=[U(:) V(:)];
Y=ic.y(V);X=ic.x(U);
xyuvav(:,1:2)=[X(:) Y(:)]+1;
xyuvav(:,5:6)=0.1;
% [Y,X]=meshgrid((ic.y(81)+1):-ndx*res:(ic.y(end-68)+1),(ic.x(81)+1):ndx*res:(ic.x(end-85)+1));
% xyuvav=zeros(length(X(:)),6);xyuvav(:,1:2)=[X(:) Y(:)];
end

param.trk.max_estimated_spd=20;%20;%Maximum estimation of the flow speed in the region, unit: m/yr
param.trk.mpp=2; %GSD of the image, unit: m/pixel. ex) Panchromatic Landsat = 15m/px

tic
%invoke the main program to start processing
MIMC2_split2([param.path.subset,'/',filename_i0],[param.path.subset,'/',filename_i1],xyuvav,param)
% velocity of i1 w.r.t. i0
toc

% save t1.mat xyuvav

% visualize the processes result
filename_vmap=['vmap_',filename_i0(1:14),'_',filename_i1(1:14),'.mat'];
filename_pp_result=['pp_result_',filename_i0(1:14),'_',filename_i1(1:14),'_raw.mat'];

load([param.path.vmap,'/',filename_vmap]);
load([param.path.result,'/',filename_pp_result]);

% dpf3 displacement in pixels 
dist=sqrt(dpf3(:,1).^2+dpf3(:,2).^2);
err=sqrt(err3(:,1).^2+err3(:,2).^2);
drms=std(dist(~isnan(dist))); %std, with mean removed
sigma=rms(err(~isnan(err))); %sigma estimated from residuals (err3)

id=err<3*sigma&abs(dist-nanmean(dist))<3*drms; %control parameters

if 0
figure;hold all
imagesc(ic.x,ic.y,ic.z(:,:,1))
view(0,90)
quiver(xyuvav(id,1),xyuvav(id,2),dpf3(id,1)*res,dpf3(id,2)*res,0,'k');
quiver(xyuvav(id,1),xyuvav(id,2),err3(id,1)*res,err3(id,2)*res,0,'r');
% coordinates in km, displacement in pixel.

[ny,nx]=size(vmap.vx);
vx=reshape(dpf3(:,1),[nx,ny])';
vy=reshape(dpf3(:,2),[nx,ny])';
vxi=reshape(xyuvav(:,1),[nx,ny])';
vyi=reshape(xyuvav(:,2),[nx,ny])';
% idm=reshape(id,[nx,ny])';
% vx(~idm)=NaN;vy(~idm)=NaN;

figure;hold all
imagesc(vxi(1,:)*1e-3,vyi(:,1)*1e-3,vx);colorbar
view(0,90)
title('vx')
colormap jet
% caxis([-2 2]);

figure;hold all
imagesc(vxi(1,:)*1e-3,vyi(:,1)*1e-3,vy);colorbar
view(0,90)
title('vy')
colormap jet
% caxis([-2 2]);
end

dx=mean(dpf3(id,1))*res;
dy=mean(dpf3(id,2))*res;

if isnan(dx) ||isnan(dy)
    dx=0;dy=0;
end

return
end

