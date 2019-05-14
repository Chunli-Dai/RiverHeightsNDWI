function  [Mstrip]=multispecstrip(stripmetafile,pzxy,wm,fis,dX4Sg)
%given strip meta file, finding all image 1 multispectral imageries, orthrectiying, get water mask
% get the shoreline from multispectral images using NDWI index
% Need files: orthorectified imagery, xml files
% coregistration: strip meta.txt file 
% Input: 
%        stripmetafile, strip meta file
%        wm, a priori water mask for getting NDWI statistics; wm.x, wm.y, wm.z, with wm.z = 1 land, 0 water;
% Output: Mstrip, merged M matrix for the whole strip: 1, water, 0 non water, NaN void data.
%         Mstrip.coast, the coastline;
%	v6: save space, Mstrip.z to int8; 1, water, 0 non water, -1 void data.

% Msv: a list of saved water masks of mono images; It could be from multispectral or panchromatic images; 
% fis: name list of Msv;
% dX4Sg: the offset applied to water mask, applied it back to mosaic the strip water mask.

if 0
addpath(genpath([macdir,'/data/chunli/scripts/']));
addpath(genpath([macdir,'/data/chunli/coastline']));
%macdir='/Users/chunlidai/surge/';
%macdir='';%

%control parameters
width=2e3; %buffer width of the a priori coastline, e.g., 2km.
threshold=0.3; %0.3; %threshold of NDWI for water classification.
orthworkdir=[macdir,'/data/chunli/coastline/orthorectwork/'];
end % 0

%load and set up the parameter
constant

flagplot=0;

[~,filename,~]=fileparts(stripmetafile);
satname=filename(1:4);

% to do: try statistics (maximum likelihood ) of training sites

% stripmetafile=[macdir,'/data3/ArcticDEM/region_31_alaska_south/strips/2m/WV03_20170414_104001002B2B4700_104001002C651800_seg6_2m_meta.txt'];
%stripmetafile=[macdir,'/data3/ArcticDEM/region_31_alaska_south/strips/2m//WV02_20160304_1030010052B75A00_1030010053B3BC00_seg3_2m_meta.txt'];
%stripmetafile=[macdir,'/data3/ArcticDEM/region_31_alaska_south/strips/2m/WV02_20120930_103001001A971B00_103001001C13FE00_seg1_2m_meta.txt'];
%get image 1 filenames, and dx2;
c=textread(stripmetafile,'%s','delimiter','\n');
rs=find(~cellfun(@isempty,strfind(c,'scene')));
nsce=length(rs)-1; %number of scenes
rmsmeta=zeros(nsce,1);idd=[];dzxyd=zeros(nsce,3);mfile=cell(nsce,1);
for i=1:nsce
c1=c{rs(1)+i};r1=strfind(c1,'dem.tif');c1(1:(r1+6))='';
tmp=sscanf(c1, '%g', 4);
rmsmeta(i)=tmp(1);dzxyd(i,1:3)=tmp(2:4);

c1=c{rs(1+i)+3};r1=strfind(c1,'Image 1=');if(isempty(r1)) Warning('Image 1 not found') ; end
r1=strfind(c1,'/');c1(1:(r1(end)))='';

c11=deblank(c1);
satname=c11(1:4);%use this, since strip could be W1W2_20110825_1020010014850E00_103001000C5BD600_seg1_2m_meta.txt

if strcmp(satname,'WV01')
    mfile{i}=deblank(c1);
else
    mfile{i}=deblank(strrep(c1,'P1BS','M1BS'));%e.g.WV02_20160304214247_1030010052B75A00_16MAR04214247-M1BS-500641617080_01_P009.tif
end

end
%find the unique images, delete the repeat file names, average the dzxy, and rms
[un idx_last idx] = unique(mfile,'stable');
unique_idx = accumarray(idx(:),(1:length(idx))',[],@(x) {(x)});
mfile=un;nsce=length(idx_last);
rmsmetau=zeros(nsce,1); dzxydu=zeros(nsce,3);
for i=1:nsce
rmsmetau(i)=mean(rmsmeta(unique_idx{i}));
dzxydu(i,1:3)=mean(dzxyd(unique_idx{i},:),1);
end

Mstrip=struct(); Mstrip.x=[];Mstrip.y=[];Mstrip.z=[];Mstrip.coast=[];
% %orthorectify usnig gdal
for is=1:nsce
ntffile=strrep(mfile{is},'.tif','.ntf');
ntffile(end-3:end)='';%delete the extension.
rs=find(~cellfun(@isempty,strfind(fis,ntffile)));

if ~isempty(rs) %if water mask file is found.
    rs=rs(1);
%   data=Msv(rs);
    
    infile=fis{rs};
    rangeov=[];

    %check whether the data is multispectral image
    name=ntffile;
    flagfmt2=0; %panchromatic band
    r1=strfind(name,'M1BS');
    if ~isempty(r1); flagfmt2=1;end

    ndwisv=struct('x',[],'y',[],'z',[]);%Initialize
    if flagfmt2==1 %multispectral image
        clear data
        data=multispecmono(infile,wm,ndwisv,rangeov); %given strip meta file, finding all image 1 multispectral imageries, orthrectiying, get water mask.

    else %stereo orthorectified panchromatic image; or WV01 mono panchromatic image
%         infile=  strrep([demdir,'/',f{i}],'meta.txt','dem.tif');
        data=maskentropy(infile,wm,ndwisv,rangeov);
    end
 
 
    data.z=int8(data.z);%avoid bug 1, Error using max,
    M=data.z;
    tzxy3=dX4Sg(rs,1:3);

    %fixing error in line 'M2.x=Mstrip.x(Mstrip.x>=min(xt)&Mstrip.x<=max(xt));':Error using >=Matrix dimensions must agree.
    if isempty(data.z)
	warning(['Scene Water mask file ',ntffile,' empty!'])
        continue
    end
else 
	warning(['Water mask file ',ntffile,' not found!'])
	continue
end

dx2=dzxydu(is,2:3);%old method, not accurate, since the mono orthoimage and strip orthoimage use different DEMs.
dx2=tzxy3(2:3)-pzxy(2:3);%use the coregistration results(mono to strip orthoimage), accurate within 0.2m 
%if in pxpymono.m, the mono image of strip file use its own strip orthoimage as reference.

% collect overlapping M, choose value from (NaN, value), and choose 1 from
% (value, 1); <-> choose the larger value of M
if is==1 || isempty(Mstrip.z)
%     Medgs=(data.z(:,:,1) == 0);
    Mstrip.x=data.x-dx2(1);Mstrip.y=data.y-dx2(2);Mstrip.z=data.z;
    resx=mean(data.x(2:end)-data.x(1:end-1));resy=mean(data.y(2:end)-data.y(1:end-1));
    % image 1, res=1.8312; image2, res=1.8097
else
    Mpre=Mstrip;
    Mstrip.x=min([Mpre.x(:);data.x(:)]):resx:max([Mpre.x(:);data.x(:)]);%coverage for all images
    Mstrip.y=max([Mpre.y(:);data.y(:)]):resy:min([Mpre.y(:);data.y(:)]);
    if 0 %slow, 20sec each interpolation for two images, 54sec
        M1z = interp2(Mpre.x,Mpre.y,Mpre.z,Mstrip.x,Mstrip.y','*nearest',NaN);
        M2z = interp2(data.x-dx2(1),data.y-dx2(2),M,Mstrip.x,Mstrip.y','*nearest',NaN);
        Mstrip.z=max(M1z,M2z);%6sec
    else %faster, only interpolate at each image coverage and compare at the overlapping area, 30 sec
	% maybe union.m functin also works.
%       tic
        Mstrip.z=-1*ones(length(Mstrip.y),length(Mstrip.x),'int8'); %initialize the mosaiced strip mask
        xt=data.x-dx2(1);yt=data.y-dx2(2);
        M2.x=Mstrip.x(Mstrip.x>=min(xt)&Mstrip.x<=max(xt));M2.y=Mstrip.y(Mstrip.y<=max(yt)&Mstrip.y>=min(yt));
        M2.z=interp2(xt,yt,M,M2.x,M2.y','*nearest',-1);
        Mstrip.z(Mstrip.y<=max(yt)&Mstrip.y>=min(yt),Mstrip.x>=min(xt)&Mstrip.x<=max(xt))=M2.z;

        %do the same for Mpre
        xt=Mpre.x;yt=Mpre.y;
        M1.x=Mstrip.x(Mstrip.x>=min(xt)&Mstrip.x<=max(xt));
        M1.y=Mstrip.y(Mstrip.y<=max(yt)&Mstrip.y>=min(yt));
        flag1=0;
        if length(M1.x)==length(Mpre.x) && length(M1.y)==length(Mpre.y)
            df=max(abs([M1.x(:)-Mpre.x(:);M1.y(:)-Mpre.y(:)]));
            if df < 1e-9
                M1.z=Mpre.z;
                flag1=1;
            end
        end
        if flag1==0
            M1.z=interp2(Mpre.x,Mpre.y,Mpre.z,M1.x,M1.y','*nearest',-1);
        end
        Mstrip.z(Mstrip.y<=max(yt)&Mstrip.y>=min(yt),Mstrip.x>=min(xt)&Mstrip.x<=max(xt))=M1.z;
        % pick the larger values at the overlapping area.
        rangref=[min(M1.x) max(M1.x) min(M1.y) max(M1.y)];
        rangtar=[min(M2.x) max(M2.x) min(M2.y) max(M2.y)];
        rangeov=[max(rangtar(1),rangref(1)),min(rangtar(2),rangref(2)), max(rangtar(3),rangref(3)),min(rangtar(4),rangref(4))];
        M1z=M1.z(M1.y<=rangeov(4)&M1.y>=rangeov(3),M1.x>=rangeov(1)&M1.x<=rangeov(2));
        M2z=M2.z(M2.y<=rangeov(4)&M2.y>=rangeov(3),M2.x>=rangeov(1)&M2.x<=rangeov(2));
        Mstrip.z(Mstrip.y<=rangeov(4)&Mstrip.y>=rangeov(3),Mstrip.x>=rangeov(1)&Mstrip.x<=rangeov(2))=max(M1z,M2z);
%       toc
    end
end

if flagplot==1
figure;
hold all
set(gcf,'Color','white');set(gca,'FontSize', 18);set(gcf, 'PaperPosition', [0 0 6 4]);set(gcf, 'PaperSize', [ 6 4]);
imagesc(Mstrip.x,Mstrip.y,double(Mstrip.z));
title('Water Mask')
colorbar
caxis([-1 1])
% title(['NDWI Image ',num2str(is)])
% caxis([-0.3 0.1])

%plot true color of the image
clear ic
maxz=255;%max(max([rhog(:,:,iRed),rhog(:,:,iGreen),rhog(:,:,iBlue)]));
ic.z=rhog(:,:,iRed)/(maxz)*255;
ic.z(:,:,2)=rhog(:,:,iGreen)/(maxz)*255;
ic.z(:,:,3)=rhog(:,:,iBlue)/(maxz)*255;
% ic.z=uint8(ic.z);
figure,imshow(ic.z);

clear ic
maxz=max(max([data.z(:,:,iRed),data.z(:,:,iGreen),data.z(:,:,iBlue)]));
ic.z=double(data.z(:,:,iRed))/double(maxz)*255;
ic.z(:,:,2)=double(data.z(:,:,iGreen))/double(maxz)*255;
ic.z(:,:,3)=double(data.z(:,:,iBlue))/double(maxz)*255;
ic.z=uint8(ic.z);
figure,imshow(ic.z);

end

if 0
[LAT,LON]=polarstereo_inv(X(M)-dx2(1),Y(M)-dx2(2),[],[],70,-45);
output=[LAT(:),LON(:)];
% ofile=[''];
save coastndwi.dat output -ascii 
end

end % for is=1:nsce

if isempty(Mstrip.z)
return
end

%[M,~]=mask2coast(Mstrip);
[M,Modfil]=mask2coast(Mstrip);
%Use the mosaiced mask, do not alter
%if isfield(wm, 'buf') %same control; do not apply the removal of lakes until last step.
%Mstrip.z=int8(Modfil); %remove clusters and lakes, for better calculation of river width out of water mask.
%end

[X,Y]=meshgrid(Mstrip.x,Mstrip.y);

if flagplot==1
figure;imagesc(Mstrip.x,Mstrip.y,double(Modfil))
hold on;plot(X(M),Y(M),'ro')
end

Mstrip.coast=M;%[X(M),Y(M)]; %merged coast for the strip
close all

end
