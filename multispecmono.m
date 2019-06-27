function  [Mstrip]=multispecmono(stripmetafile,wm,ndwisv,rangeov)
%given mono image xml file, orthrectiying, get water mask
% get the shoreline from multispectral images using NDWI index
% Need files: orthorectified imagery, xml files
% Input: 
%        stripmetafile, mono image xml file
%        wm, a priori water mask for getting NDWI statistics; wm.x, wm.y, wm.z, with wm.z = 1 land, 0 water;
% 	rangeov, the overlapping range of the data range and coastal band
% Output: Mstrip, merged M matrix for the whole strip: 1, water, 0 non water, NaN void data.
%         Mstrip.coast, the coastline;
%	v6: save space, Mstrip.z to int8; 1, water, 0 non water, -1 void data.
% ndwisv: saved data of ndwi matrix. The .x, .y coordinates are coregistered to a reference.
% Mstrip.n: to save the NDWI matrix.

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

if isempty(wm.z);step=1; %no apriori water mask
else; step =2; %there is an apriori water mask
end

%flagplot=0;
% to do: try statistics (maximum likelihood ) of training sites

% stripmetafile=[macdir,'/data3/ArcticDEM/region_31_alaska_south/strips/2m/WV03_20170414_104001002B2B4700_104001002C651800_seg6_2m_meta.txt'];
%stripmetafile=[macdir,'/data3/ArcticDEM/region_31_alaska_south/strips/2m//WV02_20160304_1030010052B75A00_1030010053B3BC00_seg3_2m_meta.txt'];
%stripmetafile=[macdir,'/data3/ArcticDEM/region_31_alaska_south/strips/2m/WV02_20120930_103001001A971B00_103001001C13FE00_seg1_2m_meta.txt'];
%/data1/pgc_projects/dai_aleutians_multi_mono/imagery//WV02/WV02_20130923221601_1030010026BD6B00_13SEP23221601-M1BS-500127380110_01_P001.xml 
%get image 1 filenames, and dx2;
    
nsce=1; %mono image
[mdir2,filename,ext]=fileparts(stripmetafile);
for i=1:nsce
mfile{i}=[filename,ext];
end
for i=1:nsce
rmsmetau(i)=0;
dzxydu(i,1:3)=0;
end

if isempty(ndwisv.z) %no saved NDWI matrix.
%ex3, image1, scene 3 5, orthorectify without dem
% infile=[macdir,'/home/chunli/coastline/orthorect1/WV03_20170414222453_104001002BA9FA00_17APR14222453-M1BS-501325596050_01_P006t1.tif'];
% refers to Maglione 2014
%Coastal 0.400 - 0.450; Blue 0.450 - 0.510; Green 0.510 - 0.580; Yellow 0.585 - 0.625
%Red 0.630 - 0.690; Red-Edge 0.705 - 0.745; NIR1 0.770 - 0.895; NIR2 0.860 - 0.940
%GEOEYE-1: 1, 2, 3, 4 BLUE GREEN RED NIR1
[~,filename,~]=fileparts(stripmetafile);
satname=filename(1:4);
if strcmp(satname,'GE01')|| strcmp(satname,'QB02') ||strcmp(satname,'IK01')
  iBlue=1;iGreen=2;iRed=3;iNIR1=4;
  ig=[iBlue,iGreen,iRed,iNIR1];
  bandstr={'BAND_B','BAND_G','BAND_R','BAND_N'};
  abscalfactor=zeros(4,1);effectivebandwith=zeros(4,1);
  flagb=1;
elseif strcmp(satname,'WV02') ||strcmp(satname,'WV03') 
    iNIR2=8;iCoastal=1;iRed=5;iGreen=3; iNIR1=7;
    iBlue=2;iYellow=4;iRE=6;
    ig=[iCoastal,iGreen,iRed,iNIR2,iNIR1,iBlue,iYellow,iRE];
    bandstr={'BAND_C','BAND_B','BAND_G','BAND_Y','BAND_R','BAND_RE','BAND_N','BAND_N2'};
    abscalfactor=zeros(8,1);effectivebandwith=zeros(8,1);
    flagb=2;
end
[GAIN,OFFSET,Esun]=readgainoffset(satname); %read Gain OFFset data, GainOffset.txt.

Mstrip=struct(); Mstrip.x=[];Mstrip.y=[];Mstrip.z=[];Mstrip.coast=[];Mstrip.n=[];
% %orthorectify usnig gdal
for is=1:nsce
ntffile=strrep(mfile{is},'.xml','.ntf');
tiffile=strrep(mfile{is},'.xml','.tif');
%[status , cmdout ]=system(['ls ',multidir,'/*/',ntffile]);%status =0 if found; but too strict on multidir.
[status , cmdout ]=system(['find ',multidir,' -name ',ntffile]); %returns empty for cmdout when not found. status=0 always
[status2 , cmdout2 ]=system(['find ',multidir,' -name ',tiffile]);
if  ~isempty(cmdout) %if .ntf file is found, tif file is produced and stored in orthworkdir.
    ntffile=deblank(cmdout);
end
if 0
    [dir,name,ext] =fileparts(ntffile);
%   tiffile=deblank([dir,'/',mfile{is}]);
    tiffile=deblank([orthworkdir,'/',mfile{is}]); %to generate the tif file
    tiffile=strrep(tiffile,'.xml','.tif');
end
if ~isempty(cmdout2) %if .tif file is found.
    tiffile=deblank(cmdout2);	
else; 
	warning(['Multispectral file ',tiffile,' not found!'])
	continue
end
if 0&& ~exist(tiffile,'file')
display(['time gdalwarp -rpc -et 0.01 -co tiled=yes -co compress=lzw -t_srs EPSG:3413 ',ntffile,' ',tiffile])
system(['time gdalwarp -rpc -et 0.01 -co tiled=yes -co compress=lzw -t_srs EPSG:3413 ',ntffile,' ',tiffile]);
end
if exist(tiffile,'file')
    infile=tiffile;
    try
    if ~isempty(rangeov)
    data=readGeotiff(infile,'map_subset', rangeov);
    else
    data=readGeotiff(infile); %% to keep long profiles
    end
    catch e
    %tif provided by pgc might not be compatable with readGeotiff.
      fprintf('There was an error! The message was:\n%s',e.message);
      %fprintf('Read tif file using importdata instead of readGeotiff.')
      continue
    end
    [~,~,nb]=size(data.z);
    if nb == max(ig)
	%good, do nothing
    elseif nb < max(ig) && (nb==4&&length(ig)==8)
	%use the classic 4 bands instead of the 8 bands for WV02 WV03
	flagb=1;
	iBlue=1;iGreen=2;iRed=3;iNIR1=4;
        ig=[iBlue,iGreen,iRed,iNIR1];
        bandstr={'BAND_B','BAND_G','BAND_R','BAND_N'};
        abscalfactor=zeros(4,1);effectivebandwith=zeros(4,1);

	GAIN(1)=GAIN(2);GAIN(2)=GAIN(3);GAIN(3)=GAIN(5);GAIN(4)=GAIN(7);
	OFFSET(1)=OFFSET(2);OFFSET(2)=OFFSET(3);OFFSET(3)=OFFSET(5);OFFSET(4)=OFFSET(7);
	Esun(1)=Esun(2);Esun(2)=Esun(3);Esun(3)=Esun(5);Esun(4)=Esun(7);
    else
       warning(['Bands are different as anticipated for ',tiffile])
       continue 
    end
else
    fprintf([tiffile,' does not exist. \n'])
    continue
end

% read .xml meta file
% metafile=[macdir,'/data3/ArcticDEM/region_31_alaska_south/strips/2m/WV03_20170414_104001002CD4E700_104001002BA9FA00_seg2_2m_meta.txt'];
%/data3/ArcticDEM/region_31_alaska_south/strips/2m/WV03_20170414_104001002CD4E700_104001002BA9FA00_seg2_2m_meta.txt 
infile=mfile{is};
[status , cmdout ]=system(['find ',multidir,' -name ',infile]);
if  ~isempty(cmdout) %if  file is found
    metafile=deblank(cmdout);
else
    warning(['Metafile ',metafile,' not found!'])
    continue
end
%metafile=deblank(strrep(ntffile,'.ntf','.xml'));
%get ABSCALFACTOR EFFECTIVEBANDWIDTH for each band, and MEANSUNEL
c=textread(metafile,'%s','delimiter','\n');
% fid=fopen(metafile);
% n = linecount(fid);
% c=cell(n,1);
% for i=1:n
%   c{i} = fgetl(fid);
% end
r=find(~cellfun(@isempty,strfind(c,'MEANSUNEL')));
if isempty(r)
warning(['XML format is different as anticipated.',stripmetafile])
Mstrip=struct(); Mstrip.x=[];Mstrip.y=[];Mstrip.z=[];Mstrip.coast=[];Mstrip.n=[];
return
end
c2=c{r};
r1=strfind(c2,'>');r2=strfind(c2,'</');c2([1:r1(1),r2(1):end])='';
meanSunEl=sscanf(c2, '%g', 1);
for i=ig % bands 1 to 8
r=find(~cellfun(@isempty,strfind(c,bandstr(i))));
%/data1/pgc_projects/coastline/imagery/WV02_20120930_103001001A971B00_103001001C13FE00/WV02_20120930220825_103001001A971B00_12SEP30220825-M1BS-052903623070_01_P001.xml
if isempty(r)
warning(['Bands are different as anticipated.'])
return;
end %
c1=c(r(1):(r(2)));
str='ABSCALFACTOR';
r=find(~cellfun(@isempty,strfind(c1,str)));
c2=c1{r};r1=strfind(c2,'>');r2=strfind(c2,'</');c2([1:r1(1),r2(1):end])='';
% Xbs=deblank(strrep(c1{r(1)},str,''));
abscalfactor(i)= sscanf(c2, '%g', 1);
str='EFFECTIVEBANDWIDTH';
r=find(~cellfun(@isempty,strfind(c1,str)));
c2=c1{r};r1=strfind(c2,'>');r2=strfind(c2,'</');c2([1:r1(1),r2(1):end])='';
effectivebandwith(i)=sscanf(c2, '%g', 1);
end

Theta=(90.-meanSunEl)*pi/180.;
dES=1.; %AU
% sun-earth distance polynomial function coefficients
% doy = day of the year: 
year=sscanf(filename(6:9), '%g', 1); month=sscanf(filename(10:11), '%g', 1); day=sscanf(filename(12:13), '%g', 1);
doy=juliandate(year,month,day)-juliandate(year,1,1)+1;
C = [1.8739e-26,-3.4455e-23,2.7359e-20,-1.2296e-17,3.0855e-15,-2.2412e-13,-5.8744e-11,6.9972e-10,2.5475e-06,-1.6415e-05,0.9833];
dES = polyval(C,doy); % earth-sun distance

rhog=zeros(size(data.z));
for i=ig %band id 
    DN=double(data.z(:,:,i));
    L=GAIN(i)*DN*(abscalfactor(i)/effectivebandwith(i))+OFFSET(i);
    rho=L*dES^2*pi/(Esun(i)*cos(Theta));
    rhog(:,:,i)=rho;
end

if strcmp(satname,'GE01')|| strcmp(satname,'QB02') ||strcmp(satname,'IK01')||flagb==1
% Green=double(data.z(:,:,2));NIR1=double(data.z(:,:,4)); 
Green=rhog(:,:,iGreen);NIR1=rhog(:,:,iNIR1);
NDWI=double(Green-NIR1)./double(Green+NIR1);%McFeeters 1996
elseif strcmp(satname,'WV02') ||strcmp(satname,'WV03') 
% Green=double(data.z(:,:,3));NIR2=double(data.z(:,:,8)); 
% Coastal=double(data.z(:,:,1)); Red=double(data.z(:,:,5));
Green=rhog(:,:,iGreen);NIR2=rhog(:,:,iNIR2);Coastal=rhog(:,:,iCoastal);
% NDWI=double(Green-NIR2)./double(Green+NIR2);%McFeeters 1996
NDWI=double(Coastal-NIR2)./double(Coastal+NIR2);
%[mean(NDWI(mp)), mean(NDWI(mpd)) , mean(NDWI(mp))-mean(NDWI(mpd))]= 1.1 for radiance, 0.9 for digital number;
end
%rhog can be 13GB
clear rhog 

NDWI(data.z(:,:,iGreen) == 0)=nan;%set the image edge of M to -1.

dx2=dzxydu(is,2:3);

% collect overlapping M, choose value from (NaN, value), and choose 1 from
% (value, 1); <-> choose the larger value of M
if is==1 || isempty(Mstrip.z)
%     Medgs=(data.z(:,:,1) == 0);
    Mstrip.x=data.x-dx2(1);Mstrip.y=data.y-dx2(2);Mstrip.n=NDWI; %Mstrip.z=M;
    resx=mean(data.x(2:end)-data.x(1:end-1));resy=mean(data.y(2:end)-data.y(1:end-1));
    % image 1, res=1.8312; image2, res=1.8097
else %ignore this
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

if 0% flagplot==1
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

else % used the saved NDWI matrix.

NDWI=ndwisv.z; data.x=ndwisv.x;data.y=ndwisv.y;
Mstrip.x=data.x;Mstrip.y=data.y;Mstrip.n=NDWI; %applied the coregistration
end

if isempty(Mstrip.n)
Mstrip=struct(); Mstrip.x=[];Mstrip.y=[];Mstrip.z=[];Mstrip.coast=[];Mstrip.n=[];
return
end

if step ==2   %step 2, use the refined threshold.
%retrieving statistics of NDWI over a priori water mask.
% wm: 1 water, 0 land, 
%resrc=abs(wm.x(2)-wm.x(1));
resx=mean(data.x(2:end)-data.x(1:end-1));resy=mean(data.y(2:end)-data.y(1:end-1));
resrc=mean([abs(resx),abs(resy)]);
% Md1 = imdilate(wm.z, ones(widthstat/resrc)); %make sure water is absolutely water.
Md1=wm.z==1;

wm1 = interp2(wm.x,wm.y,Md1,data.x,data.y','*nearest',0); %out of region: 0, non water.
wm1sv=wm1;
ndwil=NDWI(wm1==1); %water
mean1=nanmean(ndwil);std1=nanstd(ndwil);cnt1=sum(~isnan(ndwil));
b=ndwil(~isnan(ndwil));

if 0 %old method
Md1 = imdilate(wm.z~=0, ones(600/resrc)); %make sure land is absolutely land. 300 m width expansion
% Me1=imerode(wm.z, ones(widthstat/resrc)); %land
wm1 = interp2(wm.x,wm.y,Md1,data.x,data.y','*nearest',1); %out of region: 1,  water.
else
  %widthriv=600; %to test
  widthriv=100; %to test
  if isfield(wm, 'w50') %new way: use rivers for water; land exclude all water bodies.
    Md1=wm.w50~=0;
  else %old way; treat rivers and lakes the same.
    Md1=wm.z~=0; %first interpolate; then dilate to save computation time.
  end
  wm1 = interp2(wm.x,wm.y,Md1,data.x,data.y','*nearest',1); %out of region: 1,  water.
  %Only use a 300m width land on both sides of the river
  Med=imdilate(isnan(NDWI),ones(round(widthriv*3/resrc))); %image edge
  Md1 = imdilate(wm1~=0, ones(round(widthriv/resrc))); %make sure land is absolutely land. 300 m width expansion
  Mdd = imdilate(Md1, ones(round(widthriv/resrc)));
  wm1=~((Mdd-Md1)& ~Med);
end

ndwil=NDWI(wm1==0); %land
mean2=nanmean(ndwil);std2=nanstd(ndwil);cnt2=sum(~isnan(ndwil));
a=ndwil(~isnan(ndwil));

%Oct 8, 2018; Simple Adaptive threshold 
%Modification Sept. 25, 2018: Delete scenes with bad STD.
%istdthres=min([stdthres,(mean1 - mean2)/2]);
istdthres=stdthres;
badflag=0;
if (isnan(mean1)||cnt1<cntmin) && (isnan(mean2)||cnt2<cntmin)
    badflag=1;ithreshold =threshold;
elseif (isnan(mean1)||cnt1<cntmin) && (mean2<0.2&&std2<istdthres) || (isnan(mean2)||cnt2<cntmin)&&(mean1>0.8&&std1<istdthres)
        ithreshold =threshold;
elseif (std1<istdthres && std2<istdthres && (mean1-mean2)>dmthres)
        ithreshold=(mean1+mean2)/2;
        if isnan(ithreshold);warning('Estimated threshold is NAN');ithreshold=threshold;end
else %discard the image.
    badflag=1; ithreshold =threshold;
end
fprintf('\n %f %f %f %f %f %d (NDWI mean std thres) %s \n',mean1,std1,mean2,std2,ithreshold,badflag,stripmetafile);
if 0
fid = fopen('./ndwistats.dat','a'); %statistics
fprintf(fid,' %f %f %f %f %f %d (NDWI mean std thres) %s \n',mean1,std1,mean2,std2,ithreshold,badflag,stripmetafile);
fclose(fid);
end

elseif step==1 %step 1, use the general threshold.
    badflag=0;
    ithreshold =threshold;
    mean1=0; std1=0;mean2=0;std2=0;
end


if flagplot==1
figure;
hold all
set(gcf,'Color','white');set(gca,'FontSize', 18);set(gcf, 'PaperPosition', [0 0 6 4]);set(gcf, 'PaperSize', [ 6 4]);
imagesc(Mstrip.x*1e-3,Mstrip.y*1e-3,double(Mstrip.n))
colormap jet;colorbar
%hold on;plot(X(M)*1e-3,Y(M)*1e-3,'ro')
title(['NDWI ',filename(6:13)]);
axis equal
xlabel('x (km)');ylabel('y (km)')
caxis([-1 1])
print('-dpng','-r400',filename)
saveas(gcf,[filename],'fig')
end

if badflag==1    
Mstrip=struct(); Mstrip.x=[];Mstrip.y=[];Mstrip.z=[];Mstrip.coast=[];Mstrip.n=[];
close all
return
end

M=int8(NDWI>ithreshold);%0.3; %threshold
M(isnan(NDWI))=-1;%set the image edge of M to -1.

Mstrip.z=M;

%apply a priori river buffer zone.
if isfield(wm, 'buf')
buf= interp2(wm.x,wm.y,wm.buf,Mstrip.x,Mstrip.y','*nearest',0); %out of region: 0, void
Mstrip.z(~buf)=-1; 
%Do not remove clusters for fair comparison of river stages.
%Since Removing small clusters can remove different lakes from different masks, depending on threshold.
%Decide To remove the clusters, 1) for better calculation of river width; 
%2) only small lakes (associate with lower stage) will be removed, so the lowest stage will still be corrected selected.
%Update: do not apply the buffer until last step; since it will affect the water probability map.
end

    
%coast output not used.
%[M,~]=mask2coast(Mstrip);
[M,Modfil]=mask2coast(Mstrip);
if isfield(wm, 'buf') %same control; do not apply the removal of lakes until last step.
Mstrip.z=int8(Modfil); %remove clusters and lakes, for better calculation of river width out of water mask.
end

[X,Y]=meshgrid(Mstrip.x,Mstrip.y);

if flagplot==1
hold on;plot(X(M)*1e-3,Y(M)*1e-3,'k.')
title(['NDWI ',filename(6:13)]);
axis equal
print('-dpng','-r400',filename)
saveas(gcf,[filename],'fig')
end

Mstrip.coast=M;%[X(M),Y(M)]; %merged coast for the strip

close all

end


