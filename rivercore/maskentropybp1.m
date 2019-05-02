function [m] = maskentropy(ifile,wm,ndwisv,rangeov)
% Shoreline Detection.
%
% [m] = maskentropy(demFile) returns the shoreline mask structure m.x, m.y, m.z for the
%dem in demFile. Mask uses the DEM, orthoimage and matchtag files, which
%are assumed to be in the same path with the same basename and standard
%suffixes. Mask combines match point-density and orthoimage entropy to
%locate large water bodies and uses the match point-density and the
%standard deviation of slope to locate both water and cloads

% Preparation: put orthoimage, metafile, matchtag, dem files in the same directory with the same basename 
%		and standard suffixes.
% Input: 
% demFile,   filename.
% Output: 
% m, mask structure (m.x,m.y,m.z), 
%	m.z matrix, the pixel is shoreline if it is 1 and it is not shoreline if it is 0.
%  m.z, pixel is water if it is 1, land if 0.
%  m.n, brightness.
% Ian Howat  06-Apr-2017 
% Chunli Dai July 2017
% Chunli Dai December  2017
% Nov 2018: to work with both scene and mono images.

constant

if isempty(wm.z);step=1; %no apriori water mask
else; step =2; %there is an apriori water mask
end

m.coast=[];

    
[demdir,name,ext] =fileparts([(ifile)]);
mon=str2num(name(10:11));
name1=[name(end-3:end),ext];
 
if isempty(ndwisv.z) %no saved NDWI matrix. 
   % To check whether the data is mono or scene files
   
  if strcmp(name1,'meta.txt')
	flagfmt=3; %stereo files
    % file names to be used
    mtFile = strrep(ifile,'meta.txt','matchtag.tif');
    orFile = strrep(ifile,'meta.txt','ortho.tif');
    metaFile = strrep(ifile,'meta.txt','meta.txt');
    demFile= strrep(ifile,'meta.txt','dem.tif');

  elseif strcmp(ext,'.xml')
    flagfmt=1; %xml files mono
    orFile = strrep(ifile,'.xml','.tif');
    metaFile = strrep(ifile,'.xml','.xml');
  end

if mon <=4 || mon>=11 ; frozen =6;end
frozen=7;%0; % Control flag to choose different algorithms.

% There is some sort of scaling difference between the wv_correct and non
% wv_correct - applied othos that screws up the conversion to uint8 in
% entropyfilt (which uses im2uint8) unless I convert to uint8 first. So
% need to check if wvc applied.
if 0 % convert to uint8 always
wvcFlag=false;

% check if meta file exists
if ~exist(metaFile,'file')
   warning('no meta file found, assuming no wv_correct applied\n');
else
    % exists, so read meta file
    c=textread(metaFile,'%s','delimiter','\n');
    
    %find SETSM version
    str='Image 1 wv_correct=';
    r=find(~cellfun(@isempty,strfind(c,str)));
    if isempty(r) % added by chunli for some missing data files.
        warning('no meta file found, assuming no wv_correct applied\n');
    else
    value=deblank(strrep(c{r(1)},str,''));
    wvcFlag = logical(str2num(value));
    
    if wvcFlag
        fprintf('wv_correct applied\n')
    else
        fprintf('wv_correct not applied\n');
    end
    end
end
end

%% Data density image

% read the matchtag and ortho files
if 0 % not used
mt=readGeotiff(mtFile,'map_subset', rangeov);

res=mt.info.map_info.dx;
resr=2.;dsr=res/resr;
% dsr2=0.1;resr2=res/dsr2; %for faster computation 
dsr2=0.25;resr2=res/dsr2; %for faster computation 

% calculate data density map using default (11x11 kernel)
P = DataDensityMap(mt.z);

% resize to 8m - resizing the mt directly results in point loss
P = imresize(P,dsr);

Mden=(P > 0.9);
end

%% Entropy Image

% read ortho
or=readGeotiff(orFile,'map_subset', rangeov);
res=or.info.map_info.dx;
resr=2.;dsr=res/resr;
dsr2=0.25;resr2=res/dsr2; %for faster computation 
dsr2=1;resr2=res/dsr2; %for better resolution

mt.x=or.x;mt.y=or.y; mt.z=ones(size(or.z)); P=ones(size(or.z));

orsv=or;
or = or.z; % get rid of map data

% resize ortho to 8m  -> To check difference between 8m and 2 m
or = imresize(or,dsr2);
% or=uint16(or); %chunli 2017/6/7
or=double(or); %chunli 11/21/2018

% subtraction image
or_subtraction =  movmax(or,5) - movmin(or,5);
 
or_subtraction(imdilate(or == 0,ones(9))) = 0;

if 0
    if ~wvcFlag; or_subtraction = uint8(or_subtraction); end %
else
    % in entropyfilt.m, it uses im2uint8, which for input of uint16, it is
    % divided by 257, rounding, and then casting to uint8
    M=~imdilate(or == 0,ones(9));
    I=or_subtraction;I(~M)=0;
    min1=min(I(M));max1=max(I(M));
    I(M)=((I(M) - min1) / (max1 - min1)*255);
    or_subtraction=uint8(I);
    clear I
end

 % entropy image
J = entropyfilt(uint8(or_subtraction),ones(5)); %or_subtraction best to be uint8
J = imresize(J,size(P),'nearest');
or = imresize(orsv.z,size(P),'nearest');
clear orsv or_subtraction

% Get the reflectance 
data.z=or;
or=DN2reflpan(data,ifile);

%% Entropy & Density Mask

%% Sigma Slope & Density Mask

if flagfmt==3&&0 %DEM exist
% read dem
z = readGeotiff(demFile,'map_subset', rangeov);

% x = z.x; only needed for hillshade 
% y = z.y; only needed for hillshade
z = z.z;

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
Mrho=(rhosd<0.5);    

else %DEM not exist
        Mrho=true(size(P));
end

% sigma slope threshold function
% rhosdthresh= (1.*P).^4;

% apply theshold
% M1 = rhosd <= rhosdthresh;
if  frozen==2  % first test on WV01_20161017_102001005ACECD00_1
M = (J >3)  ; 
M1 = rhosd >= 0.25; % water tiles
M1 = bwareaopen(M1,1000*dsr);
M=M|M1;
Md1 = imdilate(M, ones(3));
M=logical(Md1-M);
M=M&Mden;
elseif frozen==1 %frozen river, test
%Entropy Mask
M = ( J <1)& ~(imdilate(or == 0,ones(round(9*2*8/resr))));
%M = bwareaopen(M, 1000*5); %Use Entropy only
M = bwareaopen(M, 1000*50*dsr2); % Use Entropy to get statistics for DEM and ortho

meandem=mean(z(M));stddem=std(z(M));
Mdem=z>meandem-stddem&z<meandem+stddem;

meanor=mean(double(or(M)));stdor=std(double(or(M)));
Mor=or>meanor-3*stdor;

Mod = bwareaopen(Mor&Mdem, 1000*10);

% get a rough boundary for river path
cf = 0.1;%0.5; %boundary curvature factor (0= point boundary, 1 =conv hull)
MJr = imresize(M,0.1);
[idy,idx]=find(MJr==1);
k = boundary(idx,idy,cf);
[nyov,nxov]=size(MJr);
Mb = poly2mask(idx(k),idy(k),nyov,nxov);
MJbd = imresize(Mb,size(mt.z),'nearest');
% expand to ensure boundary/coast coverage to the width of a river 400m
MJbd = imdilate(MJbd, ones(200)); 

Modj=(Mor&Mdem&MJbd);
Modj= bwareaopen(Modj, 1000*1000);
Modfil = bwareaopen(~Modj, 1000*5);
Modfil=~Modfil;
M=Modfil;

elseif frozen==0 %running Blyde River, works the  best.
% M=(J>1)&(rhosd<0.3)&(or>4000); %running water with clouds
% M=(J>1)&(rhosd<0.2)&(or>4000); %running water with clouds, 2m res %Clear city,WV02_20160725_1030010059A50600_103001005ABCBF00_500849896010_01_P001_500852789090_01_P001_2_meta
% M = bwareaopen(M, 1000*500);
% M=(J>2)&(rhosd<0.5)&(or>9000); %WV01_20161017 2m
% M = ( J <1)&(or>max(or(:))*0.3); %WV02_20160904_; res,res2=[2, 8m]
% M=(J>1)&(rhosd<0.5)&(or>max(or(:))*0.2); %WV02_20161119
% M=(J<2)&(rhosd<0.5);%running dark river; 20130526
M=(J<0.5)&(rhosd<0.5);%running dark river; 20111008
% M = bwareaopen(M, 1000*5);
M = bwareaopen(M, 1000*200);%20121011

meanor=mean(double(or(M)));stdor=std(double(or(M)));
% Mor=(or>meanor-stdor)&(or<meanor+stdor); %running dark river; 20130526
Mor=(or>meanor-3*stdor)&(or<meanor+3*stdor); %20121011

cf = 0.1;%0.5; %boundary curvature factor (0= point boundary, 1 =conv hull)
MJr = imresize(M,dsr2);
[idy,idx]=find(MJr==1);
k = boundary(idx,idy,cf);
[nyov,nxov]=size(MJr);
Mb = poly2mask(idx(k),idy(k),nyov,nxov);
MJbd = imresize(Mb,size(mt.z),'nearest');
% expand to ensure boundary/coast coverage to the width of a river 400m
MJbd = imdilate(MJbd, ones(200));

% Modj=(Mor&MJbd);
Modj=(Mor);%20121011
% Modj= bwareaopen(Modj, 1000*1000);
% Modj= bwareaopen(Modj, 1000*5);
% Modj= bwareaopen(Modj, 1000*50); %running dark river; 20130526
Modj= bwareaopen(Modj, 1000*500); %running dark river; 20111008
Modfil = bwareaopen(~Modj, 1000*5);
Modfil=~Modfil;

Md1 = imdilate(Modfil, ones(3));
M=logical(Md1-Modfil);
% M=M&Mden&~(imdilate(or == 0,ones(round(9*2*8/resr))));
M=M&~(imdilate(or == 0,ones(round(9*2*8/resr)))); %%running dark river; 20130526 coastal forest

elseif frozen==3 %i=2681;%20160514 %  Yukon
%     M = ( J <3)& ~(imdilate(or == 0,ones(round(9*2*8/resr))));
% M = bwareaopen(M, 1000*500*dsr2);
%i=1328;%20160812
M = ( J <0.4)& ~(imdilate(or == 0,ones(round(9*2*8/resr))));
M = bwareaopen(M, 1000*300*dsr2);

Modj=M;
Modfil = bwareaopen(~Modj, 1000*50);
Modfil=~Modfil;

Md1 = imdilate(Modfil, ones(3));
M=logical(Md1-Modfil);
M=M&Mden&~(imdilate(or == 0,ones(round(9*2*8/resr))));
elseif frozen ==4 % i=2710%  Yukon Branches

M = ( J <1.5)& ~(imdilate(or == 0,ones(round(9*2*8/resr)))); %i=2710 i=3314;
% M = ( J <2.5)& ~(imdilate(or == 0,ones(round(9*2*8/resr)))); %i 2719
M = bwareaopen(M, 1000*5);

meanor=mean(double(or(M)));stdor=std(double(or(M)));
Mor=(or>meanor-stdor)&(or<meanor+3*stdor);
% mp1=roipoly;Mor(~mp1)=0;% i=3314;

Modj=Mor;
% Modj= bwareaopen(Modj, 1000*1000);
Modj= bwareaopen(Modj, 1000*5);%i 2733
% Modj= bwareaopen(Modj, 1000*50);%i 3314
Modfil = bwareaopen(~Modj, 1000*5);
Modfil=~Modfil;

Md1 = imdilate(Modfil, ones(3));
M=logical(Md1-Modfil);
M=M&Mden&~(imdilate(or == 0,ones(round(9*2*8/resr))));

elseif frozen==5 %Use prior river boundary
    
    if exist('TananaUB.mat') && 0
        load('TananaUB.mat'); %mb.x, mb.y mb.z
        MJbd = interp2(mb.x,mb.y,mb.z,mt.x,mt.y','*nearest',0);
        MJbd = imresize(MJbd,size(P),'nearest');
    elseif exist('TananaUB2.mat')
	load('TananaUB2.mat') %rivns; x,y in meter
	[X,Y]=meshgrid(mt.x,mt.y);
	MJbd = inpolygon(X,Y,rivns(:,1),rivns(:,2));
    end
        [hpdf,edges]=histcounts(or(MJbd&or~=0),'Normalization','pdf');
%         hpdf(edges(2:end)<max(or(:))*0.2)=0;
       [hmax,idm]=max(hpdf);
        meanor=(edges(idm)+edges(idm+1))/2;
        hh=hmax*0.2;
        id2=idm+1:length(hpdf);
        id2i=find((hpdf(id2)-hpdf(id2-1))>0&hpdf(id2)<hmax/2.);%find the curve turning point
        if ~isempty(id2i); id2=id2(1):id2(id2i(1));end        
        [~,id2i]=min(abs(hpdf(id2)-hh));
        hw2=id2(id2i(1));
        
        id1=1:min(idm-1,length(hpdf));
        id1i=find((hpdf(id1+1)-hpdf(id1))<0&hpdf(id1)<hmax/2.);%find the curve turning point
        if ~isempty(id1i); id1=id1(id1i(end)):id1(end);end        
        [~,id1i]=min(abs(hpdf(id1)-hh));
        hw1=id1(id1i(end));
   %    hw=edges(hw2)-edges(hw1);
    if 0;figure;plot(edges(2:end),hpdf)
    hold on;plot([edges(hw1+1),edges(hw2+1)],[hh,hh],'r-')
    end

Mor=(or>edges(hw1+1))&(or<edges(hw2+1));

Modj=(Mor&MJbd);
% Modj= bwareaopen(Modj, 1000*1000);
Modj= bwareaopen(Modj, 1000*20);
Modfil = bwareaopen(~Modj, 1000*5);
Modfil=~Modfil;

Md1 = imdilate(Modfil, ones(3));
M=logical(Md1-Modfil);
M=M&~(imdilate(or == 0,ones(round(9*2*8/resr))));

elseif frozen==6 %centerline +-10meter %for frozen river
    rivw=20;%meter
    if exist('tancl.mat')
	load('tancl.mat') % ct nx2;lon,lat
	[ctx,cty]=polarstereo_fwd(ct(:,2),ct(:,1),[],[],70,-45);
idl=isnan(ctx)|isnan(cty);
ctx(idl)=[];cty(idl)=[];
%Revisit for multipile functions
ctyi=mt.y;ctxi = interp1(cty,ctx,ctyi); %interpolate the centerlines
ctx=ctxi;cty=ctyi;
idl=isnan(ctx)|isnan(cty);
ctx(idl)=[];cty(idl)=[];

idx=round((ctx-mt.x(1))/res)+1;
idy=round(-(cty-mt.y(1))/res)+1;
[ny,nx]=size(mt.z);
idl=idx<1|idx>nx|idy<1|idy>ny;
idx(idl)=[];idy(idl)=[];
id=(idx-1)*ny+idy;
MJbd= zeros(size(mt.z));
MJbd(id)=1;
        MJbd = imresize(MJbd,size(P),'nearest');
    else
	warning('Centerline file not found.')
    end

Modfil=imdilate(MJbd,ones(ceil(rivw/resr)));

M=Modfil&~(imdilate(or == 0,ones(round(9*2*8/resr))));

elseif frozen==7
%     &Mrho;
M=(J<0.5)&~(imdilate(or <= 0,ones(round(9*2*8/resr))));%running dark river; 20111008
clear J
% M = bwareaopen(M, 1000*5);
% M = bwareaopen(M, 1000*200);%20121011

meanor=nanmean(double(or(M)));stdor=nanstd(double(or(M)));
% Mor=(or>meanor-stdor)&(or<meanor+stdor); %running dark river; 20130526

end
clear P

or(or<=0)=nan;%set the image edge of M to -1.
m.n = imresize(or,size(mt.z),'nearest');
m.x = mt.x;
m.y = mt.y;

else % use the saved 

or=ndwisv.z; mt.x=ndwisv.x;mt.y=ndwisv.y;mt.z=ones(size(or));
m.x=mt.x;m.y=mt.y;m.n=or; %use the input coordinates
end


if step==1 %step 1, use the general threshold. Notice: cannot work using the saved ndwisv.
badflag=0;   
% Mor=(or>meanor-3*stdor)&(or<meanor+3*stdor); %20121011
Mor=(or<meanor+stdor); %20121011

elseif step ==2   %step 2, use the refined threshold.
%retrieving statistics of NDWI over a priori water mask.
% wm: 1 water, 0 land, -1 void
resrc=abs(wm.x(2)-wm.x(1));
% Md1 = imdilate(wm.z, ones(widthstat/resrc)); %make sure water is absolutely water.
Md1=wm.z==1;
%Remove small lakes, which may have different reflectance than river water
%Md1= bwareaopen(wm.z, round(lakearea/resrc/resrc)); %remove small clusters

wm1 = interp2(wm.x,wm.y,Md1,mt.x,mt.y','*nearest',0); %out of region: 0, non water.
ndwil=double(or(wm1==1)); %water
mean1=nanmean(ndwil);std1=nanstd(ndwil);cnt1=sum(~isnan(ndwil));

Md1 = imdilate(wm.z~=0, ones(600/resrc)); %make sure land is absolutely land. 300 m width expansion
%Me1=imerode(wm.z==0, ones(widthstat/resrc)); %land
wm1 = interp2(wm.x,wm.y,Md1,mt.x,mt.y','*nearest',1); %out of region: 1,  water.
ndwil=double(or(wm1==0)); %land
mean2=nanmean(ndwil);std2=nanstd(ndwil);cnt2=sum(~isnan(ndwil));
clear wm wm1

ithreshold=0;
badflag=0;
if (isnan(mean1)||cnt1<cntmin) 
    badflag=1;
else
meanor=mean1;stdor=std1;
Mor=(or>meanor-stdor)&(or<meanor+stdor); %20121011
end

fid = fopen('./panstats.dat','a'); %statistics
fprintf(fid,' %f %f %f %f %f %d (NDWI mean std thres) %s \n',mean1,std1,mean2,std2,ithreshold,badflag,ifile);
fclose(fid);

end

if badflag==1    
Mstrip=struct(); Mstrip.x=[];Mstrip.y=[];Mstrip.z=[];Mstrip.coast=[];Mstrip.n=[];
m=Mstrip;
return
end

M=int8(Mor);
%M(or==0)=-1;%set the image edge of M to -1.
%M(or<=0)=-1;%set the image edge of M to -1. %Reflectance
M(isnan(or))=-1;%set the image edge of M to -1. %Reflectance

if 0
% Modj=(Mor&MJbd);
Modj=(Mor);%20121011
% Modj= bwareaopen(Modj, 1000*1000);
% Modj= bwareaopen(Modj, 1000*5);
% Modj= bwareaopen(Modj, 1000*50); %running dark river; 20130526
Modj= bwareaopen(Modj, 1000*500); %running dark river; 20111008
Modfil = bwareaopen(~Modj, 1000*5);
Modfil=~Modfil;

Md1 = imdilate(Modfil, ones(3));
M=logical(Md1-Modfil);
M=M&~(imdilate(or <= 0,ones(round(9*2*8/resr)))); %%running dark river; 20130526 coastal forest
end
clear or

m.z = imresize(M,size(mt.z),'nearest');

% m.info = mt.info;
% m.Tinfo = mt.Tinfo;

Mstrip=m;
if 1
%coast output not used.
Medgs=isnan(Mstrip.z(:,:));%(Mstrip.z==-1);%isnan(Mstrip.z(:,:));
Med=imdilate(Medgs,ones(4));
Medgs=Med; clear Med
Modj=Mstrip.z;Modj(Medgs)=0;
whos %check memory usage
resx=mean(Mstrip.x(2:end)-Mstrip.x(1:end-1));resy=mean(Mstrip.y(2:end)-Mstrip.y(1:end-1));resr=mean([abs(resx),abs(resy)]);
Modj= bwareaopen(Modj, round(lakearea/resr/resr)); %remove small clusters
Modfil = bwareaopen(~Modj, round(cloudarea/resr/resr)); %fill small areas: 1e4*4m^2
Modfil=~Modfil;

Md1 = imdilate(Modfil, ones(3));
M=logical(Md1-Modfil);
% M=M&~(imdilate(data.z(:,:,1) == 0,ones(4))); 
M=M&~Medgs; 
end %if 0

if 0 %flagplot==1
[X,Y]=meshgrid(Mstrip.x,Mstrip.y);
figure;
hold all
set(gcf,'Color','white');set(gca,'FontSize', 18);set(gcf, 'PaperPosition', [0 0 6 4]);set(gcf, 'PaperSize', [ 6 4]);
imagesc(Mstrip.x,Mstrip.y,double(Mstrip.n))
colormap gray
hold on;plot(X(M),Y(M),'ro')
title(['Ortho Pan-band ',name(6:13)]);
axis equal
%saveas(gcf,[name],'fig') %too big
print('-dpdf','-r400',name)
end

m.coast=M;%[X(M),Y(M)]; %merged coast for the strip

close all
return
end

 
