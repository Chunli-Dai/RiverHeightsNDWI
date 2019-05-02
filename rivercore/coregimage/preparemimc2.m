function [OutName1,OutName2]=preparemimc2(varargin)
% Usage [OutName1,OutName2]=preparemimc2(orFile1, orFile2,rangeov);
% read files, crop common regions, write tiff files.
%orFile1, orFile2,rangeov
%change the resoluton to be 15m to fit mimc2.
resm=15;%mimc2 prefered image resolution.
if nargin==3
    orFile1=varargin(1);orFile2=varargin(2);
    i=varargin{3};
    or1=readGeotiff(orFile1);
    or2=readGeotiff(orFile2);
    
    rangtar=[min(or1.x) max(or1.x) min(or1.y) max(or1.y)];
    rangref=[min(or2.x) max(or2.x) min(or2.y) max(or2.y)];
    rangeov=[max(rangtar(1),rangref(1)),min(rangtar(2),rangref(2)), max(rangtar(3),rangref(3)),min(rangtar(4),rangref(4))];
elseif nargin==4
    orFile1=varargin{1};orFile2=varargin{2};
    i=varargin{3};
    rangeov=cell2mat(varargin(4));
    or1=readGeotiff(orFile1,'map_subset', rangeov);
    or2=readGeotiff(orFile2,'map_subset', rangeov);
    
elseif nargin==5
    orFile1=varargin{1};orFile2=varargin{2};
    i=varargin{3};
    rangeov=cell2mat(varargin(4));
    wm=cell2mat(varargin(5));
    or1=readGeotiff(orFile1,'map_subset', rangeov);
    or2=readGeotiff(orFile2,'map_subset', rangeov);

end

if 1 %use mask in dem strip files
    infile1= strrep(orFile1,'ortho.tif','dem.tif');
    dem1=readGeotiff(infile1,'map_subset', rangeov);
    or1.z(dem1.z==-9999)=0;
    or2=readGeotiff(orFile2,'map_subset', rangeov);
    infile2= strrep(orFile2,'ortho.tif','dem.tif');
    dem2=readGeotiff(infile2,'map_subset', rangeov);
    or2.z(dem2.z==-9999)=0;
end

outdir='./work/';
if ~exist(outdir,'dir')
    mkdir outdir
end

%crop the overlapping data
refdem=[];tardem=[];
idx=find(or2.x>=rangeov(1) & or2.x<=rangeov(2));
idy=find(or2.y>=rangeov(3) & or2.y<=rangeov(4));
refdem.x=or2.x(idx);refdem.y=or2.y(idy);
refdem.z=or2.z(idy,idx);
idx=find(or1.x>=rangeov(1) & or1.x<=rangeov(2));
idy=find(or1.y>=rangeov(3) & or1.y<=rangeov(4));
tardem.x=or1.x(idx);tardem.y=or1.y(idy);
tardem.z=or1.z(idy,idx);

% if two images have the same grids, leave it, if not the same grids, interpolate
flag1=0;
if length(tardem.x)==length(refdem.x) && length(tardem.y)==length(refdem.y)
df=max(abs([tardem.x(:)-refdem.x(:);tardem.y(:)-refdem.y(:)]));
if df < 1e-9
    flag1=1; % two imges have the same grid
end
end
if flag1==0 
    fprintf ('Two images are not on same grids for feature tracking, do interpolation.');
    z1=interp2(refdem.x,refdem.y,double(refdem.z),tardem.x,tardem.y','*linear',0);
    refdem.x=tardem.x;refdem.y=tardem.y;refdem.z=z1; %image 2 on image 1's grid
end

% to 15 meter resolution
x15=min(tardem.x):resm:max(tardem.x);
y15=max(tardem.y):-resm:min(tardem.y);
z1=interp2(tardem.x,tardem.y,double(tardem.z),x15,y15','*linear',0);
tardem.x=x15;tardem.y=y15;tardem.z=z1; %image 1
z1=interp2(refdem.x,refdem.y,double(refdem.z),x15,y15','*linear',0);
refdem.x=x15;refdem.y=y15;refdem.z=z1; %image 2 
    
taredge=tardem.z==0;refedge=refdem.z==0; %save the edge
comedge=taredge|refedge;

% transform DN to reflectance for panchromatic band only.
if 1
for k=1:2
    if k==1
        metafile=strrep(orFile1,'ortho.tif','meta.txt');
    elseif k==2
        metafile=strrep(orFile2,'ortho.tif','meta.txt');
    end
[~,filename,~]=fileparts(metafile);
%assume in each strip, all scenes and each image in pairs have the same
%effectivebandwith and abscalfactor.
c=textread(metafile,'%s','delimiter','\n');
%strg={['_abscalfact='],['_effbw='],['_Mean_sun_elevation=']};
strg={['abscalfact'],['effbw'],['Mean_sun_elevation']};
for j=1:3
str=strg(j);
r=find(~cellfun(@isempty,strfind(c,str)));
str=c{r(1)};r1=strfind(str,'=');str(1:r1)='';
% Xbs=deblank(strrep(c1{r(1)},str,''));
t1(j)= sscanf(str, '%g', 1);
end
abscalfactor=t1(1);effectivebandwith=t1(2);meanSunEl=t1(3);
Theta=(90.-meanSunEl)*pi/180.;
% sun-earth distance polynomial function coefficients
% doy = day of the year: 
year=sscanf(filename(6:9), '%g', 1); month=sscanf(filename(10:11), '%g', 1); day=sscanf(filename(12:13), '%g', 1);
doy=juliandate(year,month,day)-juliandate(year,1,1)+1;
C = [1.8739e-26,-3.4455e-23,2.7359e-20,-1.2296e-17,3.0855e-15,-2.2412e-13,-5.8744e-11,6.9972e-10,2.5475e-06,-1.6415e-05,0.9833];
dES = polyval(C,doy); % earth-sun distance

satname=filename(1:4);
[GAIN,OFFSET,Esun]=readgainoffset(satname); %read Gain OFFset data, GainOffset.txt.
if k==1
DN=double(tardem.z(:,:));
L=GAIN(end)*DN*(abscalfactor/effectivebandwith)+OFFSET(end);
rho=L*dES^2*pi/(Esun(end)*cos(Theta));
tardem.z=rho;
elseif k==2
DN=double(refdem.z(:,:));
L=GAIN(end)*DN*(abscalfactor/effectivebandwith)+OFFSET(end);
rho=L*dES^2*pi/(Esun(end)*cos(Theta));
refdem.z=rho;
end
end %k=1:2
end % if transform DN to reflectance

%remove snow 
mp=tardem.z>10|refdem.z>10.;
comedge=comedge|mp;

%apply ocean map; keeping only coastal band zone. 
if 0
   mp=interp2(wm.x,wm.y,double(wm.z),tardem.x,tardem.y','*nearest',0); %coastal band
   comedge=~logical(mp)|comedge;% more edges after applying the coastal map
end

%apply the comedge
tardem.z(comedge)=nan;  refdem.z(comedge)=nan;

idx=find(any(comedge==0,1));
idy=find(any(comedge==0,2));
idx=min(idx):max(idx);idy=min(idy):max(idy);
tardem.x=tardem.x(idx);tardem.y=tardem.y(idy);tardem.z=tardem.z(idy,idx);
refdem.x=refdem.x(idx);refdem.y=refdem.y(idy);refdem.z=refdem.z(idy,idx);

if or1.Tinfo.GeoDoubleParamsTag(1) > 0
    projstr='polar stereo north';
else
    projstr='polar stereo south';
end

%OutName1=strrep(orFile1,'ortho.tif',['orthos',num2str(i),'.tif']);
OutName1=strrep(orFile1,'ortho.tif',['orthos',num2str(i),'.tif']);
[dir,name,ext] =fileparts(OutName1);
OutName1=[outdir,name,ext];
%writeGeotiff(OutName1,tardem.x,tardem.y,uint16(tardem.z),12,0,projstr)
writeGeotiff(OutName1,tardem.x,tardem.y,(tardem.z),5,0,projstr)

% OutName2=strrep(orFile2,'ortho.tif','orthos.tif');
OutName2=strrep(orFile2,'ortho.tif',['orthos',num2str(i),'.tif']);
[dir,name,ext] =fileparts(OutName2);
OutName2=[outdir,name,ext];
%writeGeotiff(OutName2,refdem.x,refdem.y,uint16(refdem.z),12,0,projstr)
writeGeotiff(OutName2,refdem.x,refdem.y,(refdem.z),5,0,projstr)

return
end

