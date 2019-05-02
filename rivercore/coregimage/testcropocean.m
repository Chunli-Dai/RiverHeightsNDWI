% test mimc2, using a priori coastline +- 2km band zone only for feature tracking.
%see also CoastTilev2bp4.m
macdir='/Users/chunlidai/surge/';
addpath(genpath(['../rivergithub2/coregimage/']))
resr=2;resrc=40;width=2e3;
orFile1=['WV02_20150412_10300100400AF500_1030010040C40400_seg1_2m_orthos1668.tif'];
orFile2=['WV02_20140925_1030010036BD9500_103001003727E300_seg6_2m_orthos1668.tif'];

orFile2=['WV02_20150412_10300100400AF500_1030010040C40400_seg1_2m_ortho.tif'];
orFile1=['WV02_20140925_1030010036BD9500_103001003727E300_seg6_2m_ortho.tif'];

i=16680; %croped

or1=readGeotiff(orFile1);
or2=readGeotiff(orFile2);

rangtar=[min(or1.x) max(or1.x) min(or1.y) max(or1.y)];
rangref=[min(or2.x) max(or2.x) min(or2.y) max(or2.y)];
rangeov=[max(rangtar(1),rangref(1)),min(rangtar(2),rangref(2)), max(rangtar(3),rangref(3)),min(rangtar(4),rangref(4))];
ranget=rangeov;
ranget=round(ranget/resrc)*resrc;

%get the lat lon of rangeov.
x0=[ranget(1) ranget(2) ranget(2) ranget(1) ranget(1) ];y0=[ranget(4) ranget(4) ranget(3) ranget(3) ranget(4) ];
[lat0,lon0]=polarstereo_inv(x0,y0,[], [],70,-45);
maxlat=max(lat0);minlat=min(lat0);
maxlon=max(lon0);minlon=min(lon0);

tic
shpname='./codec2/GSHHS/GSHHS_f_L1.shp';
S = shaperead(shpname);
cnt=length(S);

ids=[];
for j=1:cnt
    if any(S(j).Y >= minlat & S(j).Y<= maxlat & S(j).X >= minlon & S(j).X <= maxlon )
        ids=[ids,j];
    end
end

xw=ranget(1):resrc:ranget(2);yw=ranget(4):-resrc:ranget(3); % a priori water mask
nwx=length(xw);nwy=length(yw);
smg=false(nwy,nwx);%water mask from a priori coastline shapefiles

for k=1:length(ids)
    j=ids(k);
    [sx,sy]=polarstereo_fwd(S(j).Y,S(j).X,[], [],70,-45);
    
    dx=resrc;
    idx=round((sx-ranget(1))/dx)+1;idy=round((sy-ranget(4))/(-dx))+1;
    idt=isnan(idx)|isnan(idy);
    idx(idt)=[];idy(idt)=[];
    sm=poly2mask(idx,idy,nwy,nwx); % fast, apply to each polygon one by one.
    smg=smg|sm;
end
wm=[];wm.x=xw;wm.y=yw;wm.z=smg;clear smg; %1 land, 0 water
if(sum(wm.z(:))==0||sum(wm.z(:))==nwx*nwy);fprintf('This tile contain no coastline (all land or all ocean).');return;end % if all land or all water, i.e. no coastline.
% get the buffer region of coastline: +- width of the coastline;
Md1 = imdilate(wm.z, ones(width/resrc));
Me1=imerode(wm.z, ones(width/resrc)); % erode the box edge also; to fix: change width in CoastTileMain;
Mcb=Md1&~Me1; %coastal band
wm.z=Mcb;
% figure;imagesc(wm.x,wm.y,Mcb);colorbar %plot the coastal band

[infile1,infile2]=preparemimc2(orFile1, orFile2,i,rangeov,wm);

infile1=['./data/WV02_20140925_1030010036BD9500_103001003727E300_seg6_2m_orthos2e.tif'];
infile2=['./data/WV02_20150412_10300100400AF500_1030010040C40400_seg1_2m_orthos2e.tif'];

infile2=['./WV02_20150412_10300100400AF500_1030010040C40400_seg1_2m_orthos1668.tif'];
infile1=['./WV02_20140925_1030010036BD9500_103001003727E300_seg6_2m_orthos1668.tif'];

infile2='./WV02_20160827_103001005A735B00_103001005CAE7900_seg1_2m_ortho_15m.tif';
infile1='./WV02_20100710213820_1030010005A70000_10JUL10213820-M1BS-052618580140_01_P004_u16ns3413_15m.tif';
infile2='./WV02_20160827214547_103001005CAE7900_16AUG27214547-M1BS-500890307070_01_P009_u16ns3413_15m.tif';

[dx,dy, sigma]=vmimc(infile1,infile2);
% [dx,dy, sigma]=vmimc('./data/t1.tif','./data/ref.tif'); %Bug: name has to start with ymd




%% Coregistration
orFile2=['WV02_20150412_10300100400AF500_1030010040C40400_seg1_2m_ortho.tif'];
orFile1=['WV02_20140925_1030010036BD9500_103001003727E300_seg6_2m_ortho.tif'];

infile=strrep(orFile1,'ortho.tif','dem.tif');
data=readGeotiff(infile,'map_subset', rangeov);
%   reduce resolution
ranget=[[min(data.x) max(data.x) min(data.y) max(data.y)]/resrc];
ranget=[ceil(ranget(1)) floor(ranget(2)) ceil(ranget(3)) floor(ranget(4))]*resrc;
tx=ranget(1):resrc:ranget(2);ty=ranget(4):-resrc:ranget(3);
res=data.info.map_info.dx;
nsr=resrc/res;
if nsr <1 
 data.z(data.z == -9999) = NaN; % %convert nodata values to nan for interpolation
 tz =interp2(data.x,data.y,double(data.z),tx,ty','*linear');%toc;%4s
 tz(isnan(tz))=-9999; %return to -9999
%the following get the exactly the same as from interp2
else
idrxs=find(abs(data.x-ranget(1))<1e-3);idrxe=find(abs(data.x-ranget(2))<1e-3);
idrys=find(abs(data.y-ranget(4))<1e-3);idrye=find(abs(data.y-ranget(3))<1e-3);
idrx=idrxs:nsr:idrxe;
idry=idrys:nsr:idrye;
dd=[idrx(2:end)-idrx(1:end-1),idry(2:end)-idry(1:end-1)];
if isempty(dd)||any(abs(dd-nsr)>1e-3);warning('Resized grid is not constant spacing.');end
if length(tx)~=length(idrx)||length(ty)~=length(idry);warning('Wrong downsizing of DEM.');end
tz=data.z(idry,idrx); %0.09s
end
datar= struct();
datar.x=tx;datar.y=ty;  datar.z=tz;
data0r=datar;

infile=strrep(orFile2,'ortho.tif','dem.tif');
    data=readGeotiff(infile,'map_subset', rangeov);
%   reduce resolution
    ranget=[[min(data.x) max(data.x) min(data.y) max(data.y)]/resrc];
    ranget=[ceil(ranget(1)) floor(ranget(2)) ceil(ranget(3)) floor(ranget(4))]*resrc;
    tx=ranget(1):resrc:ranget(2);ty=ranget(4):-resrc:ranget(3);
    res=data.info.map_info.dx;
    nsr=resrc/res;
    if nsr <1 
     data.z(data.z == -9999) = NaN; % %convert nodata values to nan for interpolation
     tz =interp2(data.x,data.y,double(data.z),tx,ty','*linear');%toc;%4s
     tz(isnan(tz))=-9999; %return to -9999
    %the following get the exactly the same as from interp2
    else
    idrxs=find(abs(data.x-ranget(1))<1e-3);idrxe=find(abs(data.x-ranget(2))<1e-3);
    idrys=find(abs(data.y-ranget(4))<1e-3);idrye=find(abs(data.y-ranget(3))<1e-3);
    idrx=idrxs:nsr:idrxe;
    idry=idrys:nsr:idrye;
    dd=[idrx(2:end)-idrx(1:end-1),idry(2:end)-idry(1:end-1)];
    if isempty(dd)||any(abs(dd-nsr)>1e-3);warning('Resized grid is not constant spacing.');end
    if length(tx)~=length(idrx)||length(ty)~=length(idry);warning('Wrong downsizing of DEM.');end
    tz=data.z(idry,idrx); %0.09s
    end
    datar= struct();
    datar.x=tx;datar.y=ty;  datar.z=tz;

    rangtar=[min(datar.x) max(datar.x) min(datar.y) max(datar.y)];
    rangref=[min(data0r.x) max(data0r.x) min(data0r.y) max(data0r.y)];
    rangeov=[max(rangtar(1),rangref(1)),min(rangtar(2),rangref(2)), max(rangtar(3),rangref(3)),min(rangtar(4),rangref(4))];
    %crop the overlapping data
    refdem=[];tardem=[];
    idx=find(data0r.x>=rangeov(1) & data0r.x<=rangeov(2));
    idy=find(data0r.y>=rangeov(3) & data0r.y<=rangeov(4));
    refdem.x=data0r.x(idx);refdem.y=data0r.y(idy);
    refdem.z=data0r.z(idy,idx);
    idx=find(datar.x>=rangeov(1) & datar.x<=rangeov(2));
    idy=find(datar.y>=rangeov(3) & datar.y<=rangeov(4));
    tardem.x=datar.x(idx);tardem.y=datar.y(idy);
    tardem.z=datar.z(idy,idx);
    [n,m]=size(tardem.z);minm=min(n,m);
    if sum(size(tardem.z)~=size(refdem.z))||minm<3
        warning(['Wrong overlapping crop, check i:',num2str(idregion(j))]);
        idd=[idd;j];
    end

      datatarz=tardem.z;datatarz(datatarz== -9999) = NaN; 
      datarefz=refdem.z;datarefz(datarefz== -9999) = NaN;
      iter= 1;
    %   [dx,dy,sigma]=vmimc(infile1,infile2) 
      [z2out,p,sigma] = coregisterdems(refdem.x,refdem.y,double(datarefz),tardem.x,tardem.y,double(datatarz));
        rmsreg2(j)=sigma;
      if sum(size(z2out)~=size(refdem.z)) || sigma>10 || isnan(sigma) %control parameter
          warning(['coregistration failure',infile]); p=zeros(3,1);
          idd=[idd;j];
      else
        dx=p(2);dy=p(3); %z, x, y
      end
