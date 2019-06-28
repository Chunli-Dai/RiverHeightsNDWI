function [Co]=ChangeRiver(igage,source,texteq,range,XYbg,f,fdir,range2,XYbg2,f2,fdir2)
% Change detection
% v1 
% v2 faster computation through the overlapping of real data boundary
% v3 data cover 10% more of the subzone; reference dem cover 50% more of the subzone
%  parallel
% v4 use ArcticDEM mosaic files for reference DEM;
%      using *reg.txt for coregistration parameters
% won't work
% v5 don't do coregistration, adjust height using the average over rock 
%  v5 bp1: 
%  v5 bp2: using Ian's coregistraton code with rock controls
%  v5 bp3: all in one with options;Do not apply filtering for faster
%     bp4: Store DEM;-> reduce time from 1 hour to 30 minutes
%	  Interpolate rock on reference one time; ->if skip the interpolation,
%         for the coregistration step, reduce from 6 sec to 3 sec,
%	  so overall it doesn't matter much.
%     bp5: MJ's coregistration, 1\ crop only the overlapping area -> 30m-> 25 minutes
%	  2\ use the transformation parameter with only rock surfaces
% v6: use the actuall overlapping
%   bp1: 
%   bp2: avoid 'find' for searching overlapping polygons for faster computation
%	 in coregflag 1, replace -9999 to NaN to avoid processing -9999+meanhrock
%	 for searching reference, has to have rock
%	 add filter to exclude bad dems --> to be realized
%   bp3: avoid the repeated calculation from different groups of data
%   Copied from ChangedeIcecap_v6bp3.m
%   v1: strip files
%   v2: scene files

if 1 % use mat2.mat

constant

close all
 
        params.I = [1 0 0 0 0 0 0];
        params.C = [50 0.001 0.001 0.05 0.0001];
%         params.G = [3000 20];
        params.M = '3D';
        params.V = [10 20 10];
        params.G = [9000 20]; %Adjust the max height parameter for the 2012 Kamchatka Volcano


res='2';
if 0
% demdir=dir(['/data2/ArcticDEM/region_',regionnum,'*']);
macdir='/Users/chunlidai/surge/';
macdir='/Users/chunli/surge/';
macdir=[];
addpath([macdir,'//home/chunli/scripts/MJ_coreg'])
addpath([macdir,'//home/chunli/scripts'])
end

yr=365.25;
neq=1;
dsr=0.05;%2m to 40m; 0.2;%0.04; %200m ; %0.2; % 8m to 40m
nsr=1./dsr;
resr=2;%str2double(res)/dsr;
resrc=40.; %for coregisteration
% demdir=[macdir,'/data2/ArcticDEM/region_08_canada_baffin/tif_results/8m/'];
coregflag=3;%1 parameter (vertical), 3 parameter (Ian's); 7 parameter (MJ's)
tifflag=1; %1 tif file; 2 strip file

if 0
% Tanana River 64.674662, -148.333690
loneq=-148.333690;lateq=64.674662; % 
% filename='boundaries_reg08.dat'; %
% filename='boundaries_reg3134_2m.dat';
filename='boundaries_reg3134_strip2m.dat';
filename='boundaries_reg34_2m.dat';
[xeq,yeq]=polarstereo_fwd(lateq,loneq,6378388.0,0.0819919,70,-45);

[gaugex,gaugey]=polarstereo_fwd(64+47/60+34/3600,  - (147+50/60+20/3600),[],[],70,-45);
xeq=gaugex;yeq=gaugey;
filename='boundaries_reg34_2md.dat';
[lateq,loneq]=polarstereo_inv(xeq,yeq,6378388.0,0.0819919,70,-45);
exb=0;

% xeq=-2.72e6;yeq= 1.455e6; %Yukon RIver
xeq=-2.708e6;yeq=6.28e5; %Tanana River
[lateq,loneq]=polarstereo_inv(xeq,yeq,6378388.0,0.0819919,70,-45);

filename='boundaries_strip2m.dat';
loneq=0;lateq=90; % 
[xeq,yeq]=polarstereo_fwd(lateq,loneq,6378388.0,0.0819919,70,-45);

%saglat=69+0+57/3600;saglon=-(148+49/60+04/3600);
filename='boundaries_reg34_2mc.dat';
loneq=-(148+49/60+04/3600);lateq=69+0+57/3600; % 
[xeq,yeq]=polarstereo_fwd(lateq,loneq,[],[],70,-45);
end

%exb=6e3;
%exb=20e3;
exb=40e3; %sag reach 34 km by 6 km; work size 10km by 10km (~86 images) is managable for memory;
exb=10e3;
loneq=source(1);lateq=source(2);
[xeq,yeq]=polarstereo_fwd(lateq,loneq,[], [],70,-45);
formatSpec = '%6.1f';
Co=[0 0];

rang0=[xeq-exb xeq+exb yeq-exb yeq+exb ];
% rang0=[-3467 -3455 110 124 ]*1e3;
% rang0=[-3288 -3281 352 358]*1e3;
x0=[rang0(1) rang0(2) rang0(2) rang0(1) rang0(1) ];y0=[rang0(4) rang0(4) rang0(3) rang0(3) rang0(4) ];
ranget=round(rang0/resr)*resr;rang0=ranget;
tx=ranget(1):resr:ranget(2);ty=ranget(4):-resr:ranget(3);
xout=tx;yout=ty;

x=[range(:,1) range(:,2) range(:,2) range(:,1) range(:,1) ];y=[range(:,4) range(:,4) range(:,3) range(:,3) range(:,4) ];
% id=find(range(:,1)>xeq-exb & range(:,2)<xeq+exb & range(:,3)>yeq-exb & range(:,4)<yeq+exb);
idregion=find(range(:,2)>xeq-exb & range(:,1)<xeq+exb & range(:,4)>yeq-exb & range(:,3)<yeq+exb);
% id=1:length(range(:,1));

x2=[range2(:,1) range2(:,2) range2(:,2) range2(:,1) range2(:,1) ];y2=[range2(:,4) range2(:,4) range2(:,3) range2(:,3) range2(:,4) ];
idregion2=find(range2(:,2)>xeq-exb & range2(:,1)<xeq+exb & range2(:,4)>yeq-exb & range2(:,3)<yeq+exb);

ck1=clock;

% length(id)

% idgage=[];
% plot orthoimages
for j=1:0 %length(id)
    i=id(j);
   	demdir=fdir{i};
    
infile= strrep([demdir,'/',f{i}],'meta.txt','ortho.tif');
if ~exist(infile);continue;end
clear data
data=readGeotiff(infile);
if 0
        ranget=[[min(data.x) max(data.x) min(data.y) max(data.y)]/resr];
        ranget=[ceil(ranget(1)) floor(ranget(2)) ceil(ranget(3)) floor(ranget(4))]*resr;
        tx=ranget(1):resr:ranget(2);ty=ranget(4):-resr:ranget(3);
res=data.info.map_info.dx;
nsr=resr/res;
        idrxs=find(abs(data.x-ranget(1))<1e-3);idrxe=find(abs(data.x-ranget(2))<1e-3);
        idrys=find(abs(data.y-ranget(4))<1e-3);idrye=find(abs(data.y-ranget(3))<1e-3);
        idrx=idrxs:nsr:idrxe;
        idry=idrys:nsr:idrye;
        dd=[idrx(2:end)-idrx(1:end-1),idry(2:end)-idry(1:end-1)];
        if isempty(dd)||any(abs(dd-nsr)>1e-3);warning('Resized grid is not constant spacing.');end
        tz=data.z(idry,idrx); %0.09s
        datar= struct();
        datar.x=tx;datar.y=ty;  datar.z=tz;
data=datar;
end
figure;set(gcf,'Color','white');set(gca,'FontSize', 12);set(gcf, 'PaperPosition', [0.25 2.5 6 4]);hold all;
imagesc(data.x*1e-3,data.y*1e-3,data.z);colormap gray; colorbar; shading interp; axis tight; view(0,-90)
title(['Orthoimage ',f{i}(6:13)]);
plot(xeq*1e-3,yeq*1e-3,'r*','Markersize',12);
hold on;plot(xriv(:,1),yriv(:),'b-','linewidth',4);
axis equal
view(0,90)
saveas(gcf,['/data2/saturationTest/TananaGage/',f{i}(6:13),'i',num2str(i),'Ortho'],'fig')

infile= strrep([demdir,'/',f{i}],'meta.txt','dem.tif');
data=readGeotiff(infile);
data.z(data.z== -9999) = NaN;
hills=hillshade(double(data.z),data.x,data.y,'plotit');
set(gcf,'Color','white');set(gca,'FontSize', 12);set(gcf, 'PaperPosition', [0.25 2.5 4 3]);
colorbar;%caxis([0 250])
hold on
plot(xeq,yeq,'r*','Markersize',12)
hold on;plot(xriv,yriv,'b-','linewidth',4)
title([f{i}(6:13)]);
axis equal
saveas(gcf,['/data2/saturationTest/TananaGage/',f{i}(6:13),'i',num2str(i),'DEM2m'],'fig')

z2n = interp2(data.x' ,data.y,data.z , xeq,yeq,'*linear');
    if ~isnan(z2n)
%     idgage=[idgage;i];
    end

end
% save idgage.mat idgage

%get water mask from a priori shapefile
if flagsect==1 %section of river
sectionname='sagExtent.shp';
S = shaperead(sectionname);
cnt=length(S); %figure;mapshow(S);

xw=rang0(1):resrc:rang0(2);yw=rang0(4):-resrc:rang0(3); % a priori water mask
nwx=length(xw);
%  %poly2mask, fast 0.5 sec.
smg=false(nwx,nwx);%water mask from a priori coastline shapefiles
for j=1:cnt
    [sx,sy]=polarstereo_fwd(S(j).Y,S(j).X,[], [],70,-45);
    dx=resrc;
    idx=round((sx-rang0(1))/dx)+1;idy=round((sy-rang0(4))/(-dx))+1;
    %Oct 10, 2018:fix bug 14; separate polygons using the NaN;
    M=[idx(:),idy(:)];
    idx = any(isnan(M),2);
    idy = 1+cumsum(idx);
    idz = 1:size(M,1);
    C = accumarray(idy(~idx),idz(~idx),[],@(r){M(r,:)});
    for k=1:length(C)
    idx=C{k}(:,1);idy=C{k}(:,2);   
    sm=poly2mask(idx,idy,nwx,nwx); % fast, apply to each polygon one by one.
    smg=smg|sm;
    end
end
wm=[];wm.x=xw;wm.y=yw;wm.z=smg;clear smg; 
if(sum(wm.z(:))==0||sum(wm.z(:))==nwx*nwx);fprintf('This tile contain no coastline (all land or all ocean).');return;end % if all land or all water, i.e. no coastline.
Mcb=wm.z; % %1 water, 0 land
end % if flagsect

fprintf ('\n Step 1: getting the real Polygon boundary for all files over the output zone.')

enl=50; %meter %100m per node
x0si=[xeq-enl xeq+enl xeq+enl xeq-enl xeq-enl];
y0si=[yeq+enl yeq+enl yeq-enl yeq-enl yeq+enl];


%use stereo images
%load all reg.txt and meta.txt files in the coverage
XYb2=cell(size(idregion2));dzxy2=cell(size(idregion2));count=0;
flagcb2=zeros(size(idregion2));
flagcb2s=zeros(size(idregion2));
for j=1:length(idregion2)
    i=idregion2(j);
 % get the Polygon boundary for actual data
        demdir=fdir2{i};
        infile= [demdir,'/',f2{i}];
        %[XYbi,rangei]=imagebd(infile);
        XYbi=XYbg2{i};
        Xb=XYbi(:,1);Yb=XYbi(:,2);
        XYb2{j}=XYbi;
        
	if flagsect==1 %1 work on a section of river; 0 work on a gage
        % check whether this polygon intersect with the river section
        idx=round((Xb-wm.x(1))/resrc)+1;
        idy=round((Yb-wm.y(1))/(-resrc))+1;
        Mb = poly2mask(idx,idy, nwx,nwx); % build polygon mask       
        overl=Mb&Mcb;
        if(sum(sum(overl))~=0);flagcb2s(j)=0;
	    fprintf(['\n section ifile:',infile])
        else
	    flagcb2s(j)=1;
	    fprintf(['\n nonsection ifile:',infile])
        end
	end

        in = inpolygon(x0si,y0si,Xb,Yb); %whether the gage is inside the polygon
        if any(in==0) %any subtile corners not in the boundary
            flagcb2(j)=1;
            %print all for sag
            if abs(lateq - 69.01583)<1e-4 % if sag station      32    69.015833  -148.817778    
            fprintf(['\n no gage ifile:',infile])
            end
        else
            fprintf(['\n ifile:',infile])
        end
end

% only keep the strips that cover the gage.
if flagsect==1 %section of river
idd=flagcb2s==1; %not on section
else
idd=flagcb2==1; %not have gage.
end
idregion2(idd)=[];XYb2(idd)=[];dzxy2(idd)=[];
%XYb2 dzxy2 not really used.

idregion2all=idregion2;
%only keep summer scenes.
id=idregion2;mon=zeros(length(id),1);
for j=1:length(id)
	ymd=f2{id(j)}(6:13);i=id(j);
        mon(j)=str2num(ymd(5:6));
end
idd=~(mon>=mons&mon<=mone); 
%idd=[]; %keep all strips; won't work since only strips with water masks are used.
idregion2(idd)=[];XYb2(idd)=[];dzxy2(idd)=[];

%When count: get rid of the strips have the same date and same sensor
count=zeros(2,1);text1={'all','summer'};
for k=1:2
if k==1
id=idregion2all;str=cell(length(id),1);
elseif k==2
id=idregion2;str=cell(length(id),1);
end
for j=1:length(id)
            % get the Polygon boundary for actual data
            i=id(j);
            str{j}=f2{id(j)}(1:13);
end
[un idx_last idx] = unique(str(:));
id1=1:length(id);idd=id1(~ismember(id1,idx_last));

count(k)=length(id)-length(idd);
fprintf(['\n Number of images cover the gage: ',num2str(count(k)), ' ',text1{k}])
end

%if(isempty(idregion2));fprintf('No images along the gage.');return;end % 

% Co=count;return;
%End of use stereo images


if 1  %||flagmono==1 %1 use mono images only; 2 use stereor; Always use mono for water mask, since strip orthoiamge has problem for water masking.
%load all reg.txt and meta.txt files in the coverage
XYb=cell(size(idregion));dzxy=cell(size(idregion));count=0;
flagcb=zeros(size(idregion));
flagcbs=zeros(size(idregion));
for j=1:length(idregion)
    i=idregion(j);
 % get the Polygon boundary for actual data
        demdir=fdir{i};
        infile= [demdir,'/',f{i}];
        %[XYbi,rangei]=imagebd(infile);
        XYbi=XYbg{i};
        Xb=XYbi(:,1);Yb=XYbi(:,2);
        XYb{j}=XYbi;
        
	if flagsect==1 %1 work on a section of river; 0 work on a gage
        % check whether this polygon intersect with the river section
        idx=round((Xb-wm.x(1))/resrc)+1;
        idy=round((Yb-wm.y(1))/(-resrc))+1;
        Mb = poly2mask(idx,idy, nwx,nwx); % build polygon mask       
        overl=Mb&Mcb;
        if(sum(sum(overl))~=0);flagcbs(j)=0; %
	    fprintf(['\n section ifile:',infile])
        else
	    flagcbs(j)=1;
	    fprintf(['\n nonsection ifile:',infile])
        end
    end
    
        in = inpolygon(x0si,y0si,Xb,Yb); %whether the gage is inside the polygon
        if any(in==0) %any subtile corners not in the boundary
          flagcb(j)=1;
	  %print all for sag
	  if abs(lateq - 69.01583)<1e-4 % if sag station      32    69.015833  -148.817778    
	    fprintf(['\n no gage ifile:',infile])
	  end
	else
	    fprintf(['\n ifile:',infile])
	end
end

% only keep the strips that cover the gage.
if flagsect==1 %section of river
idd=flagcbs==1; %not on section
else
idd=flagcb==1; %not have gage.
end
idregion(idd)=[];XYb(idd)=[];dzxy(idd)=[];

idregionall=idregion;
%only keep summer scenes.
id=idregion;mon=zeros(length(id),1);
for j=1:length(id)
	ymd=f{id(j)}(6:13);i=id(j);
        mon(j)=str2num(ymd(5:6));
end
idd=~(mon>=mons&mon<=mone); 
%idregion=id(mon>=5&mon<=10); %summer season.
idregion(idd)=[];XYb(idd)=[];dzxy(idd)=[]; %Maybe use both summer and winter, NDWI adaptive thresholding may be able to delete winter scenes.

%When count: get rid of the strips have the same date and same sensor
count=zeros(2,1);text1={'all','summer'};
for k=1:2
if k==1
id=idregionall;str=cell(length(id),1);
elseif k==2
id=idregion;str=cell(length(id),1);
end
for j=1:length(id)
            % get the Polygon boundary for actual data
            i=id(j);
            str{j}=f{id(j)}(1:13);
end
[un idx_last idx] = unique(str(:));
id1=1:length(id);idd=id1(~ismember(id1,idx_last));

count(k)=length(id)-length(idd);
fprintf(['\n Number of images cover the gage: ',num2str(count(k)), ' ',text1{k}])
end

%if(isempty(idregion));fprintf('No images along the gage.');return;end % 

Co=count;
end

if 0 %use stereo iamges;
    idregion=idregion2;XYb=XYb2;dzxy=dzxy2;
    range=range2;x=x2;y=y2;XYbg=XYbg2;f=f2;
    fdir=fdir2;
end % if


%End of loading

odir=['./gage',num2str(igage)];
if ~exist(odir,'dir')
  mkdir(odir)
end

id=idregion;
%Plot the coverage
figure
set(gcf,'Color','white')
set(gca,'FontSize', 12);
set(gcf, 'PaperPosition', [0.25 2.5 6 4.5]); %Default [0.25 2.5 8.0 6.0];Prefer 6 by 4.5 or 4 by 3
hold all
% plot(x(id), y(id) ,'.-')
for i=1:length(id)
    plot(x(id(i),:), y(id(i),:),'b>-','linewidth',4)
%     hold off
    title(['j=',num2str(i),';i=',num2str(id(i)),';',f{id(i)}(6:13)]);
t(i)=str2num(f{id(i)}(10:11));
if t(i)>= 5 && t(i)<= 10 % Tanana running season
    plot(x(id(i),:), y(id(i),:),'g>-','linewidth',4);
end
end
hold on
plot(xeq,yeq,'r*','Markersize',12)
saveas(gcf,'rangelist','fig')

id=idregion2;
figure
set(gcf,'Color','white')
set(gca,'FontSize', 12);
set(gcf, 'PaperPosition', [0.25 2.5 6 4.5]); %Default [0.25 2.5 8.0 6.0];Prefer 6 by 4.5 or 4 by 3
hold all
for i=1:length(id)
    plot(x2(id(i),:), y2(id(i),:),'b>-','linewidth',4)
    title(['j=',num2str(i),';i=',num2str(id(i)),';',f2{id(i)}(6:13)]);
t(i)=str2num(f2{id(i)}(10:11));
if t(i)>= 5 && t(i)<= 10 % Tanana running season
    plot(x2(id(i),:), y2(id(i),:),'g>-','linewidth',4);
end
end
hold on
plot(xeq,yeq,'r*','Markersize',12)
saveas(gcf,'rangeliststrip','fig')

id=idregion;
[lat,lon]=polarstereo_inv(x,y,[],[],70,-45);
% lon(lon>=0)=lon(lon>=0)-360;
figure
set(gcf,'Color','white')
set(gca,'FontSize', 12);
set(gcf, 'PaperPosition', [0.25 2.5 6 4.5]); %Default [0.25 2.5 8.0 6.0];Prefer 6 by 4.5 or 4 by 3
hold all
% plot(lon', lat' ,'-')
for i=1:length(id)
    plot(lon(id(i),:), lat(id(i),:),'b>-','linewidth',4)
%     hold off
end
plot(loneq, lateq ,'r*','Markersize',12)

fprintf ('\n Step 1.2: Searching for overlappings of Polygons at regular big grids.')
tic
%sub-zones for searching overlappings
dx=800;%400;%3600;%7e3;%%17e3; %3600; %80; %3600; 
dxov=dx/4;%200;
rang0t=round(rang0/dxov)*dxov;
rangx=rang0t(1):dx:rang0t(2);nx=length(rangx)-1;
rangy=rang0t(3):dx:rang0t(4);ny=length(rangy)-1;
ns=nx*ny;rang0s=zeros(ns,4);
novlp=zeros(ny,nx);baovlp=zeros(ny,nx);idg=cell(ns,1);%idg{ns}=[];
enl=0;%;0.3;
for ix=1:nx
    for iy=1:ny
        ixy=iy+(ix-1)*ny;
        rang0s(ixy,:)=[rangx(ix) rangx(ix+1) rangy(iy) rangy(iy+1) ];
        id=find(range(:,1)<=rangx(ix)-enl*dx & range(:,2)>=rangx(ix+1)+enl*dx & range(:,3)<=rangy(iy)-enl*dx  & range(:,4)>=rangy(iy+1)+enl*dx);
%         novlp(iy,ix)=length(id);
        % Refinement of the overlapping count, checking the actual data
        % coverage overlapping with the subzone.
        x0si=[rang0s(ixy,1)-enl*dx rang0s(ixy,2)+enl*dx rang0s(ixy,2)+enl*dx rang0s(ixy,1)-enl*dx rang0s(ixy,1)-enl*dx ];
        y0si=[rang0s(ixy,4)+enl*dx rang0s(ixy,4)+enl*dx rang0s(ixy,3)-enl*dx rang0s(ixy,3)-enl*dx rang0s(ixy,4)+enl*dx ];
        idd=[]; str=cell(length(id),1);
        for j=1:length(id)
            % get the Polygon boundary for actual data
            i=id(j); 
            str{j}=f{id(j)}(1:13);

            M=idregion==i;nt=sum(M); % only work on the strips that are selected (cover coastline and successfully coregistered).
            if nt==0 ;idd=[idd;j];continue;end
            Xb=XYb{M}(:,1);
            Yb=XYb{M}(:,2);

            % new method
            in = inpolygon(x0si,y0si,Xb,Yb); %whether subtiles are inside the polygon
            if any(in==0) %any subtile corners not in the boundary
               idd=[idd;j];
            end
        end % j=1:length(id)
        id(idd)=[];str(idd)=[];

        %get rid of the strips have the same date and same sensor
        [un idx_last idx] = unique(str(:));
        id1=1:length(id);idd=id1(~ismember(id1,idx_last));
        id(idd)=[];str(idd)=[];

        novlp(iy,ix)=length(id);
        idg{ixy}=sort(id);   %finding the DEMs at each zones.    

    end
end
display(['Counting overlapping...']);
toc

% [X,Y]=meshgrid(rangx(1:end-1),rangy(1:end-1));
[X,Y]=meshgrid((rangx(1:end-1)+rangx(2:end))/2,(rangy(1:end-1)+rangy(2:end))/2);
figure;set(gcf,'Color','white');set(gca,'FontSize', 18);set(gcf, 'PaperPosition', [0.25 2.5 4 3]); 
surf(X*1e-3,Y*1e-3,novlp);colormap jet;colorbar; shading interp; axis equal; view(0,90)  
set(gca,'FontSize', 18);
hold on
plot3(xeq*1e-3, yeq*1e-3,1e2 ,'r*','Markersize',12)
xlabel('{\itx} (km)')
ylabel('y (km)')
saveas(gcf,'OverlappingCount_Alaska','fig')

[LAT,LON]=polarstereo_inv(X,Y,[],[],70,-45);
figure;set(gcf,'Color','white');set(gca,'FontSize', 18);set(gcf, 'PaperPosition', [0.25 2.5 4 3]); 
surf(LON,LAT,novlp);colormap jet;colorbar; shading interp; axis equal; view(0,90)  
set(gca,'FontSize', 18);
hold on
plot3(loneq, lateq,1e2 ,'r*','Markersize',12)
saveas(gcf,'OverlappingCount','fig')

% Co=count;return;
% exit

close all
% Initialize
nsuby=length(yout);nsubx=length(xout);
prob=255*ones(nsuby,nsubx,'uint8'); %probability of the final water mask at each pixel. 255 no data
%oflagc=255*ones(nsuby,nsubx,'int32');
Medgsib=false(nsuby,nsubx); %collect of all images edge lines.

% Preallocate 
isv=idregion;
datarsv(length(isv))=struct('x',[],'y',[],'z',[],'n',[]);
datamtrsv(length(isv))=struct('x',[],'y',[],'z',[]);

%Get the road data for this gage region.
%groad=roadgage(rang0,odir);

% %%%% coregister DEM strips to a reference DEM tile.
[dX4Sg2,idregion2,data0r]=coreg(rang0,idregion2,XYbg2,f2,fdir2);
save RefDEM.mat data0r -v7.3
%coregister mono images to a reference DEM tile.
%dX4Sg=zeros(length(idregion),3); 
%collect the reference mono images that are components of strip files in idregion2.
%And to coregister all mono images
id=idregion; fis=f(id);fdiris=fdir(id);XYbis=XYbg(id);
id=idregion2;f2is=f2(id);fdir2is=fdir2(id);XYb2is=XYbg2(id);
if ~exist('dX4Sg.mat','file')
[dX4Sg]=pxpymono(dX4Sg2,f2is,fdir2is,XYb2is,fis,fdiris,XYbis); %12 hours for 86 images;14 hours for 86 images
save dX4Sg.mat dX4Sg
else
load dX4Sg.mat %hi
end

% %%%% Get initial water mask step 1: run the multispectral images only %%%
% %%%% Get refined water mask step 2; using the water mask in step 1 as a priori mask; using multispectral images only %%%
% %%%% Get refined water mask step 3; using the water mask in step 2 as a priori mask; using both multispectral and pan images 
id=idregion;
% get coastline/watermask for each file, and save it.
t=zeros(length(id),1);
%load mat1.mat
%load dX4Sg.mat %hi
for step=1:3   
    fprintf(['Step ',num2str(step), '...\n'])
    Msuma=zeros(nsuby,nsubx,'int32');%maximum repeat at a pixel; %repeats of non void data.
    Msum=zeros(nsuby,nsubx,'int32');%maximum repeat of water pixel.

    if step==1
        wm.x=[];wm.y=[];wm.z=[];
        wm1.x=[];wm1.y=[];wm1.z=[];
    elseif step>=2
	[wm]=prepwm(xout,yout,jump,prob,step);
        if step<3 ||flagbuf==0
           wm= rmfield(wm,'buf'); 
        end
	wm1=wm;%To be compatible with older versions.
    end
    
%   [~,idsort]=sort(t);id=id(idsort); %sort the id based on tim
    idd=[];
	for j=1:length(id)
        display(['Working on strip ',num2str(j), ': ',f{id(j)}])
        ymd=f{id(j)}(6:13);i=id(j); 
        t(j)=datenum(ymd,'yyyymmdd');
        iisv=find(isv==i);%iisv ==j
        demdir=fdir{i};
        satname=f{i}(1:4);
        
        % find the overlapping range of data range and box.
        rang0sov=range(i,:);
        rangeov=[max(rang0sov(1),rang0(1)),min(rang0sov(2),rang0(2)), max(rang0sov(3),rang0(3)),min(rang0sov(4),rang0(4))];
        
        infile= [demdir,'/',f{i}];
        [demdir,name,ext] =fileparts([(infile)]);
        name1=[name(end-3:end),ext];
        % To check whether the data is mono or scene files
	flagfmt=0;
        if strcmp(name1,'meta.txt')
        flagfmt=3; %stereo files
        elseif strcmp(ext,'.xml')
        flagfmt=1; %xml files mono
        end
	
	%check whether the data is multispectral image
	flagfmt2=0; %panchromatic band
	r1=strfind(name,'M1BS');
	if ~isempty(r1); flagfmt2=1;end
        
    if isempty(datarsv(iisv).x) % non exist ; step 1,
        ndwisv=struct('x',[],'y',[],'z',[]);%Initialize
        data=struct('x',[],'y',[],'z',[],'n',[]);%Initialize
              
        flagproc=1;
        if flagfmt2==1 %multispectral image
%             infile= strrep([demdir,'/',f{i}],'.xml','.xml');
            clear data        
            data=multispecmono(infile,wm,ndwisv,rangeov); %given strip meta file, finding all image 1 multispectral imageries, orthrectiying, get water mask.
   
        elseif step==3  %stereo orthorectified panchromatic image; or WV01 mono panchromatic image
%             infile=  strrep([demdir,'/',f{i}],'meta.txt','dem.tif');
            data=maskentropy(infile,wm1,ndwisv,rangeov);
        else % nothing processed.
            flagproc=0;
            idd=[idd;j]; continue
        end
        display(['Getting water mask of this strip ',num2str(j)]);
        
        if isempty(data.z)&&flagproc==1
	%Modification Sept. 25, 2018: save the empty results to avoid repeat calculation of NDWI.
	if 0 %NDWI is already saved. Better to reprocess the water mask with new threshold.
        datar= struct();
        datar.x=1;datar.y=[];  datar.z=[]; datar.n=[];
        datarsv(iisv)=datar;
	end

            idd=[idd;j]; continue
        end  
    if 0
        figure;
        imagesc(data.x,data.y,data.z)
        hold all;plot(x0st, y0st,'r-','linewidth',6)
        plot(x0, y0,'g-','linewidth',4)
        plot(xeq, yeq ,'r*','Markersize',12)
        title(['isel=',num2str(isel),f{i}])
        plot(Cl(1,:),Cl(2,:),'.-','Color','m','linewidth',2) %r+-
        axis(rang0)
        colormap gray
        colorbar
        ofile=[ymd,f{i}(1:4),'isel',num2str(isel),'Ortho']; %20130414WV01isel2Ortho
        saveas(gcf,ofile,'fig')
    end
        
        %not neccesary
        ranget=[[min(data.x) max(data.x) min(data.y) max(data.y)]/resr];
        ranget=[ceil(ranget(1)) floor(ranget(2)) ceil(ranget(3)) floor(ranget(4))]*resr;
        tx=ranget(1):resr:ranget(2);ty=ranget(4):-resr:ranget(3);
        tz = interp2(data.x,data.y,data.z,tx,ty','*nearest',-1);
        tn = interp2(data.x,data.y,data.n,tx,ty','*nearest',nan);%ndwi
%         tz = interp2(data.x,data.y,data.z,tx,ty','*nearest',-1); %NDWI or brightness.
%       tz(isnan(tz))=-1;tz=int8(tz); %DOUBLE TO INT8
        datar= struct();
        datar.x=tx;datar.y=ty;datar.z=tz; datar.n=tn;%datar.coast=data.coast;

        datarsv(iisv)=datar; %possible crash with parallel; original coordinates

    dzxyt=dX4Sg(idregion==i,:); % possible wrong match
    dX4S=dzxyt(2:3); 

	%collect all image edge lines
	tzedg=tz==-1; %datar.z==-1;
	tzd1=imdilate(tzedg, ones(3));
	tze1=imerode(tzedg,ones(3));
	tzl=tzd1&(~tze1); %edge lines
	%[X,Y]=meshgrid(datar.x,datar.y);
	tz = interp2(datar.x-dX4S(1),datar.y-dX4S(2),tzl,xout,yout','*nearest',0); %check time
	Medgsib=Medgsib|tz; %edge 1; non edge 0
      
	else %load water mask data for collecting demg; step 2
        datar=datarsv(iisv);

        %Modification Sept. 25, 2018:skip empty results.
        if isempty(datar.z)
                idd=[idd;j];
                continue
        end
        
        if step>=2
            ndwisv.x=datar.x;ndwisv.y=datar.y; ndwisv.z=datar.n;% The .x, .y coordinates are already coregistered to a reference in step 1.
            
%       if flagmono==1&&~strcmp(satname,'WV01')&&flagfmt==1 %mono multispectral image
        if flagfmt2==1 %multispectral image
            data=multispecmono(infile,wm,ndwisv,rangeov); %given strip meta file, finding all image 1 multispectral imageries, orthrectiying, get water mask.
        else %stereo orthorectified panchromatic image; or WV01 mono panchromatic image
            data=maskentropy(infile,wm1,ndwisv,rangeov);
        end
        datar.x=data.x;datar.y=data.y;datar.z=data.z; datar.n=data.n;
        datarsv(iisv)=datar; %update the water mask for step 2;

        if isempty(data.z) %Using a priori information, the water mask can be void if quality bad.
	%Modification Sept. 25, 2018: save the empty results to avoid repeat calculation of NDWI.
	if 0 %NDWI is already saved. Better to reprocess the water mask with new threshold.
        datar= struct();
        datar.x=1;datar.y=[];  datar.z=[]; datar.n=[];
        datarsv(iisv)=datar;
	end
            idd=[idd;j];
            continue
        end

        else %step 1
             % do nothing
        end

    end %if
    
    %Apply the offset when mosaicing the masks.
        %coregister to a reference (e.g. ICESat). reg.txt file or feature tracking.
%         dzxyt=dzxyd(idregion==i,:);%dzxy{idregion==i};%lots of them missing along coast.
    dzxyt=dX4Sg(idregion==i,:); % possible wrong match
    dX4S=dzxyt(2:3); 

    %make coordinates exactly on (resr*n) xout yout grids.
    datar.x=datar.x-dX4S(1);datar.y=datar.y-dX4S(2);
    
    %assign DEM to the box.
    tz = interp2(datar.x,datar.y,datar.z,xout,yout','*nearest',-1);
    M=tz~=-1; % data pixel;
    Msuma(M)=Msuma(M)+1; %repeats of non void data.
    M=tz==1; %water pixel;%1, water, 0 non water, -1 void data.
    Msum(M)=Msum(M)+1;
    
	
    if 0 %flagplot==1&&step>=2
    [X,Y]=meshgrid(data.x,data.y);
    M=data.coast;
    figure;
    hold all
    set(gcf,'Color','white');set(gca,'FontSize', 18);set(gcf, 'PaperPosition', [0 0 6 4]);set(gcf, 'PaperSize', [ 6 4]);
    imagesc(datar.x*1e-3,datar.y*1e-3,double(datar.n))
    colormap jet;colorbar
    hold on;plot(X(M)*1e-3,Y(M)*1e-3,'ko')
    if flagfmt2==1 %multispectral image
    title(['NDWI ',name(6:13)]);
    elseif flagfmt2==0
    title(['Pan-band Brightness ',name(6:13)]);
    end
    axis equal
    saveas(gcf,[name],'fig')
    end
	close all
    
	%Oct 9, 2018
	%in case of two many repeats, use only the novmax measurements
% 	if j-length(idd) >= novmax;idd=[idd(:);[j+1:length(id)]'];  break;end

	end %for j
% % end of loading
% 	id(idd)=[];demg(:,:,idd)=[];t(idd)=[];
%save work1.mat demg
iprobthre=probthre;
if step==2 % Let the threshold be high, so it can be used to get Pan-brightness/NDWI over water area.
   iprobthre=90; %try 80; considearing clouds make part of river having less water probability.
		 % The first run is not accurate since the ithreshold is not adaptive, and bad images are not filtered.
		 %Second run should be accurate.
end

prob=uint8(double(Msum)./double(Msuma)*100);prob(Msuma==0)=255;%
prob(Msuma==1)=uint8(0.1*100);        %0.1, only one image;

jump=-1*ones(nsuby,nsubx,'int8'); %%-1,non value;1 water; 0 non-water
jump(prob>=iprobthre&prob~=255)=1;
jump(prob<iprobthre)=0;

%save wprob.mat xout yout prob jump Msum Msuma -v7.3
ofile=[odir,'/wprob',num2str(step),'.mat']; %step2 jump is using 90% threshold.
save(ofile,'xout','yout','prob','jump','Msum','Msuma','-v7.3')
whos %check memory usage

jump90=jump; %save the mask with 90% probability

if step ==2 %Output the final water mask, with 50% probability
save mat1.mat -v7.3

iprobthre=probthre;
jump=-1*ones(nsuby,nsubx,'int8'); %%-1,non value;1 water; 0 non-water
jump(prob>=iprobthre&prob~=255)=1;
jump(prob<iprobthre)=0;
jumpsv=jump;

%Plot water mask of 50% probability.
[X,Y]=meshgrid(xout,yout);
xmin=min(X(jump==1));xmax=max(X(jump==1));
ymin=min(Y(jump==1));ymax=max(Y(jump==1));
idx=xout>=xmin&xout<=xmax;idy=yout>=ymin&yout<=ymax;
%save files

ofile1=[odir,'/coast_v1.0.shp'];
ofile2=[odir,'/prob_v1.0.tif'];
ofile3=[odir,'/nov_v1.0.tif']; %number of overlapping coverages
projstr='polar stereo north';
OutName=ofile2;
xouto=xout(idx);youto=yout(idy);probo=prob(idy,idx);novlpf=Msuma(idy,idx);
if 0
writeGeotiff(OutName,xouto,youto,probo,1,255,projstr)
writeGeotiff(ofile3,xouto,youto,uint8(novlpf),1,255,projstr)
end

if 0 % figure too big; 
[LATa,LONa]=polarstereo_inv(X(idy,idx),Y(idy,idx),[],[],70,-45);
H=figure;
set(gcf,'Color','white')
set(gca,'FontSize', 12);
set(gcf, 'PaperPosition', [0.25 2.5 6 2.5]);
hold all;
surf(LONa,LATa,jump(idy,idx)); % too slow
shading interp;
colorbar;colormap jet;view(0,90)
plot(loneq, lateq ,'rs','Markersize',12)
view(0,90)
box on
hl=xlabel('Longitude ($^{\circ}$)');
set(hl, 'Interpreter', 'latex');
hl=ylabel('Latitude ($^{\circ}$)');
set(hl, 'Interpreter', 'latex');
ofile=[odir,'/watermask'];
%saveas(gcf,ofile,'fig')
%savefig(H,ofile,'compact')
try 
saveas(gcf,ofile,'fig')
catch e
     fprintf('There was an error! The message was:\n%s',e.message);
    try 
	savefig(H,ofile,'compact')
    catch e
	fprintf('There was an error! The message was:\n%s',e.message);
	try 
	   print('-dpdf','-r400',ofile)
	end
    end
end
close all
end %if 0
clear X Y %release memory
end %if step==2
jump=jump90; %use the mask with 90% probability

end % if step
fprintf(['Collecting the final water mask ']);

save mat2.mat -v7.3

end %if 0
%load('mat2.mat')

% % %
% Get the profiles of images at gage.
fprintf ('\n Step 2: Retrieve river shorelines location and elevation.')
id=idregion;
%Msv=datarsv; %
%dX4Sgis=dX4Sg;
fis=f(id);
fdiris=fdir(id);
XYbis=XYbg(id);
id=idregion2;
f2is=f2(id);
fdir2is=fdir2(id);
XYb2is=XYbg2(id);
%dX4Sg2is=dX4Sg2;

[Co]=riverprofsub(odir,datarsv,XYbis,fis,fdiris,dX4Sg,XYb2is,f2is,fdir2is,dX4Sg2);
%[Co]=riverprofsub(odir,wm,datarsv,XYbis,fis,fdiris,dX4Sg,XYb2is,f2is,fdir2is,dX4Sg2);

%Use the jump data in step2
if exist('jumpsv','var')
jump=jumpsv; %use the chosen mask for output and centerline
else 
  try
  ofile=[odir,'/wprob2.mat']; %caution: water mask with probability <= 90% may create too narrow rivers (even broken rivers).
  load(ofile)

  %get the 50% mask.
  iprobthre=probthre;
  [nsuby,nsubx]=size(prob);
  jump=-1*ones(nsuby,nsubx,'int8'); %%-1,non value;1 water; 0 non-water
  jump(prob>=iprobthre&prob~=255)=1;
  jump(prob<iprobthre)=0;

  end
end

%Get river centerline from water mask
if exist('clsv2.mat','file')
  load clsv2.mat %use the saved one %hi
else
%40m resolution -> river can be too narrow, so a 40 m resoltuion may break the rivers into sections.
%2m resolution cropped
res1=resr;
%data.x=xout(idx);data.y=yout(idy);tz=jump(idy,idx); %idx idy can be lost
data.x=xout;data.y=yout;data.z=int8(jump==1);
data=cropmatrix(data,data.z);
tz=data.z;

tz(tz==-1)=0;
Modj= bwareaopen(tz, round(lakearea/res1/res1)); %remove small clusters
Modfil = bwareaopen(~Modj, round(cloudarea/res1/res1)); %fill small areas: 1e4*4m^2
Modfil=~Modfil;
data.z=Modfil;
data=mask2river(data);data.z(data.z==-1)=0; %remove lakes, no fill

npt=sum(sum(data.z==1));
if npt>0
	try
	[c]=mask2centerline(data);
	save cl1.mat c
	catch e
                fprintf('There was an error! The message was:\n%s',e.message);
		fprintf('\n')
        end

	%adjust the centerline to go uphill for ProcessTananaFairbanks.m.
	width=500; %m 
	pt=[c.X(1) c.Y(1);c.X(end) c.Y(end)];
        [clx,cly]=polarstereo_fwd(pt(:,2),pt(:,1),[], [],70,-45); %lateq,loneq
	%use data0r to get the height
	pth = interp2(data0r.x',data0r.y,data0r.z, clx,cly,'*nearest');
	if pth(2)<pth(1) %downhill
        fprintf(['\n Flip centerline to go uphill. Heights of two ends:',num2str(pth')]);
	c.X=flip(c.X);
	c.Y=flip(c.Y);
	end
else
   fprintf(['No water pixels in the water mask. ']); 
end %if npt>0
end %if exist 

% % % Filtering of river profiles 
fprintf ('\n Step 3: Filtering of river profiles.')
%Notice: centerline cannot have NaNs; center line monotonic
M=isnan(c.X)|isnan(c.Y);
c.X(M)=[];c.Y(M)=[];
% Requires the input centerline to go uphill.
Co=ProcessTananaFairbanks(odir,c,lateq,loneq);
gageheights

whos
ck2=clock;
dt=ck2-ck1;tsec=dt(4)*3600+dt(5)*60+dt(6);
display(['Loading files and NDWI calculation take: ',num2str(tsec),' sec'])
clear datarsv data 

fprintf ('\n Step 3b: Width profile.')
Co=getwidthall(odir,c,lateq,loneq);
%co=widthgageprof(odir,c,lateq,loneq);


fprintf ('\n Step 4: Height time series analysis at the gage.')
%gageheights
[co]=gagewidth(odir);

fprintf ('\n Step 5: Discharge Time series analysis at gage.')
%stagedischarge

close all

return
end


