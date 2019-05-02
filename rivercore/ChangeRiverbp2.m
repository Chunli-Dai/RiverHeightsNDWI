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
resr=str2double(res)/dsr;
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
exb=40e3; %sag reach 34 km by 6 km
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
id=find(range(:,2)>xeq-exb & range(:,1)<xeq+exb & range(:,4)>yeq-exb & range(:,3)<yeq+exb);
% id=1:length(range(:,1));

x2=[range2(:,1) range2(:,2) range2(:,2) range2(:,1) range2(:,1) ];y2=[range2(:,4) range2(:,4) range2(:,3) range2(:,3) range2(:,4) ];
idregion2=find(range2(:,2)>xeq-exb & range2(:,1)<xeq+exb & range2(:,4)>yeq-exb & range2(:,3)<yeq+exb);

% length(id)
riv=load('FAirport2.gmt');
[xriv,yriv]=polarstereo_fwd(riv(:,2),riv(:,1),[],[],70,-45);

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
imagesc(data.x,data.y,data.z);colormap gray; colorbar; shading interp; axis tight; view(0,-90)
title(['Orthoimage ',f{i}(6:13)]);
plot(xeq,yeq,'r*','Markersize',12);
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
wm=[];wm.x=xw;wm.y=yw;wm.z=smg;clear smg; %1 land, 0 water
if(sum(wm.z(:))==0||sum(wm.z(:))==nwx*nwx);fprintf('This tile contain no coastline (all land or all ocean).');return;end % if all land or all water, i.e. no coastline.
Mcb=wm.z; % %1 water, 0 land
end % if flagsect

fprintf ('\n Step 1: Preparing the grouping of common overlapping pieces.')
fprintf ('\n Step 1.1: getting the real Polygon boundary for all files over the output zone.')
if flagmono==1 %1 use mono images only; 2 use stereor
%load all reg.txt and meta.txt files in the coverage
idregion=id;XYb=cell(size(idregion));dzxy=cell(size(idregion));count=0;
flagcb=zeros(size(idregion));
flagcbs=zeros(size(idregion));
enl=50; %meter %100m per node
x0si=[xeq-enl xeq+enl xeq+enl xeq-enl xeq-enl];
y0si=[yeq+enl yeq+enl yeq-enl yeq-enl yeq+enl];
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
        if(sum(sum(overl))~=0);flagcbs(j)=1;
	    fprintf(['\n not on section ifile:',infile])
        else
	    fprintf(['\n on section ifile:',infile])
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
idd=flagcb==1; %not have gage.
idregion(idd)=[];XYb(idd)=[];dzxy(idd)=[];

idregionall=idregion;
%only keep summer scenes.
id=idregion;mon=zeros(length(id),1);
for j=1:length(id)
	ymd=f{id(j)}(6:13);i=id(j);
        mon(j)=str2num(ymd(5:6));
end
idd=~(mon>=5&mon<=10); 
%idregion=id(mon>=5&mon<=10); %summer season.
idregion(idd)=[];XYb(idd)=[];dzxy(idd)=[];

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

if(isempty(idregion));fprintf('No images along the gage.');return;end % 

Co=count;
end % if

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
        if(sum(sum(overl))~=0);flagcb2s(j)=1;
	    fprintf(['\n not on section ifile:',infile])
        else
	    fprintf(['\n on section ifile:',infile])
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
idd=flagcb2==1; %not have gage.
idregion2(idd)=[];XYb2(idd)=[];dzxy2(idd)=[];

idregion2all=idregion2;
%only keep summer scenes.
id=idregion2;mon=zeros(length(id),1);
for j=1:length(id)
	ymd=f{id(j)}(6:13);i=id(j);
        mon(j)=str2num(ymd(5:6));
end
idd=~(mon>=5&mon<=10); 
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
            str{j}=f{id(j)}(1:13);
end
[un idx_last idx] = unique(str(:));
id1=1:length(id);idd=id1(~ismember(id1,idx_last));

count(k)=length(id)-length(idd);
fprintf(['\n Number of images cover the gage: ',num2str(count(k)), ' ',text1{k}])
end

if(isempty(idregion2));fprintf('No images along the gage.');return;end % 

Co=count;return;


%End of loading

odir=['gage',num2str(igage)];
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

[lat,lon]=polarstereo_inv(x,y,6378388.0,0.0819919,70,-45);
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
%try inpolygon %hi
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

[LAT,LON]=polarstereo_inv(X,Y,6378388.0,0.0819919,70,-45);
figure;set(gcf,'Color','white');set(gca,'FontSize', 18);set(gcf, 'PaperPosition', [0.25 2.5 4 3]); 
surf(LON,LAT,novlp);colormap jet;colorbar; shading interp; axis equal; view(0,90)  
set(gca,'FontSize', 18);
hold on
plot3(loneq, lateq,1e2 ,'r*','Markersize',12)
saveas(gcf,'OverlappingCount','fig')

Co=count;return;
exit

% Initialize
nsel=npc;nsuby=length(yout);nsubx=length(xout);
%jump=zeros(nsuby,nsubx);jumpstd=zeros(nsuby,nsubx);
iselop=zeros(nsuby,nsubx);

% Preallocate 
isv=unique(idregion);
datarsv(length(isv))=struct('x',[],'y',[],'z',[]);
datamtrsv(length(isv))=struct('x',[],'y',[],'z',[]);

% %%%% coregister DEM strips to a reference DEM tile.
[dX4Sg]=coreg(rang0,idregion2,XYbg2,f2,fdir2);

% %%%% Get initial water mask step 1 %%%


% %%%% Get initial water mask step 2; using the water mask in step 1 as a priori mask %%%


fprintf ('\n Step 1.4: loading of DEM data and assigning chosen output pixels for each group.')
% DEM collection and time series for each zone
%Find the pair closest to epicenter
% [iy,ix]=find((((X-xeq).^2+(Y-yeq).^2))<5e3.);
%Find the maximum overlapping pairs closest to epicenter
[iys,ixs]=find( novlp>=0);
% [iys,ixs]=find( baovlp==1);
clear Xt Yt; for i=1:length(ixs);Xt(i)=X(iys(i),ixs(i));Yt(i)=Y(iys(i),ixs(i));end
[dist,im]=min((Xt-xeq).^2+(Yt-yeq).^2);

isel=im;
display(['nsel=',num2str(nsel)])
tcpu1 = cputime;
ck0=clock;
output=[];
%tag=readGeotiff([macdir,'/home/chunli/scripts/Barnesrock.tif']);
%load([macdir,'/home/chunli/scripts/Barnesrock.mat']);
% load([macdir,'/home/chunli/scripts/BarnesrockRGI.mat']); %tag
ck1=clock;
%parpool(4)
%parfor isel=1:60
for isel=1 %:npc  %nsel 
display(['isel=',num2str(isel)]);
% ix=ixs(isel);iy=iys(isel);
    
% ixy=iy+(ix-1)*ny;%isel=isel+1;
%id=find(range(:,1)<=rangx(ix) & range(:,2)>=rangx(ix+1) & range(:,3)<=rangy(iy)  & range(:,4)>=rangy(iy+1));
% id=idg{ixy};
id=strread(un{isel},'%d');
id=idregion
for j=1:length(id)%0%length(id)
% display(['Overlapping Files Date: ',f{id(j)}(6:13)])
display(['Overlapping Files Date: ',f{id(j)}])
end
%if length(id)<2;continue;end

rang0sov=[max(range(id,1)) min(range(id,2)) max(range(id,3)) min(range(id,4))];
ranget=[rang0sov/resr];
rang0sov=[ceil(ranget(1)) floor(ranget(2)) ceil(ranget(3)) floor(ranget(4))]*resr;
rang0st=[max(rang0sov(1),rang0(1)) min(rang0sov(2),rang0(2)) max(rang0sov(3),rang0(3)) min(rang0sov(4),rang0(4))]; %overlapping grid
x0st=[rang0st(1) rang0st(2) rang0st(2) rang0st(1) rang0st(1) ];y0st=[rang0st(4) rang0st(4) rang0st(3) rang0st(3) rang0st(4) ];

%plot the overlapping zone and two image boundaries
if 0
figure  (1) 
set(gcf,'Color','white')
set(gca,'FontSize', 12);
set(gcf, 'PaperPosition', [0.25 2.5 6 4.5]); %Default [0.25 2.5 8.0 6.0];Prefer 6 by 4.5 or 4 by 3
plot(x0, y0,'g-','linewidth',4)
hold all
% plot(x0s', y0s','r-','linewidth',2)
plot(x0st, y0st,'r-','linewidth',6)
plot(x(id,:)', y(id,:)' ,'b>-') 
plot(xeq, yeq ,'r*','Markersize',12)
axis equal
plot(X(idun==isel),Y(idun==isel),'ro') %overlapping points

hold all
plot3(x0st*1e-3, y0st*1e-3,100*ones(size(x0st)),'r-','linewidth',6)
end

%find the overlapping zone for this piece
mx0=find(xout>=rang0st(1) & xout<=rang0st(2) );
my0=find(yout>=rang0st(3) & yout<=rang0st(4) );
nyj=length(my0);nxj=length(mx0);
demg=-9999*ones(nyj,nxj,length(id));

% load of dem files and save the reduced-size version
t=zeros(length(id),1);
	for j=1:length(id)
	ymd=f{id(j)}(6:13);i=id(j); 
   	demdir=fdir{i};
	t(j)=datenum(ymd,'yyyymmdd');
	iisv=find(isv==i);
	if isempty(datarsv(iisv).x) % non exist
	infile= strrep([demdir,'/',f{i}],'meta.txt','dem.tif');
	data=readGeotiff(infile);
    [M,ratio]=mask(infile);M=~M.z;data.z(M)=-9999; % to do, map out set z(M)=-9999;
%    [M,ratio]=maskentropy(infile);
    	if ratio>=18;
	data.z=-9999*ones(size(data.z)); 
	display(['Bad DEM:',infile])
	end %2nd, the bad points ratio has to be less than 2%.
        %downsizing to 40 m from 8 m.
%       tz = imresize(data.z,dsr); % hard to control coodrinates %0.4s
        %make sure coordinates are exactly on 40m*n for easy faster computation
        ranget=[range(i,:)/resr];
        ranget=[ceil(ranget(1)) floor(ranget(2)) ceil(ranget(3)) floor(ranget(4))]*resr;
        tx=ranget(1):resr:ranget(2);ty=ranget(4):-resr:ranget(3);
%         data.z(data.z == -9999) = NaN; % %convert nodata values to nan for interpolation
%        tic; tzt =interp2(data.x,data.y,double(data.z),tx,ty','*linear');toc;%4s
%         tz(isnan(tz))=-9999; %return to -9999
        %the following get the exactly the same as from interp2
        idrxs=find(abs(data.x-ranget(1))<1e-3);idrxe=find(abs(data.x-ranget(2))<1e-3);
        idrys=find(abs(data.y-ranget(4))<1e-3);idrye=find(abs(data.y-ranget(3))<1e-3);
        idrx=idrxs:nsr:idrxe;
        idry=idrys:nsr:idrye;
        dd=[idrx(2:end)-idrx(1:end-1),idry(2:end)-idry(1:end-1)];
        if isempty(dd)||any(abs(dd-nsr)>1e-3);warning('Resized grid is not constant spacing.');end
        tz=data.z(idry,idrx); %0.09s
        datar= struct();
        datar.x=tx;datar.y=ty;  datar.z=tz;
        datarsv(iisv)=datar; %possible crash with parallel

        infile= strrep(infile,'dem.tif','matchtag.tif');
        data=readGeotiff(infile);
%         tz = imresize(data.z,dsr);
        tz=data.z(idry,idrx);
        datarmt= struct();
        datarmt.x=tx;datarmt.y=ty; datarmt.z=tz;
        datamtrsv(iisv)=datarmt; % save it 
	else %load data
        datar=datarsv(iisv);
	end %if
        %assign DEM to the overlapping zone
        idx=find(datar.x>=rang0st(1) & datar.x<=rang0st(2));
        idy=find(datar.y>=rang0st(3) & datar.y<=rang0st(4));
        if (length(idy)~=nyj || length(idx)~=nxj);warning(['Size mismatch, isel:',num2str(isel)]);end
        demg(1:nyj,1:nxj,j)=datar.z(idy,idx);
    
	end %for j
% end of loading

	[~,~,ni]=size(demg);
    	demp=zeros(ni,1);

        for jx=1:nxj
            for jy=1:nyj
                demp(:)=demg(jy,jx,:); % DEM at a point
                idnn=find(demp~=-9999);idout=[];
                idkp= idnn(~ismember(idnn,idout)); % B(~ismember(B,A)) excluding A from B
                if length(idkp) <=2; continue;end
                tkp=t(idkp);[~,idsort]=sort(tkp);
               % % Detection of event time
               epoch=t/yr;  
               eqall=epoch(idkp(idsort));
               dt=max(eqall)-min(eqall); % in year % maybe for later version 
%                if isempty(dt) || dt<1; continue;end
		lenf=length(idkp);
                [~,nct]=size(lenfc{my0(jy),mx0(jx)}); % my0 mx0: id for the output grid yout, xout
                lenfc{my0(jy),mx0(jx)}(:,nct+1)=[lenf dt isel]';
            end %jy
        end %jx
end % for isel
ck2=clock;
dt=ck2-ck1;tsec=dt(4)*3600+dt(5)*60+dt(6);
display(['Loading files take: ',num2str(tsec),' sec'])

% % % Find the reference DEM for absolute reference (ArcticDEM mosaics) % % %


% % % Find the reference DEM or the lowest stage DEM % % % 

% % % End of finding reference DEM.

fprintf ('\n Step 2: Height time series analysis at the gage.')
for isel=1 %1:npc  %nsel 
%   id=strread(un{isel},'%d');%id of overlapping files
% Unique grouping does not work well for rivers. It increases calculation time (clouds detection) when group size is too small (too many images).
% It's designed for optimum grouping for non-redundant coregistration and region analysis (clouds detection).
% It works well when images are not too many, and each group has a decent large area.
    id =idregion;
    fprintf(['\n Processing isel=',num2str(isel)])

    t=zeros(length(id),1);    
for j=1:length(id)%0%length(id)
	ymd=f{id(j)}(6:13);i=id(j); 
	t(j)=datenum(ymd,'yyyymmdd');
display(['Overlapping Files Date: ',f{id(j)}])
end
[~,idsort]=sort(t);id=id(idsort); %sort the id based on time

rang0sov=[max(range(id,1)) min(range(id,2)) max(range(id,3)) min(range(id,4))];
ranget=[rang0sov/resr];
rang0sov=[ceil(ranget(1)) floor(ranget(2)) ceil(ranget(3)) floor(ranget(4))]*resr;
rang0st=[max(rang0sov(1),rang0(1)) min(rang0sov(2),rang0(2)) max(rang0sov(3),rang0(3)) min(rang0sov(4),rang0(4))]; %overlapping grid
x0st=[rang0st(1) rang0st(2) rang0st(2) rang0st(1) rang0st(1) ];y0st=[rang0st(4) rang0st(4) rang0st(3) rang0st(3) rang0st(4) ];

%display('Step 2.1: searching for the best reference DEM.')
%find the reference DEM: 1, in the overlapping zone maximum good points
%find the reference DEM: 1, in the maximum good points 
% 2, Distance between center of overlapping and center of reference DEM shortest
ymd=[]; jref=[];
idd=[];
% nptsub=zeros(length(id),1);nptall=zeros(length(id),1);
nptsubraw=zeros(length(id),1);nptallraw=zeros(length(id),1);dist=nptsubraw;
score=zeros(length(id),3);rockn=zeros(length(id),1);
enlref=0.5;
tic 
p1=[(rang0sov(1)+rang0sov(2))/2 (rang0sov(3)+rang0sov(4))/2 ];
id =idregion;
for j=1:length(id)
ymd=f{id(j)}(6:13);i=id(j); 
%t(j)=datenum(ymd,'yyyymmdd');
data=datarsv((isv==i));
%checking rock coverage
rangeov=[min(data.x) max(data.x) min(data.y) max(data.y)];
if 0
idx=find(tag.x>=rangeov(1) & tag.x<=rangeov(2)); %assume rock data in 40m resolution
idy=find(tag.y>=rangeov(3) & tag.y<=rangeov(4));
rocktag=tag.z(idy,idx);
rockfilter=rocktag&data.z~=-9999;%&abs(dz)>1e-4;
rockn(j)=sum(rockfilter(:));
end
% Checking the overlapping coverage
idx=find(data.x>=rang0sov(1) & data.x<=rang0sov(2));
idy=find(data.y>=rang0sov(3) & data.y<=rang0sov(4));
demo=data.z(idy,idx);
%nptsubrt=sum(sum(demo~=-9999))*(str2double(res))^2/(dx*dx*4);
nptsubrt=sum(sum(demo~=-9999))/(length(demo(:)));

% ratio of good DEMs for whole scene
XYbi=XYbg{i};
Xb=XYbi(:,1);Yb=XYbi(:,2);
    res=data.info.map_info.dx;
	idx=round((Xb-data.x(1))/res)+1;
	idy=round(-(Yb-data.y(1))/res)+1;
	in= poly2mask(idx,idy, length(data.y),length(data.x)); % faster than inpolygon
    	nptsubrt=sum(sum(data.z~=-9999))/(sum(in(:)))
nptsubrt=sum(sum(data.z~=-9999))/(sum(in(:)));

nptsubraw(j)=nptsubrt; %in cases, all got skipped
nptallraw(j)=sum(sum(data.z~=-9999));
p2=[(data.x(1)+data.x(end))/2 (data.y(1)+data.y(end))/2 ];
dist(j)=sqrt((p1(1)-p2(1))^2+(p1(2)-p2(2))^2);
score(j,1)=nptsubrt*10; %scale 0 to 10
%data=readGeotiff(infile);
if 0
hills=hillshade(double(data.z),data.x,data.y,'plotit');
set(gcf,'Color','white');set(gca,'FontSize', 12);set(gcf, 'PaperPosition', [0.25 2.5 4 3]);
colorbar;%caxis([0 250])
[X,Y]=meshgrid(data.x,data.y);
hold on
plot(x0st, y0st,'r-','linewidth',6)
plot(xeq, yeq ,'r*','Markersize',12)
axis equal
% plot(X(M),Y(M),'r.')
% %saveas(gcf,['refDEMj',num2str(j)],'fig')
% saveas(gcf,'tt','fig')
end
%No edge within the subzone
% 
%if nptsubrt > 0.99 && ratio<1; jref=j;break;end  %hi, change 8% to 1%

% lowest stage 

end
% in cases, all got skipped
%nptsubraw(rockn<20)=0;%if no rock, not consider this dem
idmaxsub=find(abs(nptsubraw-max(nptsubraw)) <0.1); %0.01
% [~,jref]=min(dist(idmaxsub));
% [~,jref]=max(nptallraw(idmaxsub));
score(:,2)=(1-dist/max(dist))*10; %scale 0 to 10;
score(:,3)=(nptallraw/max(nptallraw))*10; %scale 0 to 10;
[~,idj]=max(score(idmaxsub,2)+score(idmaxsub,3));
jref=idmaxsub(idj);

if isempty(jref);warning(['No reference DEM found; isel=',num2str(isel)]);continue;end
j=jref;  iref=id(j);  %store i before idd
id(idd)=[];%nptsub(idd)=[];nptall(idd)=[];
display('Searching reference DEM ...');
toc % 18 minutes
if length(id)<2;continue;end
% jref=1; %datestr(ts)

% 3rd considering the overlapping area with respect to others.
% if refflag==0;warning(['No good reference DEM found; isel=',num2str(isel)]);continue;end

% % %
fis=f{id};
fdiris=fdir{id};
XYbis=XYbg{id};
[Co]=riverprofsub(odir,isel,XYbis,fis,fdiris,xeq,yeq);

% jref=2;%2 2012 Kamchatka Volcano% 1 Volcano %3 2015/10/17 landslide; %1 eq %2 %landslide
i=iref;
if 0
dzxyt=dzxy{idregion==i};
if isempty(dzxyt)
regflag=0; %0 bad 1 good
dxyzref=zeros(1,3);
else
dxyzref=[dzxyt(2:3) dzxyt(1)];
regflag=1;
end
end
%if refflag==0 
data0r=struct();
data0r=datarsv((isv==i));
if 0 %plot the reference DEM
hills=hillshade(double(data0r.z),data0r.x,data0r.y,'plotit');
set(gcf,'Color','white');set(gca,'FontSize', 12);set(gcf, 'PaperPosition', [0.25 2.5 4 3]);
colorbar;%caxis([0 250])
hold on
plot(x0st, y0st,'r-','linewidth',6)
plot(xeq, yeq ,'r*','Markersize',12)
axis equal
title('RefDEM')
saveas(gcf,'refDEM','fig')
end
%find the overlapping zone for this piece
mx0=find(xout>=rang0st(1) & xout<=rang0st(2) );
my0=find(yout>=rang0st(3) & yout<=rang0st(4) );
%reference dem at this piece
datadsx=data0r.x; datadsy=data0r.y ;
idx=find(datadsx>=rang0st(1) & datadsx<=rang0st(2));
idy=find(datadsy>=rang0st(3) & datadsy<=rang0st(4));
nyj=length(idy);nxj=length(idx);
if (length(my0)~=nyj || length(mx0)~=nxj);warning(['Size mismatch, isel:',num2str(isel)]);end
%xout(1:nxj,isel)=datadsx(idx)';     yout(1:nyj,isel)=datadsy(idy)'; 
yout(my0)=data0r.y(idy);xout(mx0)=data0r.x(idx);
%collecting reference dem 
[X,Y]=meshgrid(datadsx(idx),datadsy(idy));
[LAT,LON]=polarstereo_inv(X,Y,6378388.0,0.0819919,70,-45);
% wmp1 = roipoly;data0r.z(wmp1)=NaN; data0r.z(wmp1)=-9999; 
% temporary output of dem
demo=data0r.z(idy,idx);
% output=[output;LAT(:),LON(:),double(demo(:))];
% save -ascii demTyndalbig.dat output
%  continue
latout(my0,mx0,1)=LAT;latout(my0,mx0,2)=LON;latout(my0,mx0,3)=double(demo);

continue

%fprintf ('\n Step 2.2: Coregistration.')
tic
t=zeros(length(id),1);%clear demg demgmt
demg=-9999*ones(nyj,nxj,length(id));demgmt=zeros(nyj,nxj,length(id));
idd=[];
% for j=[3,7,9,2] %2015 10 17 landslide
for j=1:length(id)%[3,5,24,26] % 2012 Kamchatka Volcano 
% for j=1:length(id)  if j==jref;continue;end
    i=id(j);
%     figure (1)
%     hold all
%     plot(x(i,:), y(i,:),'c>-','linewidth',2)
%     plot(x(i,:), y(i,:),'c>-','linewidth',4)
%     axis equal
%     pause
% end
	datar= struct();
	datar=datarsv((isv==i));
    
    % steps for preprocessing (getting edgemask)
    % end of filtering edgemask
    
    idcom=[];out=[];
    if ~isempty(idcom) %exist
%      iter= dzsv{idcom}.iter;
%      out{1,1}.z=dzsv{idcom}.ztar;
%      out{1,1}.x=dzsv{idcom}.x;
%      out{1,1}.y=dzsv{idcom}.y;
    elseif(i==iref)
       iter= 1;
       out{1,1}.z=data0r.z;
       [X,Y]=meshgrid(data0r.x,data0r.y);
       out{1,1}.x=X(:);
       out{1,1}.y=Y(:); 
       dz=zeros(size(out{1,1}.z));
       out{1,2}=out{1,1};out{1,2}.z=ones(size(out{1,1}.z));
       rockfilter=false(size(out{1,1}.z));
    else  % do coregistration
        % get the overlapping target and reference DEM subzones, and the
        % rock flag
        rangref=[min(data0r.x) max(data0r.x) min(data0r.y) max(data0r.y)];
        rangtar=[min(datar.x) max(datar.x) min(datar.y) max(datar.y)];
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
        if (size(tardem.z)~=size(refdem.z));warning(['Wrong overlapping crop, check isel:',num2str(isel)]);end
	if 0
%           rocktag = interp2(tag.x,tag.y,tag.z,datar.x,datar.y','*nearest',0);
        idx=find(tag.x>=rangeov(1) & tag.x<=rangeov(2)); %assume rock data in 40m resolution
        idy=find(tag.y>=rangeov(3) & tag.y<=rangeov(4));
        rocktag=tag.z(idy,idx);
        if (size(rocktag)~=size(refdem.z));warning(['Wrong overlapping crop, check isel:',num2str(isel)]);end
        end
 
        if coregflag==1 % minus rock mean
         %convert nodata values to nan for interpolation
%           refz=data0r.z;refz(refz==-9999)=NaN;%avoid -9999 minus a dem or -9999 minus -9999
          % avoid interpolation for speed
%           demapr = interp2(data0r.x,data0r.y,refz,datar.x,datar.y','*linear',NaN);
          datatarz=tardem.z;datatarz(datatarz== -9999) = NaN; 
          datarefz=refdem.z;datarefz(datarefz== -9999) = NaN;
          iter= 1;
          %avoid -9999 - meanhrock
          out{1,1}.z=datatarz;%initial
          [X,Y]=meshgrid(tardem.x,tardem.y);
          out{1,1}.x=X(:);
          out{1,1}.y=Y(:);
          out{1,2}=out{1,1};out{1,2}.z=ones(size(out{1,1}.z));
          out{1,3}=out{1,1};out{1,3}.z=refdem.z;	
          dz = datatarz-datarefz;
          rockfilter=rocktag&~isnan(dz)&refdem.z~=-9999&abs(dz)<50;%&abs(dz)>1e-4;
          dzrock=dz(rockfilter);
          if isempty(dzrock)
              meanhrock=0; iter= 49;
          else
              meanhrock=mean(dzrock);
          end
          out{1,1}.z=out{1,1}.z-meanhrock;dz=dz-meanhrock;
	elseif coregflag==2 % use reg.txt file
          iter= 1;
          out{1,1}.z=tardem.z;%initial
          [X,Y]=meshgrid(tardem.x,tardem.y);
          out{1,1}.x=X(:);
          out{1,1}.y=Y(:);
          out{1,2}=out{1,1};out{1,2}.z=ones(size(out{1,1}.z));
          out{1,3}=out{1,1};out{1,3}.z=refdem.z;
          dzxyt=dzxy{idregion==i};
          if regflag==0 || isempty(dzxyt)
             iter=49;
             dxyz=zeros(1,3);
          else
              dxyztar=[dzxyt(2:3) dzxyt(1)];
              dxyz=dxyztar-dxyzref; % read from reg.txt file
          out{1,4}.P{1,2}(1,1)=1; out{1,4}.P{1,2}(2:4,1)=0;
          out{1,4}.P{1,2}(5:7,1)=-dxyz;
	%in coregisterdems2, tardem+txyz=refdem;in coregisterdems and reg.txt, tardem-dxyz=refdem;
          [outx]=tarx(out,tardem,refdem,[resr resr],params);
          out=outx;
          end
          dz=out{1,1}.z-out{1,3}.z;
        elseif coregflag==3 %use mp1 mp2 with simple coregistration
%         mp1 = rocktag;%interp2(tag.x,tag.y,tag.z,data0r.x,data0r.y','*nearest',0);
%         mp2 = rocktag;
          datatarz=tardem.z;datatarz(datatarz== -9999) = NaN; 
          datarefz=refdem.z;datarefz(datarefz== -9999) = NaN;
          iter= 1;
%         [z2out,p,d0] = coregisterdems(refdem.x,refdem.y,double(datarefz),tardem.x,tardem.y,double(datatarz),mp1,mp2);
          [z2out,p,d0] = coregisterdems(refdem.x,refdem.y,double(datarefz),tardem.x,tardem.y,double(datatarz));
          if size(z2out)~=size(refdem.z);warning('coregistration failure');iter= 49; end
          dz = z2out-refdem.z;
          out{1,1}.z=z2out;
          [X,Y]=meshgrid(refdem.x,refdem.y);
          out{1,1}.x=X(:);
          out{1,1}.y=Y(:);
          out{1,2}=out{1,1};out{1,2}.z=ones(size(out{1,1}.z));
          out{1,3}=out{1,1};out{1,3}.z=refdem.z; 
%         rockfilter=mp1&(refdem.z~=-9999);
%        rockfilter=mp1&~isnan(data0r.z);
        elseif coregflag==7 
            if 0
            [out] = coregisterdems2(tardem, refdem, [resr resr],params);% best results
            else % use rock 
%         [idy,idx]=find(rocktag==1); %assume rock data in 40m resolution
        [idy,idx]=find(rocktag==1&refdem.z~=-9999); %assume rock data in 40m resolution
        idy=min(idy):max(idy);idx=min(idx):max(idx);
        tardemr.z=tardem.z(idy,idx);tardemr.x=tardem.x(idx);tardemr.y=tardem.y(idy);
        refdemr.z=refdem.z(idy,idx);refdemr.x=refdem.x(idx);refdemr.y=refdem.y(idy);
        [out] = coregisterdems2(tardemr, refdemr, [resr resr],params);
%         out{1,4}.P{1,2}(1)=1; out{1,4}.P{1,2}(2:4)=0;% in case using parameters from small area to large area
        %out{1,4}.P{1,2}(5:7)=-[1.1020,   39.4060 , 4.3870]';
        [outx]=tarx(out,tardem,refdem,[resr resr],params);
        out=outx;
        end 
            
            iter=out{1,4}.S{1,1};
            % control points percentage needs to > 0.3%
            [m,n]=size(out{1,2}.z);npt=m*n;
            ncs=out{1,4}.S{1,2}(1); %(npt-sum(out{1,2}.z(:)))
            RC=ncs/npt*100;
            display(['RC=',num2str(RC)])
            if RC < 0.0001;iter=49;end
            
            dz=out{1,1}.z-out{1,3}.z;
        end % if coregflag
    end %if exist
     %check the quality
    if iter==49
        idd=[idd;j];
        continue
    end
    z2out=out{1,1}.z;
   % assign to the subzone
    tmpy=reshape(out{1,1}.y,size(out{1,1}.z));tmpx=reshape(out{1,1}.x,size(out{1,1}.z));
    datadsx=tmpx(1,:)';datadsy=tmpy(:,1); % to check the space
	idx=find(datadsx>=rang0st(1) & datadsx<=rang0st(2));
	idy=find(datadsy>=rang0st(3) & datadsy<=rang0st(4));
%   z2inp = interp2(datadsx,datadsy,z2out,xoutt(mx0),youtt(my0)','*nearest');    
%   demg(my0,mx0,j)=z2inp;
if (length(idy)~=nyj || length(idx)~=nxj);warning(['Size mismatch, isel:',num2str(isel)]);end
    datadsx=xout(mx0);datadsy=yout(my0);
    demg(1:nyj,1:nxj,j)=z2out(idy,idx);
    ymd=f{i}(6:13);
    t(j)=datenum(ymd,'yyyymmdd');
    
        %matchtag file as weight
	data= struct();
	data=datamtrsv((isv==i));
    
    % interpolate to match the coregistered DEM
        mt = interp2(data.x,data.y,data.z,datadsx,datadsy','*nearest');    
        demgmt(1:nyj,1:nxj,j)=mt;%mt(idy,idx);
    
    if 0 % plot orthoimage
%     if t(j) < teq && ~any(isv==i)
        infile=strrep(metaFile,'meta.txt','ortho.tif');
        data=readGeotiff(infile);
        if j==jref % reduce the size to be the same as reference dem
        data.z = imresize(data.z,0.2);% 8m to 40m
        data.x=data0r.x;data.y=data0r.y;
        end
        h=double(data.z);%h=I.z;h=m.z;
        minh=min(min(h(h~=0)));maxh=max(max(h(h~=0)));
        meanh=mean(mean(h(h~=0)));stdh=std(h(h~=0));
        %DEM  zone
        figure;set(gcf,'Color','white');set(gca,'FontSize', 12);set(gcf, 'PaperPosition', [0.25 2.5 6 4]); 
        % surf(data.x*1e-3,data.y*1e-3,data.z);colormap gray; colorbar; shading interp; axis tight; view(0,90)  
        imagesc(data.x*1e-3,data.y*1e-3,data.z);colormap gray; colorbar; shading interp; axis tight; view(0,-90)  
        hold all
        colorbar;   caxis([max(meanh-stdh, minh) min(meanh+1*stdh,maxh)]);
        axis equal
        set(gca,'FontSize', 12);
        title(['Over lapping Orthoimage ',num2str(j),' at zone ', num2str(isel)]);
        plot(x0st*1e-3, y0st*1e-3,'r-','linewidth',4)
        axis equal
        set(gca,'FontSize', 12);
        [dir,str,ext] =fileparts(infile);
        ofile=['SAlaska/',str,'Zone', num2str(ixy),'ov',num2str(j)];
%         mp1 = roipoly;
%         mp2 = roipoly;
    %     ofile=['Zone', num2str(ixy),'DEM',num2str(j)];
    %       saveas(gcf,ofile,'fig')
%         print('-dpdf','-r300',ofile) 
      end
      if 0 % plotting DEMs
        infile= strrep(metaFile,'meta.txt','dem.tif');
        data=readGeotiff(infile);
        hills=hillshade(double(data.z),data.x,data.y,'plotit');
        set(gcf,'Color','white');set(gca,'FontSize', 12);set(gcf, 'PaperPosition', [0.25 2.5 4 3]);
        colorbar;%caxis([0 250])
        title(['Over lapping DEM ',num2str(j),' at zone ', num2str(isel)]);
        hold all;plot(x0st, y0st,'r-','linewidth',4)
        hold on %control points
plot(out{1,2}.x(~out{1,2}.z(:)),out{1,2}.y(~out{1,2}.z(:)),'r>','Markersize',2,'Linewidth',2)
      end
      if 0 %Plot DEM difference
%         figure;imagesc(data0r.x*1e-3,data0r.y*1e-3,dz,'alphadata',~isnan(dz) & abs(dz) < 10)
%         figure;imagesc(data0r.x*1e-3,data0r.y*1e-3,dz,'alphadata',~isnan(dz))
        figure;imagesc(out{1,1}.x*1e-3,out{1,1}.y*1e-3,dz,'alphadata',~isnan(dz) & abs(dz) < 100)
        set(gcf,'Color','white');set(gca,'FontSize', 12);set(gcf, 'PaperPosition', [0.3 0 3.8 3]); set(gcf, 'PaperSize', [4 3]);
        colormap jet; axis equal; view(0,-90)
        xlabel('{\itx} (km)')
        ylabel('y (km)')        
        hold all
        colorbar
        caxis([-3 3])
        plot(x0st*1e-3, y0st*1e-3,'r-','linewidth',4)
        hold on %control points
        plot(X(rockfilter)*1e-3,Y(rockfilter)*1e-3,'k.')
% plot(out{1,2}.x(~out{1,2}.z(:))*1e-3,out{1,2}.y(~out{1,2}.z(:))*1e-3,'r>','Markersize',8,'Linewidth',2)
        box on        
title([f{i}(6:13),'-',f{iref}(6:13)])
      end
        if 0
        figure;imagesc(out{1,2}.x*1e-3,out{1,2}.y*1e-3,out{1,2}.z)
        set(gcf,'Color','white');set(gca,'FontSize', 12);set(gcf, 'PaperPosition', [0.3 0 3.8 3]); set(gcf, 'PaperSize', [4 3]);
        colormap gray; axis equal; view(0,-90)
        xlabel('{\itx} (km)')
        ylabel('y (km)')        
        hold all
        colorbar
        caxis([0 1])
        plot(x0st*1e-3, y0st*1e-3,'r-','linewidth',4)
        box on    
        title(['Control points']) % 0 is control point
%         figure;imagesc(out{1,1}.x*1e-3,out{1,1}.y*1e-3,out{1,1}.z)
%         title(['Coregistered DEM'])
      end
%   isv=[isv,i];
end
% delete the bad coregistered DEMs
id(idd)=[];t(idd)=[];
demg(:,:,idd)=[];demgmt(:,:,idd)=[];
if length(id)<2;continue;end
toc  % loading 34 files take 8 minutes, too slow.
% coregistering 17 files take 29 minutes.
% end of Gathering DEMs
% pause
% close all
% h = get(0,'children');
% for i=1:length(h)
%   saveas(h(i), ['figure' num2str(i)], 'fig');
% end

%fprintf ('\n Step 2.3: Time series analysis.')
    %Preparing the time series analysis
    epochorg=zeros(length(t),1);
    epochorg(:)=t(:);
    maxcount=0;
    algorithm=2; % 1 linear; 2 constant
%   yr=365.25;
    epoch=epochorg/yr;  % Change the unit of time from day to year, to avoid singular matrix
    tm=mean(epoch);
    
    [~,~,ni]=size(demg);
    demg(isnan(demg))=-9999;
    % Time series analysis
    tic
    demp=zeros(ni,1);dempmt=ones(ni,1);
    AMa=ones(ni,1);
    stdmax=40;
%     jx=find(abs(datadsx(idx)-xeq)<1);jy=find(abs(datadsy(idy)-yeq)<1);
    jumpt2=zeros(nyj,nxj);jumpstdt2=zeros(nyj,nxj);
    for jxy=1:length(idop)
            mxp=ceil(idop(jxy)/nsuby);myp=idop(jxy)-(mxp-1)*nsuby;
            jy=find(my0==myp);jx=find(mx0==mxp);
            if length(jy)~=1 ||length(jx)~=1 ; warning('wrong index'); end
            demp(:)=demg(jy,jx,:); % DEM at a point
            dempmt(:)=demgmt(jy,jx,:);
            Pd=ones(ni,1);Pd(dempmt==0)=0.01; %lesser weight at false match points
            idnn=find(demp~=-9999);idout=[];idoutpre=[];
            idkp= idnn(~ismember(idnn,idout)); % B(~ismember(B,A)) excluding A from B
%             idkp=[2,4,7,8]'; % n by1
%             idkp=[1,3];
%             idkp=[2,3,6,7]';
            if length(idkp) <=2; continue;end
            
            %Outlier detection 
	    % Change to no iteration, in case of big jump events
            %%%  Linear fitting model: y=a + b(t-tm) 
    	    algorithm=2; % 1 linear; 2 constant
             if (algorithm == 1 ) %linear trend
                mp=2; 
             elseif (algorithm == 2) % constant
                mp=1;
             end
            
            minnpt=mp;
            epoch=epochorg/yr;
            % detecting outliers
            count=0;ndw=0;flagcond=0;
            while 1
            idkp= idkp(~ismember(idkp,idout));
            tkp=t(idkp);[ts,idsort]=sort(tkp);
            if length(idkp)<=minnpt; break;   end
            P=diag(Pd(idkp)); yobs=demp(idkp);
            
            AM=[];
            AM=zeros(length(idkp),mp);  %
            AM(:,1)=AMa(idkp,:);
            if (algorithm == 1 ) %linear trend            
            AM(:,2)=epoch(idkp)-tm;
            elseif (algorithm == 2) % constant
                %do nothing
            end

            cdA=cond(AM'*P*AM);
            if(cdA>1e5) 
display(['condtion number of AM*P*AM is large: ',num2str(cdA),';isel:',num2str(isel)]);
                flagcond=1;
                break
            end
            est=inv(AM'*P*AM)*AM'*P*yobs;
%             fit=est(1) + est(2)*(epoch(idkp)-tm);
            fit=AM*est;
            etilde=yobs-AM*est;
            sigma02hat=etilde'*P*etilde/(length(yobs)-1);
            meani=mean(fit);  %est(1); bad value for less keeped points and a bad slant angle. 
            stdi=[];
            stdi=sqrt(sigma02hat);
%             if stdi >stdmax;multi=2;else; multi=3;end
            multi=3;
            id=find(abs(etilde)>=multi*stdi);idout=[idoutpre;idkp(id)];idsign=id;
            if isempty(idsign); ndw=ndw+1;else ndw=0;end
%             id=find(abs(etilde)>=2*stdi);Pd(idkp(id))=(1/3.)^2;Pd(dempmt==0)=0.01;  %Algorithm 4:down weight the points far from mean.%
            count=count+1;
            if 0 % plotting
            figure % (1)
            set(gcf,'Color','white')
            set(gca,'FontSize', 18);
            set(gcf, 'PaperPosition', [0.25 2.5 4 3]); %Default [0.25 2.5 8.0 6.0];Prefer 6 by 4.5 or 4 by 3
            hold all
            plot(t(idkp(idsort)),demp(idkp(idsort)),'b>-','MarkerSize',12,'linewidth',4)
%                     plot(t(idkp),demp(idkp),'b>-','MarkerSize',12,'linewidth',4)
            plot(t(dempmt==0),demp(dempmt==0),'co','MarkerSize',18,'linewidth',4)
            plot(t(idout),demp(idout),'ks','MarkerSize',12,'linewidth',4)
            plot(t(idkp(idsort)),fit(idsort),'g*-','MarkerSize',12,'linewidth',4)
            plot(t(idkp(idsort)),fit(idsort)-multi*stdi,'r-','linewidth',4)
            plot(t(idkp(idsort)),fit(idsort)+multi*stdi,'r-','linewidth',4)
            datetick('x','mm/yy')
            box on
            ylabel('DEM time series') 
            end % plotting
            if isempty(idsign)&&ndw>=2; break;   end %If no outliers detected, and the weight has been adjusted once, stop the iteration
            if count>=2 %99
%                 display(['Too many iterations ',count,'; jx=',num2str(jx),'jy=',num2str(jy)]);
                break
            end
            idoutpre=idout;
%               break % no iteration applied
            end % outliers loop
            %end of outlier detection
            if count>maxcount; maxcount=count;end
            idkp= idkp(~ismember(idkp,idout));
%             idkp= idkp(~ismember(idkp,[6,13,21,23,25]));
            tkp=t(idkp);[ts,idsort]=sort(tkp);
            if flagcond ==1 || length(idkp) <=2
                continue;
            end

            %Initializing
            P=[];
            lenf=length(idkp);
            T6=demp(idkp(idsort));
	    stdi=1;
            T6std=stdi*ones(size(T6));
%             T6std=stdi*Pd(idkp(idsort)).^-0.5;
            P=diag(T6std(:).^-2); 

            % % Detection of event time
           eqall=epoch(idkp(idsort));
           eqs=datestr(min(eqall)*yr,26);eqe=datestr(max(eqall)*yr,26);
           dt=max(eqall)-min(eqall); % in year
           if isempty(dt) || dt<1; continue;end

            %Jump fitting
 %%%  Linear fitting model: y=a + b(t-tm) 
%      +d11*(0 or 1);            
            epoch=epoch(idkp(idsort));
            mpnoj=2; %no jump
	    algorithm = 1;
            
            AM=[];
            AM=zeros(lenf,mp);  %
            AM(:,1)=1.;
            if (algorithm == 1 ) %linear trend            
                AM(:,2)=epoch-tm;
            elseif (algorithm == 2) % constant
                %do nothing
            end
            yobs=T6;

            cdA=cond(AM'*P*AM);
            if(cdA>1e5); display(['condtion number of AM*P*AM is large: ',num2str(cdA),';isel:',num2str(isel)]); continue;end

%             est=var*AM'*P*yobs;
            est=AM\yobs; % equal weight
            
            var=inv(AM'*P*AM);
            %reestimate the reference variance;
            if lenf <=2 
                sigma02hat=NaN;
            else
                etilde=yobs-AM*est;
                sigma02hat=etilde'*P*etilde/(length(yobs)-1);
            end

            var=var*sigma02hat;
            T6std=T6std*sqrt(sigma02hat);
            
            for m=1:mp 
%               eststd(m:m)=sqrt(var(m,m)); 
            end

            fit=AM*est;
            
            jumpt2(jy,jx)=est(mpnoj);
            jumpstdt2(jy,jx)=sqrt(var(mpnoj,mpnoj));       

            
	    if algorithm==1
	     % m/yr
		trend=est(mpnoj);trest=sqrt(var(mpnoj,mpnoj));
        if (trest< 0.8)
          fprintf(fid2,'%12.6f %12.6f  %23.15e %23.15e\n',LAT(jy,jx),LON(jy,jx),trend,trest);
        end %if
            jump(idop(jxy))=trend;
            jumpstd(idop(jxy))=trest;       
        end

            if 0 %plotting
% % For plotting figures       
        plot(x0sov*1e-3, y0sov*1e-3,'m-','linewidth',6)

            end
    end
end %im   or isel    

fprintf ('\n Step 3: Discharge Time series analysis at gage.')


e = cputime-tcpu1
ck2=clock;
dt=ck2-ck0;tsec=dt(3)*86400+dt(4)*3600+dt(5)*60+dt(6);
fprintf (['\n Program takes: ',num2str(tsec),' sec'])
% get river map
% [X,Y]=meshgrid(-141.716:0.001:-141.,60.0833:0.001:60.25);
% figure;imagesc(X(1,:),Y(:,1),ones(size(X)))
% mp1=roipoly;
% output=[ X(mp1) Y(mp1)];
% save -ascii riv.dat output

fid4 = fopen('demTyndalbig.dat','w');
for jx=1:nsubx
for jy=1:nsuby
latt=latout(jy,jx,1);lont=latout(jy,jx,2);
fprintf(fid4,'%12.6f %12.6f %23.15e\n',latt,lont,latout(jy,jx,3));
end %jy
end %jx

figure;set(gcf,'Color','white');set(gca,'FontSize', 18);set(gcf, 'PaperPosition', [0.25 2.5 4 3]); 
hold all
xt=xout*1e-3;yt=yout*1e-3;zt=jump;
imagesc(xt,yt,zt); colormap jet;colorbar; shading interp; axis equal; view(0,90)  
set(gca,'FontSize', 18);
plot(xeq*1e-3, yeq*1e-3 ,'r*','Markersize',12)
xlabel('{\itx} (km)')
ylabel('y (km)')
caxis([-3 3])
axis(rang0*1e-3)
print('-dpdf','-r300','jump')   
saveas(gcf,'jump','fig')

figure;set(gcf,'Color','white');set(gca,'FontSize', 18);set(gcf, 'PaperPosition', [0.25 2.5 4 3]); 
hold all
xt=xout*1e-3;yt=yout*1e-3;zt=jumpstd;
imagesc(xt,yt,zt);
colormap jet;colorbar; shading interp; axis equal; view(0,90)  
set(gca,'FontSize', 18);
plot(xeq*1e-3, yeq*1e-3,'r*','Markersize',12)
xlabel('{\itx} (km)')
ylabel('y (km)')
caxis([0 1])
axis(rang0*1e-3)
print('-dpdf','-r300','jumpstd')   
saveas(gcf,'jumpstd','fig')

exit

jump=load('jump.txt');
xt=jump(:,1)*1e-3;yt=jump(:,2)*1e-3;zt=jump(:,3);ztstd=jump(:,4);
figure;set(gcf,'Color','white');set(gca,'FontSize', 18);set(gcf, 'PaperPosition', [0.25 2.5 4 3]); 
hold all
plot3(xt,yt,zt); colormap jet;colorbar; shading interp; axis equal; view(0,90)  
set(gca,'FontSize', 18);
plot3(xeq*1e-3, yeq*1e-3,1e2 ,'r*','Markersize',12)
xlabel('{\itx} (km)')
ylabel('y (km)')
print('-dpdf','-r300','jump') 


xt=xout(1:nxj,isel)*1e-3;yt=yout(1:nyj,isel)*1e-3;zt=jump(1:nyj,1:nxj,isel);
dz=zt;
M=(~isnan(dz) & abs(dz) < 200);
dz(~M)=0;%NaN;
pt= 8.3829995e+00; -30.64;%-169.6
[jy,jx]=find(abs(dz-pt)< 1e-6);
