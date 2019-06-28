function [Co]=riverprofsub(odir,Msv,XYbis,fis,fdiris,dX4Sg,XYb2is,f2is,fdir2is,dX4Sg2)
%function [Co]=riverprofsub(odir,wm,Msv,XYbis,fis,fdiris,dX4Sg,XYb2is,f2is,fdir2is,dX4Sg2)
% method 1: direct method; % has to be scene files.
% method 2: imagery-altimetry method % can be scene files or mono imageries with one reference DEM.

% Elevation extraction.
%
% Chunli Dai, July 2017
% Chunli Dai, December 2017
% Chunli Dai, November 2018
%Gage at Fairbanks
%filename='filelist';
constant
Co=[];

iref=1; %20111008 the lowest stage. 


%fprintf ('\n Step 0: geting the boundary for all files in the region.')

n=length(fis);n2strip=length(f2is);
id=1:n;

if 0 %flagcoreg==1
%control/stable surfaces for coregistration; roads database
riv=load('FAirport2.gmt');
[xap,yap]=polarstereo_fwd(riv(:,2),riv(:,1),[],[],70,-45);
park=load('park.gmt');park2=load('park2.gmt');road=load('road.gmt');
[xpark,ypark]=polarstereo_fwd(park(:,2),park(:,1),[],[],70,-45);
[xpark2,ypark2]=polarstereo_fwd(park2(:,2),park2(:,1),[],[],70,-45);
[xroad,yroad]=polarstereo_fwd(road(:,2),road(:,1),[],[],70,-45);
end % if 

% For strip files using the direct method.
datarsv2(n2strip)=struct('x',[],'y',[],'z',[]);
idd=[];
for j=1:n2strip
stripmetafile= [fdir2is{j},'/',f2is{j}];
%WV02_20160515222857_1030010054295A00_16MAY15222857-M1BS-500728930090_01_P003.tif
texttar=f2is{j}(1:19);
[datestr]=strip2date(stripmetafile);
texttar(6:19)=datestr(:);
p=dX4Sg2(j,:);
fprintf(['\n Work on j p stripmetafile ',num2str([j p]),' ',stripmetafile])

% % % Get water mask % % %
try
M=multispecstrip(stripmetafile,p,Msv,fis,dX4Sg); %Notice: the matrix may not match the DEM strip matrix; due to different resolution. 
%M=multispecstrip(stripmetafile,p,wm,fis,dX4Sg,Msv); 
catch e
     fprintf('\nThere was an error! The message was:\n%s',e.message);
     fprintf(['\nstripmetafile is ',stripmetafile])
save test1.mat stripmetafile p Msv fis dX4Sg  -v7.3
idd=[idd;j];
continue
end

%strip DEM only covers the overlapping area of two mono images.
if isempty(M.z);
 idd=[idd;j];
 continue
end
M1=M; M1= rmfield(M1,'coast');
datarsv2(j)=M1;
% M=Msv{j}; %maskentropy(infile);

M=M.coast;%M.z;
y0s=0;%y0;
ymd=str2num(f2is{j}(6:13));
satname=f2is{j}(1:4);

%p=dX4Sg2(j,:);
     
infile= strrep(stripmetafile,'meta.txt','dem.tif');
data=readGeotiff(infile);
%filtering out bad edges, e.g. 20110804
% Me=~(imdilate((data.z==-9999),ones(round(30*8)))); %Bug 23, remove too much river area
%data.z(~Me)=-9999; %April 2019, pzxy may change the value of -9999.
%edge based on boundary files.
dx=data.x(2)-data.x(1); dy=data.y(2)-data.y(1);[nwy,nwx]=size(data.z);
sx=XYb2is{j}(:,1); sy=XYb2is{j}(:,2); 
idx=round((sx-data.x(1))/dx)+1;idy=round((sy-data.y(1))/(dy))+1;
Med=~poly2mask(idx,idy,nwy,nwx); % fast, apply to each polygon one by one.
Me=~(imdilate(Med,ones(round(30*8)))); 
data.z(~Me|data.z==-9999)=NaN;

[m1,n1]=size(M);[m2,n2]=size(data.z);
if m1~=m2||n1~=n2 %Align/Interpolate the water mask to strip DEM.
Mt= interp2(M1.x,M1.y,double(M),data.x,data.y','*nearest',0);
M=logical(Mt);
end

mtFile = strrep(infile,'dem.tif','matchtag.tif');
mt=readGeotiff(mtFile);

% Align the DEM coordinates to reference DEM using the translation parameters from coregistration.
p
data.x=data.x- p(2);data.y=data.y- p(3);data.z=data.z- p(1);
[X,Y]=meshgrid(data.x,data.y);
[LAT,LON]=polarstereo_inv(X(M),Y(M),[],[],70,-45);

if 0 %flagplot==1 %too big
hills=hillshade(double(data.z(1:10:end,1:10:end)),data.x(1:10:end),data.y(1:10:end),'plotit');
hold on; plot(X(M),Y(M),'r.')
ofile=[odir,'/rivlocd',num2str(ymd),'j',num2str(j)];
title(texttar)
saveas(gcf,ofile,'fig')
end

%if method==1 || 1
%Output the extracted elevation along shorelines
output=[LAT(:),LON(:),double(data.z(M)),X(M),Y(M)];
%ofile=[odir,'/rivprof',num2str(ymd),satname,'rawj',num2str(j),'.dat'];
ofile=[odir,'/rivprof',texttar,'rawj',num2str(j),'.dat'];
save(ofile,'output','-ascii') 

%save M water mask
projstr='polar stereo north';
ofile=[odir,'/watermask',texttar,'sj',num2str(j),'.tif'];
%writeGeotiff(ofile,data.x,data.y,uint8(M),1,255,projstr) %wrong, M is just shoreline, not water mask.
writeGeotiff(ofile,M1.x- p(2),M1.y - p(3),uint8(M1.z),1,255,projstr)
fprintf(['\n Double check j p:',num2str([j p])])

% [X,Y]=meshgrid(data.x,data.y);
% data.z=z;
demp=data.z(M);
dempmt=mt.z(M);
epochorg=Y(M);
[idkp ]=outlier(odir,demp,dempmt,epochorg,ymd);
yp=epochorg(idkp);

output=[LAT(idkp),LON(idkp),double(demp(idkp)),(yp(:)-y0s)];
if isempty(output) %avoid writing empty files
continue;
end
% ofile=[odir,'rivprof',num2str(ymd),'s.dat']; %remove outlier
%ofile=[odir,'/rivprof',num2str(ymd),satname,'sj',num2str(j),'.dat'];
ofile=[odir,'/rivprof',texttar,'sj',num2str(j),'.dat'];
save(ofile,'output','-ascii') 
%output=[epoch(t:)-y0s,T6(:),T6std(:)];
%ofile=['rivprof',num2str(i),'ms.dat']; %mean and std
%save(ofile,'output','-ascii') 
%end % if method ==1
end %j
%delete strips that have empty water mask.
if length(idd)>0
    fprintf(['\n ',num2str(length(idd)),' strip DEMs have no water mask data:',f2is{idd},'\n'])
    if n2strip > length(idd) %there is strip files left.
        XYb2is(idd)=[];f2is(idd)=[];fdir2is(idd)=[];dX4Sg2(idd,:)=[];datarsv2(idd)=[];
        flagleft=1; %there is good strip DEM left
    else
        flagleft=0;% all strips bad.
    end
    n2strip=length(f2is);
end
if flagleft==1
fprintf(['\n Use these strip DEMs (have water mask data):',num2str(n2strip),' ',f2is{1:n2strip},'\n'])
elseif flagleft==0
    fprintf('\n None strip DEM has water mask. \n')
    fprintf(['\n Use these strip DEMs :',num2str(n2strip),' ',f2is{1:n2strip},'\n'])
end
% % % Find the lowest stage DEM in all strips % % % 
[flagstats1,idref1]=lowest(datarsv2,XYb2is);
if ~isempty(idref1)
idref1=idref1(1); %idref1 could be a list of sorted id.
end
% load data
if flagstats1==1
    k=idref1;
    stripmetafile= [fdir2is{k},'/',f2is{k}];
    p=dX4Sg2(k,:);     
    infile= strrep(stripmetafile,'meta.txt','dem.tif');
fprintf (['\n The lowest stage DEM in all strips:',infile])

    data=readGeotiff(infile);
    %filtering out bad edges, e.g. 20110804
%     Me=~(imdilate((data.z==-9999),ones(round(30*8)))); 
%data.z(~Me)=-9999; %April 2019, pzxy may change the value of -9999.
    dx=data.x(2)-data.x(1); dy=data.y(2)-data.y(1);[nwy,nwx]=size(data.z);
    sx=XYb2is{k}(:,1); sy=XYb2is{k}(:,2); 
    idx=round((sx-data.x(1))/dx)+1;idy=round((sy-data.y(1))/(dy))+1;
    Med=~poly2mask(idx,idy,nwy,nwx); % fast, apply to each polygon one by one.
    Me=~(imdilate(Med,ones(round(30*8)))); 
    data.z(~Me|data.z==-9999)=NaN;
    
    mtFile = strrep(infile,'dem.tif','matchtag.tif');
    mt=readGeotiff(mtFile);
    % Align the DEM coordinates to reference DEM using the translation parameters from coregistration.
    data.x=data.x- p(2);data.y=data.y- p(3);data.z=data.z- p(1);

    data0a=data;mt0a=data;mt0a.z=mt.z;
else
    fprintf (['\n The lowest stage DEM in all strips not found.'])
end
% % % End of finding reference DEM.

res=2;

% For mono images using the altimetry-imagery method
for j=1:n%0%length(id)
i=j;
texttar=fis{j}(1:19);
[~,filename,~]=fileparts([fdiris{j},'/',fis{j}]);
display(['Working on mono image ',num2str(j), ': ',fdiris{j},'/',fis{j}])

% % % Get water mask % % %
M=Msv(j); %maskentropy(infile);

if j==137
save testj137.mat M datarsv2 XYb2is -v7.3
end

if isempty(M.z)
   display(['Empty mono image ',num2str(j), ': ',fdiris{j},'/',fis{j}])
   continue
end

p=dX4Sg(j,:) 
[X,Y]=meshgrid(M.x- p(2),M.y - p(3));
%[LAT,LON]=polarstereo_inv(X(M1),Y(M1),[],[],70,-45);

%save M water mask
projstr='polar stereo north';
ofile=[odir,'/watermask',texttar,'bj',num2str(j),'.tif'];
writeGeotiff(ofile,M.x- p(2),M.y - p(3),uint8(M.z),1,255,projstr)

rang0=[min(M.x) max(M.x) min(M.y) max(M.y)];

% demdir=fdiris{j};
% infile= strrep([demdir,'/',fis{j}],'meta.txt','dem.tif');
%WV01_20080514213507_10200100021A7A00_08MAY14213507-P1BS-052804651020_01_P007_u16ns3413.xml
ymd=str2num(fis{j}(6:13));
satname=fis{j}(1:4);

mon=str2num(fis{j}(10:11));
if mon <=4 || mon>=11  %winter
    flagwin=1; %1 winter
    display(['Warning: winter scenes j=',num2str(j),';ymd=',num2str(ymd)])
else
    flagwin=0; %0 summer
end

dx=M.x(2)-M.x(1); dy=M.y(2)-M.y(1);
% % % Check whether the lowest stage DEM covers this image.% % % 
M1=double(M.z);M1(M1==-1)=0;
%to do: apply the river centerline map, make sure the coverage of DEM is calculated using the river rather than all water areas.
%done: in maskentropy.m, the river mask is the water mask within the buffer zone.
nwx=length(M.x);nwy=length(M.y);
%use the manually selected DEM;
if flaglowest==1 %use the given dem;
    iref=idlowest;
    fprintf (['\n Use the manually given DEM as lowest stage DEM:',f2is{iref},'\n'])
else %flaglowest = 0 or 2
    if flaglowest == 0 
        
if flagstats1==1
    k=idref1;
    % Use strip dem to get shoreline elevation; need water mask to determine its water stage.
    %Sometimes part of strip has no water data in it, but here it's only checking the coverage of dem.
    sx=XYb2is{k}(:,1); sy=XYb2is{k}(:,2); 
    idx=round((sx-rang0(1))/dx)+1;idy=round((sy-M.y(1))/(dy))+1;
    sm=poly2mask(idx,idy,nwy,nwx); % fast, apply to each polygon one by one.
    %sm:  1 has data; 0 has no dem data.
    sm=double(sm);
    ratio=sum(sum(sm.*M1))/sum(M1(:))*100; %river pixels that covered by the DEM over total water surface
    %
    if ratio==0
    flag1=0; % 0 no overlap; 1 overlapping 0< < 80%; 2 overlapping > 80%
    elseif ratio>=80
        flag1=2;
    else
        flag1=1;
    end
else
    flag1=0;
end

%Find the lowest stage DEM for this image.
%Gage must have a common reference strip for all images.
%For a reach of river, they may not have a common reference strip DEM.
if flag1==0|| (flagsect==1&&flag1~=2)
    %find the overlapping > 80% strips
    id=[];ratiog=zeros(n2strip,1);stageg=ones(n2strip,1); 
    %stagei, 1, stage (water level) of mono image is lower than strip.
    %        0, stage (water level) of mono image is higher than strip (good). 
    for k=1:n2strip
    fprintf(['\n \n Check strip ',num2str(k),':',f2is{k},'. '])
    sx=XYb2is{k}(:,1); sy=XYb2is{k}(:,2);
    idx=round((sx-rang0(1))/dx)+1;idy=round((sy-M.y(1))/(dy))+1;
    sm=poly2mask(idx,idy,nwy,nwx); % fast, apply to each polygon one by one.
    sm=double(sm);
    ratio=sum(sum(sm.*M1))/sum(M1(:))*100; %
    ratiog(k)=ratio;

    %check the relative stage of mono image and strip.
    datac(2)=struct('x',[],'y',[],'z',[]);
    datac(1)=datarsv2(k);XYbc(1)=XYb2is(k);
    datac(2).x=Msv(j).x;datac(2).y=Msv(j).y;datac(2).z=Msv(j).z;XYbc(2)=XYbis(j);
    [flagstatsj,idrefj]=lowest(datac,XYbc);
    stagei=0;
    if flagstatsj==1&&~isempty(idrefj)
        if idrefj(1)~=1
        %warning(['\n Image has lower stage than the DEM strip: ',fdiris{j},'/',fis{j},' VS ',fdir2is{k},'/',f2is{k}])
        stagei=1;
        end
    else %flagstatsj==0; %riverpic/riverwork2b/2/out1bp2 
        %fprintf(['\n Can not decide the relative water stage; Image and its DEM strip have no common water mask.',fdiris{j},'/',fis{j},' VS ',fdir2is{k},'/',f2is{k}])
        fprintf(['\n Can not decide the relative water stage; Image and its DEM strip have no common water mask.'])
        %Sometimes the reference just lose some part of water mask, the DEM data could still have good coverage.
        stagei=1;
    end
    stageg(k)=stagei;

    if ratio >=50 && stagei==0 %good coverage and water level of strip is lower than image.
        id=[id(:);k];
    end
    fprintf(['\n Strip ',f2is{k},' covers ',num2str(ratio),'%% of rivers in this mono image. '])
    fprintf(['\n Water level of this mono image is lower than the strip (1 yes (bad), 0 no): ',num2str(stagei)])
    end %k

if length(id)==0
   fprintf(['\n No strips cover at least 50%% of rivers in this mono image.\n'])
   continue
end
    
    %find the lowest;
    %It's possible that all selected strps have no common water area, hence can't find the lowest.
    %Because the water mask of strip does not cover the strip DEM domain exactly. Some missing masks.
    %[flagstatsj,idrefj]=lowest(datarsv2(id),XYb2is(id));
    % 90% coverage and lowest. If not, best coverage.
    id1=find(ratiog>=90&stageg==0);
    [flagstatsj,idrefj]=lowest(datarsv2(id1),XYb2is(id1));
    if ~isempty(idrefj)
      idrefj=idrefj(1);
      iref=id1(idrefj);
    end
    if flagstatsj==0
    warning('Lowest stage DEM strip not found (second try). Use the one with best coverage.')
    %continue
    [~,idrefj]=max(ratiog(id));
    iref=id(idrefj);
    end

%   iref=id(idrefj);

    k=iref;
    stripmetafile= [fdir2is{k},'/',f2is{k}];
    infile= strrep(stripmetafile,'meta.txt','dem.tif');
    fprintf (['\n The lowest stage DEM strip for this mono:',infile])
else %Use the reference for all
    iref=idref1;    
end

    elseif flaglowest==2 %use the given DEM; but to compare relative water level;
        iref=idlowest;
        fprintf (['\n Use the manually given DEM as lowest stage DEM:',f2is{iref},'\n'])
        fprintf (['\n Still compare the relative water level.\n'])
    end
    
% Check whether the water level of this image is higher than the reference.
% Skip if not.
if isempty(iref)
    warning('Lowest stage DEM strip not found.')
    continue
else
    datac(2)=struct('x',[],'y',[],'z',[]);
    datac(1)=datarsv2(iref);XYbc(1)=XYb2is(iref);
    datac(2).x=Msv(j).x;datac(2).y=Msv(j).y;datac(2).z=Msv(j).z;XYbc(2)=XYbis(j);
    [flagstatsj,idrefj]=lowest(datac,XYbc);
    %if idrefj~=1||flagstatsj==0 % []~=1 -> empty logical
    if flagstatsj==1&&~isempty(idrefj)
        if idrefj(1)~=1
        warning(['\n Image has lower stage than the reference DEM strip: ',fdiris{j},'/',fis{j},' VS ',fdir2is{iref},'/',f2is{iref}])
        continue
        end
    else %flagstatsj==0; %riverpic/riverwork2b/2/out1bp2 
        warning(['\n Can not decide the relative water stage; Image and its reference DEM strip have no common water mask.',fdiris{j},'/',fis{j},' VS ',fdir2is{iref},'/',f2is{iref}])
        %since the detected reference is the lowest stage strip DEM, still can use it at a risk.-> Don't do it.
        %Sometimes the reference just lose some part of water mask, the DEM data could still have good coverage.
        continue
    end
end
end %use the mannually selected DEM

%load lowest stage.
if isempty(iref)
     warning('Lowest stage DEM strip not found.')
elseif iref==idref1
%if iref==idref1
    data0=data0a;mt0=mt0a;
else
    k=iref;
    stripmetafile= [fdir2is{k},'/',f2is{k}];
    pr=dX4Sg2(k,:);     
    infile= strrep(stripmetafile,'meta.txt','dem.tif');
    data=readGeotiff(infile);
    %filtering out bad edges, e.g. 20110804
%     Me=~(imdilate((data.z==-9999),ones(round(30*8)))); 
%data.z(~Me)=-9999; %April 2019, pzxy may change the value of -9999.
    dx=data.x(2)-data.x(1); dy=data.y(2)-data.y(1);[nwy,nwx]=size(data.z);
    sx=XYb2is{k}(:,1); sy=XYb2is{k}(:,2); 
    idx=round((sx-data.x(1))/dx)+1;idy=round((sy-data.y(1))/(dy))+1;
    Med=~poly2mask(idx,idy,nwy,nwx); % fast, apply to each polygon one by one.
    Me=~(imdilate(Med,ones(round(30*8)))); 
    data.z(~Me|data.z==-9999)=NaN;
    
    mtFile = strrep(infile,'dem.tif','matchtag.tif');
    mt=readGeotiff(mtFile);
    % Align the DEM coordinates to reference DEM using the translation parameters from coregistration.
    data.x=data.x- pr(2);data.y=data.y- pr(3);data.z=data.z- pr(1);

    data0=data;mt0=data;mt0.z=mt.z;
end
data0.z(data0.z==-9999)=NaN;
% % % End of finding reference DEM.

%get the shores.
if 0
Medgs=(M.z==-1);%isnan(Mstrip.z(:,:));
Med=imdilate(Medgs,ones(4));
Medgs=Med;
Modj=M.z;Modj(Medgs)=0;
%Modj= bwareaopen(Modj, 1000*500); %remove small clusters
%Modfil = bwareaopen(~Modj, 1000*5); %fill small areas: 1e4*4m^2
Modj= bwareaopen(Modj, lakearea/resr/resr); %remove small clusters
Modfil = bwareaopen(~Modj, cloudarea/resr/resr); %fill small areas: 1e4*4m^2
Modfil=~Modfil;

Md1 = imdilate(Modfil, ones(3));
M=logical(Md1-Modfil);
M=M&~Medgs; 
else
%M=M.coast;%M.z;
[M,~]=mask2coast(M);
end

ymd=str2num(fis{j}(6:13));

% %second method : imagery-altimetry method
% Retrieve the shoreline location in the lowest-stage DEM coordinates system.
y0s=0; %630998;%y0;
xriv=X(M); yriv=Y(M); %data.x- p(2), in reference DEM coordinates
% Interpolate the shoreline heights from the lowest-stage DEM.
zriv= interp2(data0.x,data0.y,double(data0.z),xriv,yriv,'*linear',NaN);
dempmt= interp2(mt0.x,mt0.y,double(mt0.z),xriv,yriv,'*nearest',0);
epochorg=Y(M);
zriv(isnan(zriv))=-9999;

%Detect the outliers and mean and std
try
[idkpb]=outlier(odir,zriv,dempmt,epochorg,ymd);
catch e
     fprintf('\nThere was an error! The message was:\n%s',e.message);
     fprintf(['\nstripmetafile is ',stripmetafile])
save testoutlier.mat odir zriv dempmt epochorg ymd  -v7.3
continue
end

yp=epochorg(idkpb);

[LAT,LON]=polarstereo_inv(X(M),Y(M),[],[],70,-45);
%Write the shoreline heights
output=[LAT(idkpb),LON(idkpb),double(zriv(idkpb)),(yp(:)-y0s)];
% ofile=['rivprof',num2str(ymd),'b.dat']; %remove outlier
%ofile=[odir,'/rivprof',num2str(ymd),satname,'bj',num2str(j),'.dat'];
ofile=[odir,'/rivprof',texttar,'bj',num2str(j),'.dat'];
save(ofile,'output','-ascii') 

%plot the orthoimage and shorelines
infile=[fdiris{j},'/',fis{j}];
orFile = strrep(infile,'.xml','.tif');
if ~exist(orFile,'file')
orFile = strrep(infile,'meta.txt','ortho.tif');
end

if ~exist(orFile,'file')
    warning([orFile,'not exist for plotting'])
end

or=readGeotiff(orFile);
resrc=40;
res=or.info.map_info.dx;
nsr=resrc/res; dsr=res/resrc;

[~,~,nb]=size(or.z);
if nb==4
    iNIR1=4;
    or.z=or.z(:,:,iNIR1);
elseif nb==8
    iNIR1=7;
    or.z=or.z(:,:,iNIR1);
end

or.z = imresize(or.z,dsr);
or.x = imresize(or.x,dsr);
or.y = imresize(or.y,dsr);
[X,Y]=meshgrid(or.x,or.y);
[LATa,LONa]=polarstereo_inv(X- p(2),Y- p(3),[],[],70,-45);

figure;
set(gcf,'Color','white')
set(gca,'FontSize', 12);
set(gcf, 'PaperPosition', [0.25 2.5 6 2.5]);
hold all;
surf(LONa,LATa,or.z); % too slow
shading interp;
colorbar;colormap gray;view(0,90)
% caxis([0 10000])
hold on;plot3(LON(idkpb),LAT(idkpb),1e9*ones(size(LAT(idkpb))),'r.','Markersize',8)
 view(0,90)
box on
hl=xlabel('Longitude ($^{\circ}$)');
set(hl, 'Interpreter', 'latex');
hl=ylabel('Latitude ($^{\circ}$)');
set(hl, 'Interpreter', 'latex');
% ofile=['rivloc',num2str(ymd)]; %mean and std
%ofile=[odir,'/rivloc',num2str(ymd),'j',num2str(j)];
ofile=[odir,'/rivloc',texttar,'j',num2str(j)];
% hold on;plot( - (147+50/60+20/3600), 64+47/60+34/3600, 'bs','Markersize',12) %gauge
% hold on;plot(riv(:,1),riv(:,2),'b-','linewidth',4) %airport
% hold on;plot3(LON(idkp(in)),LAT(idkp(in)),1e9*ones(size(LAT(idkp(in)))),'b>','Markersize',8)
% hold on;plot3(LON(idkpb(inb)),LAT(idkpb(inb)),1e9*ones(size(LAT(idkpb(inb)))),'co','Markersize',8)
title(texttar)
saveas(gcf,ofile,'fig')

end %for j

close all

return
end
