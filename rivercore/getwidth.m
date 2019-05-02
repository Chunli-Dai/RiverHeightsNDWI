
function [gagewidth,widthp]=getwidth(data,infile,c,lateq,loneq);
%function [gagewidth]=getwidth(infile,lateq,loneq)
% Given water mask to get width profile
% %Better use a fixed rivercenterline for consistent width comparison.
% input: data, water mask
%	 c, river centerline, c.X, c.Y 

% loneq=-148.8177778;lateq= 69.0158333; %~/data/chunli/scripts/stations2.gmt
%xeq,yeq]=polarstereo_fwd(lateq,loneq,[],[],70,-45);

constant

lakearea=100*500*4/4;
wmin=20; %20 m%mininum river width;

gagewidth=[];
widthp.x=[]; widthp.y=[];

% infile='watermaskWV02_20160827bj80.tif';
% infile='watermaskWV01_20160713bj35.tif';
% infile='watermaskWV01_20110720bj24.tif';
% infile='watermaskWV01_20080608bj15.tif'
% infile='watermaskWV01_20080514bj14.tif'
%data=readGeotiff(infile); %for getting river width, dx dy better be equal ine
%data, water mask in 2 m resolution, translation parameters applied.
[demdir,filename,ext] =fileparts([(infile)]);

%Interpolate water mask to be on integer grids. If not, the round operator
%(centerline's polar stereographic coordinate to image coordinate) slightly
%alter the centerline, which can accumulate errors in distance up to 310m.
resr=2;
ranget=[[min(data.x) max(data.x) min(data.y) max(data.y)]/resr];
ranget=[ceil(ranget(1)) floor(ranget(2)) ceil(ranget(3)) floor(ranget(4))]*resr;
tx=ranget(1):resr:ranget(2);ty=ranget(4):-resr:ranget(3);
tz = interp2(data.x,data.y,data.z,tx,ty','*nearest',0); %1  water, 0 nonwater;
datar= struct();
datar.x=tx;datar.y=ty;datar.z=tz; 
data=datar;

resx=mean(data.x(2:end)-data.x(1:end-1));resy=mean(data.y(2:end)-data.y(1:end-1));
resr=mean([abs(resx),abs(resy)]);

% get a rough boundary for river path
if 0 %old way as in /home/dai.56/chunli/scripts/maskentropy.m
cf = 0.1;%0.5; %boundary curvature factor (0= point boundary, 1 =conv hull)
MJr = imresize(data.z,0.1);
[idy,idx]=find(MJr==1);
k = boundary(idx,idy,cf);
[nyov,nxov]=size(MJr);
Mb = poly2mask(idx(k),idy(k),nyov,nxov);
MJbd = imresize(Mb,size(data.z),'nearest');
else %new way
    datar.x= imresize(data.x,0.1);datar.y= imresize(data.y,0.1);
    datar.z= imresize(data.z,0.1);
    [X,Y]=meshgrid(datar.x,datar.y);
    idx=X(datar.z==1);idy=Y(datar.z==1);
    cf = 0.1;%0.5; %boundary curvature factor (0= point boundary, 1 =conv hull)
    k = boundary(idx,idy,cf);
% hold on;plot(idx(k),idy(k),'r-')
    bdx=idx(k);bdy=idy(k);
end


if 1 %filter out small clusters. fill holes. (Filled holes may increase the width)
Mstrip=data;
Medgs=(Mstrip.z==-1);%isnan(Mstrip.z(:,:));
Med=imdilate(Medgs,ones(4));
Medgs=Med;
Modj=Mstrip.z;Modj(Medgs)=0;
% resx=mean(Mstrip.x(2:end)-Mstrip.x(1:end-1));resy=mean(Mstrip.y(2:end)-Mstrip.y(1:end-1));resr=mean([abs(resx),abs(resy)]);
Modj= bwareaopen(Modj, round(lakearea/resr/resr)); %remove small clusters
% Modfil = bwareaopen(~Modj, round(cloudarea/resr/resr)); %fill small areas: 1e4*4m^2
% Modfil=~Modfil;
Modfil=Modj;
data.z=Modfil;
end

% [c]=mask2centerline(data,1); %data better be rivers only, with no lakes. 

%Use saved centerline
if 0
if exist('clsv2.mat','file') %better use a fixed rivercenterline for consistent width comparison.
load clsv2.mat
else
    fprintf(['River centerline file clsv2.mat not found for width calculation. ']);
    return
end
end

% figure;plot(c.X,c.Y)
[clx,cly]=polarstereo_fwd(c.Y,c.X,[], [],70,-45);
% if centerline resolution is too large, eg. 40 m, interpolate it to be 2m. Must do this!
% If the input centerline (green) has too large interval (40 m spacing), the smoothed line is too far from river in width_from_mask.
S = [0; cumsum(sqrt(diff(clx(:)).^2+diff(cly(:)).^2))];
rescl=nanmean(S(2:end)-S(1:end-1));
if rescl > 2*resr %
fprintf(['\n Densify the centerline for width calculation!'])
lat=cly;lon=clx;maxdiff=resr;
[latout,lonout] = interpm(lat,lon,maxdiff);
clx=lonout;cly=latout;
end
S = [0; cumsum(sqrt(diff(clx(:)).^2+diff(cly(:)).^2))];
rescl=nanmean(S(2:end)-S(1:end-1));%average node distances of centerline.

%polar stereographic coordinates to image coordinates.
[ny,nx]=size(data.z);
clear cl
if 0 % old %  %round cause truncation error, which can be accumulated to distance error.
cl(:,1)=round((clx-data.x(1))/resr)+1; 
cl(:,2)=round(-(cly-data.y(1))/resr)+1;
else
cl(:,1)=((clx-data.x(1))/resr)+1; 
cl(:,2)=(-(cly-data.y(1))/resr)+1;
end
% M=cl(:,1)>=1&cl(:,1)<=nx&cl(:,2)>=1&cl(:,2)<=ny;
% cl(~M,:)=[];
M = inpolygon(clx,cly,bdx,bdy); %Use only the centerline within the water mask.
id1=find(M==1);ks=id1(1);ke=id1(end);
M(ks:ke)=1; %avoid centerline being separate into pieces; %eg.watermaskWV01_20080514bj14.tif
cl(~M,:)=[];
lencl=length(cl(:,1));
% figure;plot(cl(:,1),cl(:,2),'>-')

%the coordinate of first point in cropped centerline along input
%centerline.
id1=find(M==1);k=id1(1);
s0=S(k);

% spacing = Wn/2; %spacing of centerline, in pixels.
spnode=100; % spacing in meter
nnode=round(lencl*rescl/spnode); 
if 1% nnode <=400 %17km*2 river length / 100m
spacing = max([round(spnode/rescl),1]); %100 m; spacing of centerline, in pixels (nodes of centerlines).
else
fprintf('\nLess number of nodes for faster computation of river width! \n')
nnode2=100; %less number of nodes for faster computation.
spacing=round(lencl/nnode2);
fprintf(['\nThe spacing of centerline is ',num2str(spacing*rescl),' instead of ',num2str(spnode),'m. \n'])
end

try
[Wm, SWm] = width_from_mask(data.z, cl, spacing);
catch e
     fprintf('width_from_mask error! The message was:\n%s',e.message);
     gagewidth=0;return
end
% figure;imagesc(I,'alphadata',Ibuffer);view(0,-90)

widthp.x=SWm*resr+s0; %centerline distance
widthp.y=Wm*resr; %width
fprintf(['\n Width nodes number: ',num2str(length(widthp.x)),'; estimate:',num2str(nnode),'.'])

ofile=strrep(infile,'watermask','widprof'); %watermaskWV02_20150714bj73.tif widprofWV02_20150714bj73.dat
ofile=strrep(ofile,ext,'.dat');
output=[widthp.x(:), widthp.y(:)]; %distance along centerline (m), width(m)
save(ofile,'output','-ascii')

figure
plot((SWm*resr+s0)*1e-3,Wm*resr,'r');
xlabel('streamwise distance (km)'); ylabel('width (m)')
legend('W_m_a_s_k')

%To get the width of gage
%gage coordinate to image coordinate
if exist('lateq','var')
[xeq,yeq]=polarstereo_fwd(lateq,loneq,[],[],70,-45);
gagex=((xeq-data.x(1))/resr)+1;
gagey=(-(yeq-data.y(1))/resr)+1;
M = inpolygon(xeq,yeq,bdx,bdy); %Check if gage is within the water mask boundary.
%gage coordinate to centerline coordinate
dist=sqrt( (gagey-cl(:,2)).^2 + (gagex - cl(:,1)).^2);
[distmin,k]=min(dist); % if gage is out of the image domain, k is the 1st point of centerline.
dists = sqrt(diff(cl(:,1)).^2+diff(cl(:,2)).^2); % distances between each centerline node
cumdists = [0; cumsum(dists)];
xobs=cumdists(k); %pixels along cropped centerline; %could be end points (rather than closest point) , if the gage is outside of the boundary.

%method2
dist=sqrt( (yeq-cly(:)).^2 + (xeq - clx(:)).^2);
[distmin,k]=min(dist);
xobs2=S(k);
xobs1=xobs*resr+s0;
df=xobs2-(xobs*resr+s0);
if M % gage within image boundary
fprintf(['\n Gage coordinate along input centerline, xobs1:',num2str(xobs1),', xobs2:',num2str(xobs2),', difference:',num2str(df),''])
else
fprintf(['\n Gage is out of the water mask boundary.'])
end

%1km reach along SWm
% gagereach=1e3/resr;%1KM, in pixels.
% x2=SWm;
gagereach=1e3;%1KM, 
x2=widthp.x;
% id=find(abs(x2-xobs)<=gagereach/2.);
% id=find(abs(x2-xobs)<=gagereach/2.&Wm*resr>=wmin);
id=find(abs(x2-xobs2)<=gagereach/2.&Wm*resr>=wmin);
nn=length(id);
hest=Wm(id);
xsrdx=mean(x2(2:end)-x2(1:end-1));
lenr=nn*abs(xsrdx);
gagewidth=nanmean(hest)*resr;
% lenr=max(x2(id))-min(x2(id))+abs(xsrdx);
if (lenr-gagereach)< -gagereach*0.1
    fprintf(['\n This profile does not have a good coverage of the gage within ',num2str(gagereach),' m:',infile]);
    gagewidth=0;
end

hold on;plot((SWm(id)*resr+s0)*1e-3,Wm(id)*resr,'g>');
% hold on;plot(([xobs*resr,xobs*resr]+s0)*1e-3,[0 300],'b-')
hold on;plot(([xobs2,xobs2])*1e-3,[0 300],'b-')

%title(filename(1:13))
title(filename(10:22))
ofile=['pic/',filename,'width'];
legend(filename(10:22))
% print('-dpdf','-r300',ofile) 
saveas(gcf,ofile,'fig')


%plot the profile along input river centerline.
%first point of the centerline of this image

figure (1)
hold all
h=plot((SWm*resr+s0)*1e-3,Wm*resr,'-','DisplayName',filename(10:22));
hold on;plot((SWm(id)*resr+s0)*1e-3,Wm(id)*resr,'g>');
% hold on;plot(([xobs*resr,xobs*resr]+s0)*1e-3,[0 300],'b-')
hold on;plot(([xobs2,xobs2])*1e-3,[0 300],'b-')
xlabel('streamwise distance (km)'); ylabel('width (m)')
% legend(h,filename(10:22))
% legend(gca,'off');    
legend('show');
ofile=['pic/allwidth'];
saveas(gcf,ofile,'fig')

end

close all
end


