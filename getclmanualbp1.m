%Manually select the river centerline from water probability map.

addpath(genpath(['/Users/chunlidai/surge/home/dai.56/arcticdemapp/river/rivergithub2v2/']))

constant

load wprob2.mat
figure;imagesc(xout*1e-3,yout*1e-3,prob);colorbar;title('wprob2')
caxis([0 100]);colormap jet

%get the 50% mask.
iprobthre=probthre;
[nsuby,nsubx]=size(prob);
jump=-1*ones(nsuby,nsubx,'int8'); %%-1,non value;1 water; 0 non-water
jump(prob>=iprobthre&prob~=255)=1;
jump(prob<iprobthre)=0;

%crop out edges
data.x=xout;data.y=yout;data.z=jump;
data=cropmatrix(data,data.z);

figure;imagesc(data.x*1e-3,data.y*1e-3,data.z);colorbar

resx=mean(data.x(2:end)-data.x(1:end-1));resy=mean(data.y(2:end)-data.y(1:end-1));
resr=mean([abs(resx),abs(resy)]);
res1=resr;

%preparation: remove small clusters for the downsizing, so that wont
%connect smaller clusters to the main rivers.
BW=data.z;
Medge=(BW==-1);
BW(BW==-1)=0;
Modj= bwareaopen(BW, round(lakearea/2/resr/resr/10)); %remove clusters of lakes 500m by 500m
Modj=int8(Modj);
Modj(Medge)=-1;
data.z=Modj;

%reduce resolution; purposes: connect segments broken by bridges. 
% 2m resolution too slow, 8 hours.
%40m resolution -> river can be too narrow, so a 40 m resoltuion may break the rivers into sections.
resrc=30;%10m % faster and better than 2m, e.g. remove small thin tributary and make the centerline along main stream.
clear datar
datar.x=data.x(1):resrc:data.x(end);datar.y=data.y(1):-resrc:data.y(end);
tz = interp2(data.x,data.y,data.z,datar.x,datar.y','*nearest',0);
datar.z=tz;
res1=resrc;resr=res1;

%preparation: remove small clusters
BW=datar.z;
Medge=(BW==-1);
BW(BW==-1)=0;
Modj= bwareaopen(BW, round(lakearea/2/resr/resr)); %remove clusters of lakes 500m by 500m
Modj=int8(Modj);
Modj(Medge)=-1;
datar.z=Modj;

wm1=mask2river(datar);
wm1.z(wm1.z==-1)=0;

mp=roipoly;%mannually draw polygon;
wm1.z(mp)=0;

[c]=mask2centerline(wm1);

[clx,cly]=polarstereo_fwd(c.Y,c.X,[], [],70,-45);

hold on;plot(clx*1e-3,cly*1e-3,'r-','linewidth',2)


