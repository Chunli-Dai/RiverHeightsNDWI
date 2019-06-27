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

figure;imagesc(data.x*1e-3,data.y*1e-3,data.z);colorbar;view(0,-90)

resx=mean(data.x(2:end)-data.x(1:end-1));resy=mean(data.y(2:end)-data.y(1:end-1));
resr=mean([abs(resx),abs(resy)]);
res1=resr;

%preparation: remove small clusters 
BW=data.z;
Medge=(BW==-1);
BW(BW==-1)=0;
Modj= bwareaopen(BW, round(lakearea/2/resr/resr/4)); %remove clusters of lakes 500m by 500m
Modj=int8(Modj);
Modj(Medge)=-1;
data.z=Modj;

figure;imagesc(data.x*1e-3,data.y*1e-3,data.z);colorbar;view(0,-90)
title('Remove small clusters')

%reduce resolution; purposes: connect segments broken by bridges. 
% 2m resolution too slow, 8 hours.
%40m resolution -> river can be too narrow, so a 40 m resoltuion may break the rivers into sections.
if 0
resrc=10;%10m % faster and better than 2m, e.g. remove small thin tributary and make the centerline along main stream.
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
end

mp=roipoly;%mannually draw polygon;
data.z(~mp)=0;

%crop out edges
data2=cropmatrix(data,data.z);

[c]=mask2centerline(data2);

[clx,cly]=polarstereo_fwd(c.Y,c.X,[], [],70,-45);

hold on;plot(clx*1e-3,cly*1e-3,'r-','linewidth',2)

%get river width
tic
infile='watermaskstacked';
[gagewidth,widthp]=getwidth(data2,infile,c);
%get buffer width: 80 percentile of river width times 2.
width80 = prctile(widthp.y(widthp.y>0),80);	
fprintf(['\n80 percentile river width is:',num2str(width80),'m.']);
if width80<10||isnan(width80)||isinf(width80)
   width80=200;
   fprintf(['\n Unable to get a reasonable width, use the fixed width:',num2str(width80),'m.']);
end
toc

c.widave=mean(widthp.y(widthp.y>0));

save clsv2.mat c

