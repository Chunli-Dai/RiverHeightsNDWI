addpath(genpath('/home/dai.56/arcticdemapp/river/rivergithub2/'))
constant

if 0
load('mat2.mat')
load dX4Sg.mat %hi

% % %
% Get the profiles of images at gage.
fprintf ('\n Step 2: Get the profiles of images at gage.')
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

dX4Sg

[Co]=riverprofsub(odir,datarsv,XYbis,fis,fdiris,dX4Sg,XYb2is,f2is,fdir2is,dX4Sg2);
end %if 0

loneq=-148.8177778;lateq= 69.0158333;
%get river centerline
if 0
load clsv2.mat %use the saved one %hi
else
%Merge water masks
load('1/gage1bp/wprob2.mat')
data1.x=xout;data1.y=yout;data1.z=prob;
load('2/gage2bp/wprob2.mat')
data2.x=xout;data2.y=yout;data2.z=prob;
load('3/gage3bp/wprob2.mat')
data3.x=xout;data3.y=yout;data3.z=prob;
load('4/gage4bp/wprob2.mat')
data4.x=xout;data4.y=yout;data4.z=prob;

x=[data1.x(:);data2.x(:);data3.x(:);data4.x(:);];y=[data1.y(:);data2.y(:);data3.y(:);data4.y(:);];
rang0=[min(x) max(x)  min(y) max(y)];
ranget=round(rang0/resr)*resr;rang0=ranget;
tx=ranget(1):resr:ranget(2);ty=ranget(4):-resr:ranget(3);
xout=tx;yout=ty;
nsuby=length(yout);nsubx=length(xout);

tz1 = interp2(data1.x,data1.y,data1.z,xout,yout','*nearest',255);
tz2 = interp2(data2.x,data2.y,data2.z,xout,yout','*nearest',255);
tz3 = interp2(data3.x,data3.y,data3.z,xout,yout','*nearest',255);
tz4 = interp2(data4.x,data4.y,data4.z,xout,yout','*nearest',255);
tall=tz1;tall(:,:,2)=tz2;tall(:,:,3)=tz3;tall(:,:,4)=tz4;
prob=min(tall,[],3);

iprobthre=probthre;
jump=-1*ones(nsuby,nsubx,'int8'); %%-1,non value;1 water; 0 non-water
jump(prob>=iprobthre&prob~=255)=1;
jump(prob<iprobthre)=0;
%jumpsv=jump;
save jumpc.mat xout yout jump

%Plot water mask of 50% probability.
[X,Y]=meshgrid(xout,yout);
xmin=min(X(jump==1));xmax=max(X(jump==1));
ymin=min(Y(jump==1));ymax=max(Y(jump==1));
idx=xout>=xmin&xout<=xmax;idy=yout>=ymin&yout<=ymax;
clear data
data.x=xout(idx);data.y=yout(idy);tz=jump(idy,idx);

%save M water mask
projstr='polar stereo north';
ofile=['watermask50.tif'];
writeGeotiff(ofile,data.x,data.y,uint8(tz),1,255,projstr) %wrong, M is just shoreline, not water mask.

%Get river centerline from water mask
res1=resr;
tz(tz==-1)=0;
Modj= bwareaopen(tz, round(lakearea/res1/res1)); %remove small clusters
Modfil = bwareaopen(~Modj, round(cloudarea/res1/res1)); %fill small areas: 1e4*4m^2
Modfil=~Modfil;
data.z=Modfil;

npt=sum(sum(data.z==1));
if npt>0
	try
        [c]=mask2centerline(data);
        catch e
                fprintf('There was an error! The message was:\n%s',e.message);
                fprintf('\n')
        end

       %adjust the centerline to go uphill for ProcessTananaFairbanks.m.
        width=500; %m 
        pt=[c.X(1) c.Y(1);c.X(end) c.Y(end)];
        %use data0r to get the height
        if 0
        c.X=flip(c.X);
        c.Y=flip(c.Y);
        end
	save clsv2.mat c
else
   fprintf(['No water pixels in the water mask. ']); 
end
end

odir='./'
M=isnan(c.X)|isnan(c.Y);
c.X(M)=[];c.Y(M)=[];

% Requires the input centerline to go uphill.
Co=ProcessTananaFairbanks(odir,c,lateq,loneq);

%clear datarsv data

fprintf ('\n Step 4: Height time series analysis at the gage.')
gageheights
