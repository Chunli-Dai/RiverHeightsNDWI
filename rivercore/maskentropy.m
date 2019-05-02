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

% read ortho
if ~isempty(rangeov)
or=readGeotiff(orFile,'map_subset', rangeov);
else
or=readGeotiff(orFile); % to keep long profiles
end
res=or.info.map_info.dx;
resr=2.;dsr=res/resr;
dsr2=0.25;resr2=res/dsr2; %for faster computation 
dsr2=1;resr2=res/dsr2; %for better resolution

or.z = imresize(or.z,dsr);
or.x = imresize(or.x,dsr); 
or.y = imresize(or.y,dsr);

mt.x=or.x;mt.y=or.y; mt.z=ones(size(or.z)); P=ones(size(or.z));

orsv=or;
or = or.z; % get rid of map data

% resize ortho to 8m  -> To check difference between 8m and 2 m
or = imresize(or,dsr2);
% or=uint16(or); %chunli 2017/6/7
or=double(or); %chunli 11/21/2018

% subtraction image
if 0 %Don't use entropy;  save computation time
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

if flagplot==1
figure;
hold all
set(gcf,'Color','white');set(gca,'FontSize', 18);set(gcf, 'PaperPosition', [0 0 6 4]);set(gcf, 'PaperSize', [ 6 4]);
imagesc(mt.x*1e-3,mt.y*1e-3,double(J))
colormap jet;colorbar
title(['Entropy Pan-band ',name(6:13)]);
axis equal
ofile=[name,'J'];
print('-dpng','-r300',ofile)
saveas(gcf,ofile,'fig')
end

end % if 0

% Get the reflectance 
data.z=or;
try
  or=DN2reflpan(data,ifile);
catch e
     fprintf('There was an error! The message was:\n%s',e.message);
end


if 0 %frozen==7 %Not work well
%     &Mrho;
M=(J<0.5)&~(imdilate(or <= 0,ones(round(9*2*8/resr))));%running dark river; 20111008
clear J
% M = bwareaopen(M, 1000*5);
% M = bwareaopen(M, 1000*200);%20121011
meanor=nanmean(double(or(M)));stdor=nanstd(double(or(M)));
% Mor=(or>meanor-stdor)&(or<meanor+stdor); %running dark river; 20130526
end % if frozen

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
meanor=0.3;stdor=0.3; %arbituray value
Mor=(or<meanor+stdor); %20121011
warning('Process Panchromatic image without a priori water mask!')
% mean1o=meanor;std1o=stdor;

elseif step ==2   %step 2, use the refined threshold.
%retrieving statistics of NDWI over a priori water mask.
% wm: 1 water, 0 land, -1 void
%resrc=abs(wm.x(2)-wm.x(1));%old
resx=mean(mt.x(2:end)-mt.x(1:end-1));resy=mean(mt.y(2:end)-mt.y(1:end-1));
resrc=mean([abs(resx),abs(resy)]);
% Md1 = imdilate(wm.z, ones(widthstat/resrc)); %make sure water is absolutely water.
Md1=wm.z==1;
%Remove small lakes, which may have different reflectance than river water
%Md1= bwareaopen(wm.z, round(lakearea/resrc/resrc)); %remove small clusters

wm1 = interp2(wm.x,wm.y,Md1,mt.x,mt.y','*nearest',0); %out of region: 0, non water.
wm1sv=wm1;
ndwil=double(or(wm1==1)); %water
mean1=nanmean(ndwil);std1=nanstd(ndwil);cnt1=sum(~isnan(ndwil));
b=ndwil(~isnan(ndwil));
if 0% length(b) >=cntmin 
ithres=histhres(name,b); %not working well
mean1=mean(ithres);std1=abs(ithres(2)-ithres(1))/2.;
end
mean1o=mean1;std1o=std1;

if 0
Md1 = imdilate(wm.z~=0, ones(600/resrc)); %make sure land is absolutely land. 300 m width expansion
%Me1=imerode(wm.z==0, ones(widthstat/resrc)); %land
wm1 = interp2(wm.x,wm.y,Md1,mt.x,mt.y','*nearest',1); %out of region: 1,  water.
else
%widthriv=600; %to test
widthriv=100; %to test
if isfield(wm, 'w50')
Md1=wm.w50~=0;
else
Md1=wm.z~=0; %first interpolate; then dilate to save computation time.
end
wm1 = interp2(wm.x,wm.y,Md1,mt.x,mt.y','*nearest',1); %out of region: 1,  water.
%Only use a 300m width land on both sides of the river
Med=imdilate(isnan(or),ones(round(widthriv*3/resrc))); %image edge
Md1 = imdilate(wm1~=0, ones(round(widthriv/resrc))); %make sure land is absolutely land. 300 m width expansion
Mdd = imdilate(Md1, ones(round(widthriv/resrc))); 
wm1=~((Mdd-Md1)& ~Med);
end

ndwil=double(or(wm1==0)); %land
a=ndwil(~isnan(ndwil));
mean2=nanmean(ndwil);std2=nanstd(ndwil);cnt2=sum(~isnan(ndwil));
%clear wm wm1

if 0 %flagplot==1
%Plot the PDF of land (a) and water (b)
figure
set(gcf,'Color','white')
set(gca,'FontSize', 12);
set(gcf, 'PaperPosition', [0 0 6 3]); 
set(gcf, 'PaperSize', [ 6 3]); 
hold all;
hold on
ha=histogram(a,'Normalization','pdf','Facecolor','red','Edgecolor','red'); %,'Edgecolor','green'
hb=histogram(b,'Normalization','pdf','Facecolor','blue','Edgecolor','blue');
ofile=[name,'pdflw'];
saveas(gcf,ofile,'fig')
end

dmthresp=0.02; istdthresp=0.02;
ithreshold=0;
badflag=0;
meanor=mean1;stdor=std1;
ithres=[0 0];sigma0hat=[0.5,0.5];%default
if (isnan(mean1)||cnt1<cntmin) 
    badflag=1;
%elseif (std1<istdthresp && abs(mean1-mean2)>dmthresp)
else
[ithres,badflag,sigma0hat]=adaptivethres(name,a,b);
mean1=mean(ithres);std1=abs(ithres(2)-ithres(1))/2.;

meanor=mean1;stdor=std1;
Mor=(or>meanor-stdor)&(or<meanor+stdor); %20121011

if 0
%automated searching of threshold.
[mtli]=autothres(or,meanor,stdor)
Mor=(or>meanor-mtli*stdor)&(or<meanor+mtli*stdor); %20121011
end % if 0

%else %discard the image.
%    badflag=1; 
end %if isnan(mean1)..


fprintf('\n %f %f %f %f %f %f %d %f %f (NDWI mean std thres) %s \n',mean1o,std1o,mean2,std2,ithres,badflag,sigma0hat,ifile);
if 0
fid = fopen('./panstats.dat','a'); %statistics
fprintf(fid,' %f %f %f %f %f %f %d %f %f (NDWI mean std thres) %s \n',mean1o,std1o,mean2,std2,ithres,badflag,sigma0hat,ifile);
fclose(fid);
end

end  % if step

if flagplot==1
Mstrip=m;
[X,Y]=meshgrid(Mstrip.x,Mstrip.y);

if 0
figure;
hold all
set(gcf,'Color','white');set(gca,'FontSize', 18);set(gcf, 'PaperPosition', [0 0 6 4]);set(gcf, 'PaperSize', [ 6 4]);
imagesc(Mstrip.x*1e-3,Mstrip.y*1e-3,double(Mstrip.n))
colormap gray;colorbar
hold on;plot(X(wm1sv==1&~isnan(or))*1e-3,Y(wm1sv==1&~isnan(or))*1e-3,'r.')
title(['Ortho Pan-band ',name(6:13)]);
axis equal
ofile=[name,'apriori'];
axis([-2232.7 -2231 548 549.5])
%print('-dpng','-r400',ofile)
saveas(gcf,ofile,'fig')
end

figure;
hold all
set(gcf,'Color','white');set(gca,'FontSize', 18);set(gcf, 'PaperPosition', [0 0 6 4]);set(gcf, 'PaperSize', [ 6 4]);
imagesc(Mstrip.x*1e-3,Mstrip.y*1e-3,double(Mstrip.n))
colormap gray;colorbar
%hold on;plot(X(M),Y(M),'ro')
title(['Ortho Pan-band ',name(6:13)]);
axis equal
ofile=name;
try
caxis([meanor-3*stdor, meanor+3*stdor])
end
%axis([-2232.7 -2231 548 549.5])
try 
print('-dpng','-r400',ofile)
saveas(gcf,ofile,'fig')
catch e
     fprintf('There was an error! The message was:\n%s',e.message);
	try 
	   print('-dpdf','-r400',ofile)
	end
end
%saveas(gcf,[name],'fig') %too big
end

if badflag==1    
Mstrip=struct(); Mstrip.x=[];Mstrip.y=[];Mstrip.z=[];Mstrip.coast=[];Mstrip.n=[];
m=Mstrip;
close all
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

%apply a priori river buffer zone.
if isfield(wm, 'buf')
buf= interp2(wm.x,wm.y,wm.buf,m.x,m.y','*nearest',0); %out of region: 0, void
m.z(~buf)=-1; 
%Do not remove clusters for fair comparison of river stages.
%Since Removing small clusters can remove different lakes from different masks, depending on threshold.
%Decide To remove the clusters, 1) for better calculation of river width; 
%2) only small lakes (associate with lower stage) will be removed, so the lowest stage will still be corrected selected.
%Update: do not apply the buffer until last step; since it will affect the water probability map.
end

% m.info = mt.info;
% m.Tinfo = mt.Tinfo;

Mstrip=m;
%[M,~]=mask2coast(Mstrip);
[M,Modfil]=mask2coast(Mstrip);
if isfield(wm, 'buf') %same control; do not apply the removal of lakes until last step.
m.z=int8(Modfil); %remove clusters and lakes, for better calculation of river width out of water mask.
end

if flagplot==1
hold on;plot(X(M)*1e-3,Y(M)*1e-3,'ro')
title(['Ortho Pan-band ',name(6:13)]);
axis equal
%saveas(gcf,[name],'fig') %too big
%print('-dpdf','-r400',name)
%axis([-2232.7 -2231 548 549.5])
try 
print('-dpng','-r400',ofile)
saveas(gcf,ofile,'fig')
catch e
     fprintf('There was an error! The message was:\n%s',e.message);
	try 
	   print('-dpdf','-r400',ofile)
	end
end
end

m.coast=M;%[X(M),Y(M)]; %merged coast for the strip

close all
return
end

 
