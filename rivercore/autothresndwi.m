function [mtli,badflag2,sigma0hat]=autothresndwi(datax,datay,or,meanw,meanl)
%automatically searching threshold to find the one that gives the minimum change of water area.
%Wang et al. (2014, https://doi.org/10.1016/j.rse.2014.06.004).
constant
%flagplot=1;

NDWI=or;
edge=isnan(or);
resrc=abs(datax(2)-datax(1));

filename='t1'

%automated searching of threshold.
mtl=min([meanw,meanl]):0.01:max([meanw,meanl]); nc=length(mtl);
waterarea=zeros(size(mtl));
Topt=zeros(size(mtl));

fprintf(['hi 4 ',num2str(int32(clock)),'.\n'])
whos

%figure(1)
for i=nc %1:nc
fprintf(['Calculating ',num2str(i),' out of ',num2str(nc),' candidates, at ',num2str(int32(clock)),'.\n'])
whos
mtli=mtl(i);
Mor=(or>mtli); %20121011 ; slow
waterarea(i)=sum(Mor(:));

%get the water mask for optimal thresholding method.
Md1=Mor==1; %1 water
wm1=Md1;

ndwil=NDWI(wm1==1); %water
mean1=nanmean(ndwil);std1=nanstd(ndwil);cnt1=sum(~isnan(ndwil));
b=ndwil(~isnan(ndwil)); %water
if flagplot==1
figure;
hold all
set(gcf,'Color','white');set(gca,'FontSize', 18);set(gcf, 'PaperPosition', [0 0 6 4]);set(gcf, 'PaperSize', [ 6 4]);
imagesc(datax*1e-3,datay*1e-3,double(NDWI),'alphadata',wm1==1)
colormap jet;colorbar
%title(['NDWI ',filename(6:13)]);
title(['NDWI water area']);
axis equal
xlabel('x (km)');ylabel('y (km)')
ofile=[filename,'water'];
print('-dpng','-r400',ofile)
saveas(gcf,[ofile],'fig')
end

widthriv=600; %to test
Md1=Md1|edge; %non land area(water or edge);first interpolate; then dilate to save computation time.
wm1 = Md1; %out of region: 1,  water.
%Only use a 300m width land on both sides of the river
Med=imdilate(isnan(NDWI),ones(round(widthriv*3/resrc))); %image edge
Md1 = imdilate(wm1~=0, ones(round(600/resrc))); %make sure land is absolutely land. 300 m width expansion
Mdd = imdilate(Md1, ones(round(widthriv/resrc)));
wm1=~((Mdd-Md1)& ~Med);

ndwil=NDWI(wm1==0); %land
a=ndwil(~isnan(ndwil));
mean2=nanmean(ndwil);std2=nanstd(ndwil);cnt2=sum(~isnan(ndwil));

if flagplot==1
figure;
hold all
set(gcf,'Color','white');set(gca,'FontSize', 18);set(gcf, 'PaperPosition', [0 0 6 4]);set(gcf, 'PaperSize', [ 6 4]);
imagesc(datax*1e-3,datay*1e-3,double(NDWI),'alphadata',wm1==0)
colormap jet;colorbar
%title(['NDWI ',filename(6:13)]);
title(['NDWI land area']);
axis equal
xlabel('x (km)');ylabel('y (km)')
ofile=[filename,'land'];
print('-dpng','-r400',ofile)
saveas(gcf,[ofile],'fig')
end

filenamei=[filename,'i',num2str(i)];
[ithreshold,badflag2,sigma0hat]=adaptivethresNDWI(filename,a,b);
Topt(i)=ithreshold;

close all
end%for i

waterareart=zeros(nc-1,1);
for i=1:nc-1
waterareart(i)=(waterarea(i+1)-waterarea(i))/waterarea(i)*100; % relative change percentage
end

if flagplot==1
figure;plot(mtl(1:nc-1),waterareart,'o-')
hold all;plot(mtl,Topt,'>-')
legend('water area change','Optimal Threshold')
saveas(gcf,'waterareart','fig')
end

[~,idi]=min(abs(waterareart));
mtli=mtl(idi);
Mor=(or>mtli); %20121011
%if strcmp(filename(1:13),'WV01_20110530')
%save test1.mat meanor stdor or
%end
save test1.mat -v7.3

return
end
