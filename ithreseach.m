   
% Pan chromatic images.
% Get the brightness threshold for each a priori water clusters;time consuming.


constant

BW=wm1;
data.x=mt.x;data.y=mt.y;

if 1
resr=abs(data.x(2)-data.x(1));
width2=100; %river width m
narea=round(width2*width2/resr/resr);
%remove small clusters
Modj= bwareaopen(BW, narea);BW=Modj;
end
        
CC = bwconncomp(BW);

[X,Y]=meshgrid(data.x,data.y);

figure;imagesc(data.x*1e-3,data.y*1e-3,or);colorbar

BW2=false(size(BW)); 
% get the threshold for each a priori water clusters;time consuming.
for k=1:CC.NumObjects

        BW3=BW;BW3(:)=0;
        BW3(CC.PixelIdxList{k})=BW(CC.PixelIdxList{k});

        ndwil=double(or(BW3)); %water
        mean1=nanmean(ndwil);std1=nanstd(ndwil);cnt1=sum(~isnan(ndwil));
        
        badflag=0;
        if (isnan(mean1)||cnt1<cntmin)
            badflag=1;continue
        else
            meanor=mean1;stdor=std1;
            Mor=(or>meanor-stdor)&(or<meanor+stdor); 
        end
        
        if 1
            hold on;plot(X(BW3)*1e-3,Y(BW3)*1e-3,'r.')            
            figure;imagesc(data.x*1e-3,data.y*1e-3,Mor);colorbar;title(['k=',num2str(k)])
        end
        
        BW2=BW2|Mor;

end
