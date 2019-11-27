function [mtli]=autothres(or,meanor,stdor)
%automatically searching threshold to find the one that gives the minimum change of water area.
%Wang et al. (2014, https://doi.org/10.1016/j.rse.2014.06.004).
constant

%automated searching of threshold.
mtl=0.1:0.01:100; nc=length(mtl);
waterarea=zeros(size(mtl));
%figure(1)
for i=1:nc
mtli=mtl(i);
Mor=(or>meanor-mtli*stdor)&(or<meanor+mtli*stdor); %20121011 ; slow
waterarea(i)=sum(Mor(:));
%imagesc(Mor)
%title([num2str(i)])
end
waterareart=zeros(nc-1,1);
for i=1:nc-1
waterareart(i)=(waterarea(i+1)-waterarea(i))/waterarea(i)*100; % relative change percentage
end

if flagplot==1
figure;plot(mtl(1:nc-1),waterareart)
end
[~,idi]=min(abs(waterareart));
mtli=mtl(idi);
Mor=(or>meanor-mtli*stdor)&(or<meanor+mtli*stdor); %20121011
if strcmp(name(1:13),'WV01_20110530')
save test1.mat meanor stdor or
end

return
end
