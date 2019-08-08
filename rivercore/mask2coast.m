function [M,Modfil]=mask2coast(data)
%input:  data, a water mask, -1 void, 0 land, 1 water, int8
%output: M (logical), the shorelines, a matrix the same size of data.z.
%	Modfil (int8), the filtered mask (removing lakes, small clusters, filled small areas) that generated the shoreline.
constant
Mstrip=data;
Medgs=(Mstrip.z==-1);%isnan(Mstrip.z(:,:));
Med=imdilate(Medgs,ones(4));
Medgs=Med;
Modj=Mstrip.z;Modj(Medgs)=0;

resx=mean(Mstrip.x(2:end)-Mstrip.x(1:end-1));resy=mean(Mstrip.y(2:end)-Mstrip.y(1:end-1));resr=mean([abs(resx),abs(resy)]);
Modj= bwareaopen(Modj, round(lakearea/resr/resr)); %remove small clusters
Modfil = bwareaopen(~Modj, round(cloudarea/resr/resr)); %fill small areas: 1e4*4m^2
Modfil=~Modfil;

%remove lakes
tic
data.z=int8(Modfil);
data=mask2river(data); %remove lakes, no fill
Modfil=int8(data.z);Modfil(Modfil==-1)=0; %int8
Medgs2=(data.z==-1); %keep the edge for lowest.m
fprintf('\nGetting river mask, excluding lakes.')
toc

Md1 = imdilate(Modfil, ones(3));
M=logical(Md1-Modfil);
%M=M&~Medgs; 
M=M&~(Medgs2|Medgs); 

Modfil(Medgs2|Medgs)=-1;  %keep the edge for lowest.m
end
