function [datestr]=strip2date(stripmetafile)
%given strip file name, to get date YYYYMMDDHHMMSS

c=textread(stripmetafile,'%s','delimiter','\n');
rs=find(~cellfun(@isempty,strfind(c,'scene')));
nsce=length(rs)-1; %number of scenes
rmsmeta=zeros(nsce,1);idd=[];dzxyd=zeros(nsce,3);mfile=cell(nsce,1);
for i=1 %nsce
c1=c{rs(1)+i};r1=strfind(c1,'dem.tif');c1(1:(r1+6))='';
tmp=sscanf(c1, '%g', 4);
rmsmeta(i)=tmp(1);dzxyd(i,1:3)=tmp(2:4);

c1=c{rs(1+i)+3};r1=strfind(c1,'Image 1=');if(isempty(r1)) Warning('Image 1 not found') ; end
r1=strfind(c1,'/');c1(1:(r1(end)))='';

c11=deblank(c1);
satname=c11(1:4);%use this, since strip could be W1W2_20110825_1020010014850E00_103001000C5BD600_seg1_2m_meta.txt

%WV02_20160304214247_1030010052B75A00_16MAR04214247
datestr=c1(6:19);

if strcmp(satname,'WV01')
    mfile{i}=deblank(c1);
else
    mfile{i}=deblank(strrep(c1,'P1BS','M1BS'));%e.g.WV02_20160304214247_1030010052B75A00_16MAR04214247-M1BS-500641617080_01_P009.tif
end

end


return
end
