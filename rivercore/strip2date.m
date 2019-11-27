function [datestr]=strip2date(stripmetafile)
%given strip file name, to get date YYYYMMDDHHMMSS

c=textread(stripmetafile,'%s','delimiter','\n');
rs=find(~cellfun(@isempty,strfind(c,'scene')));

if ~isempty(rs) %strip DEM
nsce=length(rs)-1; %number of scenes
for i=1 %nsce
% c1=c{rs(1)+i};r1=strfind(c1,'dem.tif');c1(1:(r1+6))='';
%new release scene name is '_dem_smooth.tif'
c1=c{rs(1)+i};
matchStr = regexp(c1,'.tif','split');c1=matchStr{end};
tmp=sscanf(c1, '%g', 4);

c1=c{rs(1+i)+3};r1=strfind(c1,'Image 1=');if(isempty(r1)) Warning('Image 1 not found') ; end
r1=strfind(c1,'/');c1(1:(r1(end)))='';

c11=deblank(c1);
satname=c11(1:4);%use this, since strip could be W1W2_20110825_1020010014850E00_103001000C5BD600_seg1_2m_meta.txt

%WV02_20160304214247_1030010052B75A00_16MAR04214247
datestr=c11(6:19);
end

else %rs is empty; scene DEM meta file
% compatible for scene DEM meta file.
check={'Image 1=','/','.tif'};
Mcheck=contains(c,check{1})&contains(c,check{2})&contains(c,check{3});
c1all=c(Mcheck);
nsce=sum(Mcheck(:));

for i=1  %:nsce
c1=c1all{i};
r1=strfind(c1,'/');c1(1:(r1(end)))='';

c11=deblank(c1);
satname=c11(1:4);%use this, since strip could be W1W2_20110825_1020010014850E00_103001000C5BD600_seg1_2m_meta.txt
datestr=c11(6:19);

end
end %if ~isempty(rs)



return
end
