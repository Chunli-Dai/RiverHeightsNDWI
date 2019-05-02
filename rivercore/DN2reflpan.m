function rho=DN2reflpan(data,metafile)
% transform DN to reflectance for panchromatic band only.
% Refers to landslide/workYukon2/Landslide.m

[~,filename,ext]=fileparts(metafile);

%Check xml file or meta file
name=filename;
name1=[name(end-3:end),ext];
satname=filename(1:4);

% To check whether the data is mono or scene files
flagfmt=0;
if strcmp(name1,'meta.txt') %Pan band only
	flagfmt=3; %stereo files
elseif strcmp(ext,'.xml') % Could be pan band or multispectral band
	flagfmt=1; %xml files mono
end

%check whether the data is multispectral image
flagfmt2=0; %panchromatic band
r1=strfind(name,'M1BS');
if ~isempty(r1); flagfmt2=1;end

if flagfmt2==0 %pan band
%assume in each strip, all scenes and each image in pairs have the same
%effectivebandwith and abscalfactor.
if flagfmt==3
c=textread(metafile,'%s','delimiter','\n');
strg={['_abscalfact='],['_effbw='],['_Mean_sun_elevation=']};
for j=1:3
str=strg(j);
r=find(~cellfun(@isempty,strfind(c,str)));
str=c{r(1)};r1=strfind(str,'=');str(1:r1)='';
% Xbs=deblank(strrep(c1{r(1)},str,''));
t1(j)= sscanf(str, '%g', 1);
end
abscalfactor=t1(1);effectivebandwith=t1(2);meanSunEl=t1(3);

elseif flagfmt==1 %mono image
  iPan=1;
  ig=[iPan];
  bandstr={'BAND_P'};
  abscalfactor=zeros(1,1);effectivebandwith=zeros(1,1);

%Refers to rivergithub2/multispecmono.m
%get ABSCALFACTOR EFFECTIVEBANDWIDTH for each band, and MEANSUNEL
c=textread(metafile,'%s','delimiter','\n');
% fid=fopen(metafile);
% n = linecount(fid);
% c=cell(n,1);
% for i=1:n
%   c{i} = fgetl(fid);
% end
r=find(~cellfun(@isempty,strfind(c,'MEANSUNEL')));
c2=c{r};
r1=strfind(c2,'>');r2=strfind(c2,'</');c2([1:r1(1),r2(1):end])='';
meanSunEl=sscanf(c2, '%g', 1);
for i=ig % bands 1 to 8
r=find(~cellfun(@isempty,strfind(c,bandstr(i))));
%/data1/pgc_projects/coastline/imagery/WV02_20120930_103001001A971B00_103001001C13FE00/WV02_20120930220825_103001001A971B00_12SEP30220825-M1BS-052903623070_01_P001.xml
if isempty(r)
warning(['Bands are different as anticipated.'])
return;
end %
c1=c(r(1):(r(2)));
str='ABSCALFACTOR';
r=find(~cellfun(@isempty,strfind(c1,str)));
c2=c1{r};r1=strfind(c2,'>');r2=strfind(c2,'</');c2([1:r1(1),r2(1):end])='';
% Xbs=deblank(strrep(c1{r(1)},str,''));
abscalfactor(i)= sscanf(c2, '%g', 1);
str='EFFECTIVEBANDWIDTH';
r=find(~cellfun(@isempty,strfind(c1,str)));
c2=c1{r};r1=strfind(c2,'>');r2=strfind(c2,'</');c2([1:r1(1),r2(1):end])='';
effectivebandwith(i)=sscanf(c2, '%g', 1);
end %for i
end %if 

Theta=(90.-meanSunEl)*pi/180.;
dES=1.; %AU
% sun-earth distance polynomial function coefficients
% doy = day of the year: 
year=sscanf(filename(6:9), '%g', 1); month=sscanf(filename(10:11), '%g', 1); day=sscanf(filename(12:13), '%g', 1);
doy=juliandate(year,month,day)-juliandate(year,1,1)+1;
C = [1.8739e-26,-3.4455e-23,2.7359e-20,-1.2296e-17,3.0855e-15,-2.2412e-13,-5.8744e-11,6.9972e-10,2.5475e-06,-1.6415e-05,0.9833];
dES = polyval(C,doy); % earth-sun distance

satname=filename(1:4);
[GAIN,OFFSET,Esun]=readgainoffset(satname); %read Gain OFFset data, GainOffset.txt.

DN=double(data.z(:,:));
L=GAIN(end)*DN*(abscalfactor/effectivebandwith)+OFFSET(end);
rho=L*dES^2*pi/(Esun(end)*cos(Theta));

if isfield(data, 'x')
tt=rho;clear rho
rho.x=data.x;rho.y=data.y;rho.z=tt;
end

elseif  flagfmt2==1 %multi band; multispecmono.m
        [~,~,nb]=size(data.z);
    if nb == 4 %the classic 4 bands        
        iBlue=1;iGreen=2;iRed=3;iNIR1=4;
        ig=[iBlue,iGreen,iRed,iNIR1];
        bandstr={'BAND_B','BAND_G','BAND_R','BAND_N'};
        abscalfactor=zeros(4,1);effectivebandwith=zeros(4,1);
        flagb=1;

    elseif nb ==8 
	% the 8 bands for WV02 WV03
        iNIR2=8;iCoastal=1;iRed=5;iGreen=3; iNIR1=7;
        iBlue=2;iYellow=4;iRE=6;
        ig=[iCoastal,iGreen,iRed,iNIR2,iNIR1,iBlue,iYellow,iRE];
        bandstr={'BAND_C','BAND_B','BAND_G','BAND_Y','BAND_R','BAND_RE','BAND_N','BAND_N2'};
        abscalfactor=zeros(8,1);effectivebandwith=zeros(8,1);
        flagb=2;
    else
       warning(['Bands are different as anticipated for ',metafile])
    end
    [GAIN,OFFSET,Esun]=readgainoffset(satname); %read Gain OFFset data, GainOffset.txt.
    % read .xml meta file
    % metafile=[macdir,'/data3/ArcticDEM/region_31_alaska_south/strips/2m/WV03_20170414_104001002CD4E700_104001002BA9FA00_seg2_2m_meta.txt'];
       %get ABSCALFACTOR EFFECTIVEBANDWIDTH for each band, and MEANSUNEL
    c=textread(metafile,'%s','delimiter','\n');

    r=find(~cellfun(@isempty,strfind(c,'MEANSUNEL')));
    c2=c{r};
    r1=strfind(c2,'>');r2=strfind(c2,'</');c2([1:r1(1),r2(1):end])='';
    meanSunEl=sscanf(c2, '%g', 1);
    for i=ig % bands 1 to 8
    r=find(~cellfun(@isempty,strfind(c,bandstr(i))));
    %/data1/pgc_projects/coastline/imagery/WV02_20120930_103001001A971B00_103001001C13FE00/WV02_20120930220825_103001001A971B00_12SEP30220825-M1BS-052903623070_01_P001.xml
    if isempty(r)
    warning(['Bands are different as anticipated.'])
    return;
    end %
    c1=c(r(1):(r(2)));
    str='ABSCALFACTOR';
    r=find(~cellfun(@isempty,strfind(c1,str)));
    c2=c1{r};r1=strfind(c2,'>');r2=strfind(c2,'</');c2([1:r1(1),r2(1):end])='';
    % Xbs=deblank(strrep(c1{r(1)},str,''));
    abscalfactor(i)= sscanf(c2, '%g', 1);
    str='EFFECTIVEBANDWIDTH';
    r=find(~cellfun(@isempty,strfind(c1,str)));
    c2=c1{r};r1=strfind(c2,'>');r2=strfind(c2,'</');c2([1:r1(1),r2(1):end])='';
    effectivebandwith(i)=sscanf(c2, '%g', 1);
    end

    Theta=(90.-meanSunEl)*pi/180.;
    dES=1.; %AU
    % sun-earth distance polynomial function coefficients
    % doy = day of the year: 
    year=sscanf(filename(6:9), '%g', 1); month=sscanf(filename(10:11), '%g', 1); day=sscanf(filename(12:13), '%g', 1);
    doy=juliandate(year,month,day)-juliandate(year,1,1)+1;
    C = [1.8739e-26,-3.4455e-23,2.7359e-20,-1.2296e-17,3.0855e-15,-2.2412e-13,-5.8744e-11,6.9972e-10,2.5475e-06,-1.6415e-05,0.9833];
    dES = polyval(C,doy); % earth-sun distance

    rhog=zeros(size(data.z));
    for i=ig %band id 
        DN=double(data.z(:,:,i));
        L=GAIN(i)*DN*(abscalfactor(i)/effectivebandwith(i))+OFFSET(i);
        rho=L*dES^2*pi/(Esun(i)*cos(Theta));
        rhog(:,:,i)=rho;
    end
    clear rho
    if isfield(data, 'x')
       rho.x=data.x;rho.y=data.y;rho.z=rhog;
    else 
       rho=rhog;
    end
    
end % flagfmt2

return
end
