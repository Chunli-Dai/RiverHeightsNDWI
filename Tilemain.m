% main program for getting coastline for each ArcticDEM tile
% Requirements: gdal software
% %%%% inputs needed
macdir='/Users/chunlidai/surge/';
macdir=[];

currentdir=pwd;
addpath(genpath(currentdir));
addpath(genpath([currentdir,'/../rivergithub2/']));
% addpath(genpath([macdir,'/data/chunli/scripts/']));
% addpath(genpath([macdir,'/data/chunli/coastline/codec2/mapformats/']));

constant 
if 0
multidir=[macdir,'/data1/pgc_projects/dai_aleutians_multi_mono/imagery/WV*/']; % directory of mono multispectral images
stripdir='/*/ArcticDEM/region*/strips/2m/';
end

% %%%% control parameters
width=2e3; %buffer width of the a priori coastline, e.g., 2km.

% Get filelist and boundaries
if 1||flagmono==1
if 0
load('alaska2.mat')
else
filename='monolist'; %'boundaries_reg31.dat';
if ~exist(filename,'file')
   str=sprintf('find  %s -name ''*[0-9].xml'' > %s',deblank(multidir),filename);
  [status, cmdout]=system(str);
end
fid = fopen(filename);
n = linecount(fid);
fid = fopen(filename);
% range=fscanf(fid, '%f', [4, n))';
range=zeros(n,4);XYbg=cell(n,1);
idd=[];
for i=1:n
   ifile=[macdir,fgetl(fid)];
   [demdir,name,ext] =fileparts([strtrim(ifile)]);
   f{i}=[name,ext];
   fdir{i}=[demdir,'/'];%working on two regions 
   satname=f{i}(1:4);

   % get the boundary from xml file
   [XYbi,rangei]=imagebd(ifile);
   range(i,1:4)=rangei;XYbg{i}=XYbi;

   r1=strfind(name,'P1BS');
   if ~strcmp(satname,'WV01') && ~isempty(r1) %File is not WV01 but is panchromatic band.
        idd=[idd;i];
   end

end
display(['demdir=',demdir])
range(idd,:)=[];f(idd)=[];fdir(idd)=[];XYbg(idd)=[];
end %if 0

else
n=1;
range=zeros(n,4);XYbg=cell(n,1);
f=[];fdir=[];
end

% %%% Preparation: get the list of strip files and boundries
filename='boundaries_regall_strip.dat'; %'boundaries_reg31.dat';
filename='striplist.dat'; %'boundaries_reg31.dat';
if ~exist(filename,'file')
   str=sprintf('find  %s -name ''*meta.txt'' > %s',deblank(stripdir),filename);
  [status, cmdout]=system(str);
end
fprintf ('\n Step 0: geting the boundary for all files in the region.\n')
%READ INPUT PARAMETERS; getting the boundaries for all files
% filename='boundaries_reg31.dat';
fid = fopen(filename);
n = linecount(fid);
fid = fopen(filename);
range2=zeros(n,4);XYbg2=cell(n,1);
for i=1:n
   ifile=[macdir,fgetl(fid)];
   [demdir,name,ext] =fileparts([strtrim(ifile)]);
   f2{i}=[name,ext];
   fdir2{i}=[demdir,'/'];%working on two regions 
   satname=f2{i}(1:4);

   % get the boundary from xml file
   [XYbi,rangei]=imagebd(ifile);
   range2(i,1:4)=rangei;XYbg2{i}=XYbi;
end

%load volcano list
file=[macdir,'/data/chunli/scripts/gagelistsv1.txt']; %
file=[macdir,'~/data/chunli/scripts/stations2.gmt']; %
file=['stations2.gmt']; %
%vol=load(file);
fid = fopen(file);
n = linecount(fid)-4;
fid = fopen(file);
for i=1:4
str=fgetl(fid);
end
% range=fscanf(fid, '%f', [4, n))';
vol=zeros(n,2);voltext=cell(n,1);
for i=1:n
   vol(i,1:2)=fscanf(fid, '%f', [2, 1])';ifile=fgetl(fid);
   voltext{i}=deblank(ifile);
end

[nv,~]=size(vol);

fid2 = fopen('stations2cnt.dat','w');
exb=20e3;
% for i=1
for i=32 %1:nv
%   loneq=vol(i,2);lateq=vol(i,1);
    loneq=vol(i,1);lateq=vol(i,2);
%   [xeq,yeq]=polarstereo_fwd(lateq,loneq,[], [],70,-45);
%   rang0=[xeq-exb xeq+exb yeq-exb yeq+exb ];
%   ts=datenum(vol(i,4),vol(i,5),vol(i,6));
%   te=datenum(vol(i,7),vol(i,8),vol(i,9));
    source(1:2)=[loneq,lateq];
    texteq=voltext{i};

    % get the data boundary, rang0, of this DEM tile 
%     rang0=[x x+dx/2 y y+dx/2];
    
%             tic
    if 0
    loneq=-(148+49/60+04/3600);lateq=69+0+57/3600; %
    [xeq,yeq]=polarstereo_fwd(lateq,loneq,[],[],70,-45);
    source(1:2)=[loneq,lateq];
    texteq='Sag';
    end
    
	if i==32 
	fprintf('sag station')
	end
    [Co]=ChangeRiver(i,source,texteq,range,XYbg,f,fdir,range2,XYbg2,f2,fdir2);
	%write count count_non_frozen
fprintf(fid2,'%12.6f %12.6f   %23.15e %23.15e\n',lateq,loneq,Co(2),Co(1));
%             fprintf(['\n Tile ',tilefile])
%             toc
            
end %xid

