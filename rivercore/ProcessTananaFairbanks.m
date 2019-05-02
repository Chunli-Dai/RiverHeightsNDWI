function Co=ProcessTananaFairbanks(odir,c,lateq,loneq)
% Filtering and fitting of river height profiles.
% Plot Fig.2b in Dai et al., GRL, 2018.
% Modified based on code by Mike Durand, https://github.com/mikedurand/SmoothRiverElevations
% Modified by Chunli Dai, 2017

% clear all
%  macdir=[''];
% addpath([macdir,'//home/chunli/scripts'])
% addpath([macdir,'//home/chunli/scripts/riverfilter/'])
Co=[];

%
if 0
c=shaperead('tan_cl_Close7.shp');
c=c(1); %One line; .X (longitude), .Y (latitude).
end

for k=1:2 % 1 direct method; 2 imagery-altimetry method.
%get data list
if k==1
    filename='proflist'; %'boundaries_reg31.dat';
    str=sprintf('find  %s -name rivprof*sj*[0-9].dat > %s',deblank(odir),filename);
    [status, cmdout]=system(str);
    fid2 = fopen([deblank(odir),'/gageft.txt'], 'w');
    fid3 = fopen([deblank(odir),'/gageslope.txt'], 'w');

elseif k==2
    filename='proflistb'; %'boundaries_reg31.dat';
    str=sprintf('find  %s -name rivprof*bj*[0-9].dat > %s',deblank(odir),filename);
    [status, cmdout]=system(str);
    fid2 = fopen([deblank(odir),'/gageftb.txt'], 'w');
    fid3 = fopen([deblank(odir),'/gageslopeb.txt'], 'w');
end

fid = fopen(filename);
n = linecount(fid);
fid = fopen(filename);
ymdg=zeros(n,1);
for i=1:n
   ifile=[fgetl(fid)];
   [demdir,name,ext] =fileparts([strtrim(ifile)]);
   f{i}=[name,ext];
   fdir{i}=[demdir,'/'];%working on two regions 
   if length(name)>=19 %rivprofWV02_20150714sj1.dat
       satname{i}=name(8:11);
       ymdg(i)=str2num(name(13:20));
   else    %rivprof20150714sj1.dat
       satname{i}='';
       ymdg(i)=str2num(name(8:15));
   end
end

slope=[];

for i=1:n
    ymd=ymdg(i);
    
    ymds=num2str(ymd);
    mon=str2num(ymds(5:6));
    if 0 & ( mon <=4 || mon>=11)  %winter
    continue
    else
    display(['Run i=',num2str(i),';ymd=',num2str(ymd)])
    end
    
    display(['i=',num2str(i),';ymd ',num2str(ymd)])
    rivprof1=[fdir{i},f{i}];
    
    %assign the loaded profiles to the cell s.
    s{i}=load(rivprof1);
    % s{2}=eval(rivprof2);
end %for i
 
%parameters
reAttach=true; %best to keep set to true
p.distmax=500; %max distance from centerline to include pixels
p.ncut=20;     %min # of pixels required for each node
p.StdMax=2.5;  %max pixel standard deviation to still include
p.pct=25;      %percentile of pixels chosen to represent that node
p.dx=.1;       %node spacing in km
%p.dx=.5;       %500m not working good.
p.x=1:p.dx:120; %120km is the total length of the centerline.
p.N=length(p.x);
p.HrefCut=1;   %cutoff for maximum elevation difference allowed from "reference" (i.e. first pass) 

if ~exist('Elevations','dir')
  mkdir('Elevations')
end

%% first pass
for i=1:n
% Href{1}=nan(1,p.N); Href{2}=nan(1,p.N);
Href{i}=nan(1,p.N);
end
[Est]=ProcessData(c,s,reAttach,p,Href,'Fairbanks','LP',lateq,loneq);

%% second pass
reAttach=true;
% Href{1}=Est{1}.Hc; Href{2}=Est{2}.Hc;
for i=1:n
    Href{i}=Est{i}.Hc;
end
[Est,Data,xobs]=ProcessData(c,s,reAttach,p,Href,'Fairbanks','SLM',lateq,loneq); % 5 km might be too short for this
fprintf(['\n Gage centerline coordinate xobs:',num2str(xobs),'\n']);
%%
% xobs=93; % from CL #7 from Rui ;  USGS gage loneq=-147.8389; lateq=64.7928
% using Vdatum, 404.93' above NAVD88 = 400.67' above EGM08
% Hobs=([21.9 18.15]+400.67).*.3048; %USGS gage data, above EGM08
dHxform=0;%-11.23; % using Vdatum, 0 m WGS84 = -11.23 m EGM08

for i=1:n
        ymd=ymdg(i);
        if strcmp(satname{i},'WV01')
            satflag=1;
        else
            satflag=0;
        end
            
if 0
%% Plotting 
xmin=min([Data{1}.lon; Data{2}.lon;]);
xmax=max([Data{1}.lon; Data{2}.lon;]);
ymin=min([Data{1}.lat; Data{2}.lat;]);
ymax=max([Data{1}.lat; Data{2}.lat;]);

%% Plot the shoreline location.
figure(1)
subplot(311)
mapshow(c); hold on;
plot(Data{1}.lon,Data{1}.lat,'.',Data{1}.lon(Data{1}.iuse),Data{1}.lat(Data{1}.iuse),'.'); hold off;
axis([xmin xmax ymin ymax]); hold off;
legend('All','Use')
title('rivp20110804')
subplot(312)
mapshow(c); hold on;
plot(Data{1}.lon,Data{1}.lat,'.',Data{2}.lon(Data{2}.iuse),Data{2}.lat(Data{2}.iuse),'.'); hold off;
axis([xmin xmax ymin ymax]); hold off;
title('rivp20121011')
subplot(313)
mapshow(c); hold on;
plot(Data{1}.lon(Data{1}.iuse),Data{1}.lat(Data{1}.iuse),'.',Data{2}.lon(Data{2}.iuse),Data{2}.lat(Data{2}.iuse),'.'); hold off;
hold on;plot( - (147+50/60+20/3600), 64+47/60+34/3600, 'bs','Markersize',12) %gauge
axis([xmin xmax ymin ymax]); hold off

%% Plot the river elevation profiles
figure(2)
subplot(411)
h311=plot(Data{1}.FDh(Data{1}.iuse),Data{1}.h(Data{1}.iuse),'.',p.x,Est{1}.Hhat,'o');
set(h311(1),'Color',[.7 .7 .7])
title('August 2011')
subplot(412)
h312=plot(Data{2}.FDh(Data{2}.iuse),Data{2}.h(Data{2}.iuse),'.',p.x,Est{2}.Hhat,'o');
set(h312(1),'Color',[.7 .7 .7])
title('October 2012')
subplot(413)
plot(p.x,Est{1}.Hhat,'o',p.x,Est{2}.Hhat,'o','LineWidth',2)
set(gca,'FontSize',14)
grid on; legend('August 2011','October 2012','Location','Best')
subplot(414)
plot(p.x,Est{1}.HhatStd,'o',p.x,Est{2}.HhatStd,'o','LineWidth',2)
grid on
set(gca,'FontSize',14)
end
%% Plot the fitting curve. Fig. 2b
figure %(3);
set(gcf,'Color','white')
hold all
han1=plot(p.x(Est{i}.Use),Est{i}.Hhat(Est{i}.Use)+dHxform,'o');
% plot(p.x(Est{2}.Use),Est{2}.Hhat(Est{2}.Use)+dHxform,'o'); hold on;
% han2=plot(p.x(Est{2}.iSolUse),Est{2}.Hc(Est{2}.iSolUse)+dHxform,'g');
plot(p.x(Est{i}.iSolUse),Est{i}.Hc(Est{i}.iSolUse)+dHxform,'r','LineWidth',2); 
% han3=plot(xobs,flip(Hobs)+1.23,'rs','LineWidth',2,'MarkerSize',8);% hold off;
y1=Est{i}.Hc(Est{i}.iSolUse)+dHxform;
if ~isempty(y1)
plot([xobs,xobs],[min(y1(:)) max(y1(:))],'g-')
end
set(gca,'FontSize',16)
% set(han1(1:2),'Color',[.7 .7 .7])
% set(han2(1),'Color',[0.93 0.69 0.13])
% set(han2(2),'Color',[0.49 0.18 0.56])
% set(han3(1),'Color',get(han2(1),'Color'))
% set(han3(2),'Color',get(han2(2),'Color'))
% set(han1(2),'Color',get(han2(1),'Color'))
% set(han1(1),'Color',get(han2(2),'Color'))
% legend([han2; han3;],'ArcticDEM 11 Oct 2012','ArcticDEM 4 Aug 2011','USGS 11 Oct 2012','USGS 4 Aug 2011','Location','Best')
xlabel('Main channel flow distance, km')
ylabel('Water surface elevation (EGM08), m')
% set(gca,'XLim',[85 94])
box on
rivprof1=[fdir{i},f{i}];
ofile=strrep(rivprof1,'.dat','ft');
title(f{i})
saveas(gcf,ofile,'fig')

ofile2=strrep(rivprof1,'.dat','ft.dat'); %rivprofWV02_20180626bj83ft.dat rivprofWV02_20180626bj83ft.fig
profx=p.x(Est{i}.iSolUse);profy=Est{i}.Hc(Est{i}.iSolUse);
output=[profx(:)*1e3,profy(:)]; %m m
save(ofile2,'output','-ascii')

% xming(i)=min(p.x(Est{1}.iSolUse));xmaxg(i)=max(p.x(Est{1}.iSolUse));
% 
% xmin=max([xmin min(p.x(Est{i}.iSolUse))]);
% xmax=min([xmax max(p.x(Est{i}.iSolUse))]);

% xminb=max([xminb min(p.x(Est{2}.iSolUse))]);
% xmaxb=min([xmaxb max(p.x(Est{2}.iSolUse))]);

if isfield(Est{i}, 'xp')
 %use fitting line for average
    gagereach=1;%1km
    if 0 %Do not use this. This keep the near ends of fitting function, which can be curved/ inaccurate.
    id=find(abs(Est{i}.xp-xobs)<=gagereach/2.);nn=length(id);
    % idb=find(abs(Est{2}.xp-xobs)<=500e-3);nnb=length(idb);
    hest=Est{i}.yp(id);heststd=Est{i}.stdres*ones(size(hest));
    % hestb=Est{2}.yp(idb);heststdb=Est{2}.stdres*ones(size(hestb));
    xsr=Est{i}.xp(id);xsrdx=Est{1}.xp(2)-Est{1}.xp(1);lenr=nn*abs(xsrdx);%length of selected reach
    if (lenr-gagereach)<-gagereach*0.1; 
    warning(['i=',num2str(i),';',num2str(ymd),' This profile does not have a good coverage of the gage within ',num2str(gagereach),' km']);
    continue;
    end
    else %Use this, which trim off the ends.
    x2=p.x(Est{i}.iSolUse);y2=Est{i}.Hc(Est{i}.iSolUse);
    id=find(abs(x2-xobs)<=gagereach/2.);nn=length(id);
    hest=y2(id);heststd=Est{i}.stdres*ones(size(hest));
    xsrdx=p.dx;lenr=nn*abs(xsrdx);
    if (lenr-gagereach)< -gagereach*0.1
    warning(['i=',num2str(i),';',num2str(ymd),' This profile does not have a good coverage of the gage within ',num2str(gagereach),' km']);
    continue;
    end
end
 hgage=nanmean(hest);hstd=Est{i}.stdres; %1./nn*sqrt(sum(heststd.^2));

% hgageb=mean(hestb);hstdb=Est{2}.stdres;  %1./nnb*sqrt(sum(heststdb.^2));

% 92.9000  100.1000; 7.2km
% xl=92.9;xr=100.1;
% id=find(Est{1}.xp>=xl&Est{1}.xp<=xr);
% idb=find(Est{2}.xp>=xl&Est{2}.xp<=xr);
slope(1,i)=[nanmean(Est{i}.Slope)];% mean(Est{2}.Slope(idb))];

fprintf(fid2,'%d %12.6f %12.6f %d \n',ymd,hgage,hstd,satflag); %1 WV01; 0 otherwise
fprintf(fid3,'%d %12.6f \n',ymd,slope(1,i));
end

close all

end %for i=1:n
% close(fid2)
% close(fid3)

end %for k

return
end

