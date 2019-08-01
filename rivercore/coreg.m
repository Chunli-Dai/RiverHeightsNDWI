function [dX4Sg,idregion,data0r]=coreg(rang0,idregion,XYbg,f,fdir)
close all
constant

% find the DEM mosaics within rang0.
xr0=[rang0(1) rang0(2) rang0(2) rang0(1) rang0(1) ];yr0=[rang0(4) rang0(4) rang0(3) rang0(3) rang0(4) ];

resrc=40;
ranget=round(rang0/resrc)*resrc;rang0=ranget;
tx=ranget(1):resrc:ranget(2);ty=ranget(4):-resrc:ranget(3);
xout=tx;yout=ty;
data0r.x=xout;data0r.y=yout;data0r.z=nan(length(yout),length(xout));

fprintf ('\n Step 1.1a: getting the tile DEM as a reference coordinate system.')
% get the data boundary, rang0, of this DEM tile 
dx=100e3;x0=-4000e3;y0=-4000e3;%xe=3400e3;ye=4000e3; %ArcticDEM Mosaic tiles coordinate reference;
% x=x0+(xid-1)*dx+(xids-1)*dx/2;y=y0+(yid-1)*dx+(yids-1)*dx/2; %bottom left
% rang0b=[x x+dx/2 y y+dx/2]; %exact tile boundary;
icount=0;clear xc yc
for i=[1,3]
    icount=icount+1;
    x=xr0(i);y=yr0(i);
    xid=floor((x-x0)/dx)+1;
    yid=floor((y-y0)/dx)+1;
    xids=floor((x-x0-(xid-1)*dx)/(dx/2))+1;
    yids=floor((y-y0-(yid-1)*dx)/(dx/2))+1;
    xc(icount)=x0+(xid-1)*dx+(xids-1)*dx/2; yc(icount)=y0+(yid-1)*dx+(yids-1)*dx/2; %bottom left
%     tilefile=sprintf('%02d_%02d_%01d_%01d_5m_v2.0_reg_dem.tif',yid,xid,xids,yids)
%     x=xc(icount);y=yc(icount);
%     rang0b=[x x+dx/2 y y+dx/2];
%     xb0=[rang0b(1) rang0b(2) rang0b(2) rang0b(1) rang0b(1) ];
%     yb0=[rang0b(4) rang0b(4) rang0b(3) rang0b(3) rang0b(4) ];
%     hold on;plot(xb0*1e-3,yb0*1e-3,'.-')
end
icount=0;
for x=min(xc):dx/2:max(xc)
    for y=min(yc):dx/2:max(yc)
        icount=icount+1;
    xid=floor((x-x0)/dx)+1;
    yid=floor((y-y0)/dx)+1;
    xids=floor((x-x0-(xid-1)*dx)/(dx/2))+1;
    yids=floor((y-y0-(yid-1)*dx)/(dx/2))+1;
    tilefile=sprintf('%02d_%02d_%01d_%01d_2m_v3.0_reg_dem.tif',yid,xid,xids,yids);
    %check
%     xp=x0+(xid-1)*dx+(xids-1)*dx/2;yp=y0+(yid-1)*dx+(yids-1)*dx/2; %bottom left
%     rang0b=[xp xp+dx/2 yp yp+dx/2];
%     xb0=[rang0b(1) rang0b(2) rang0b(2) rang0b(1) rang0b(1) ];
%     yb0=[rang0b(4) rang0b(4) rang0b(3) rang0b(3) rang0b(4) ];
%     hold on;plot(xb0*1e-3,yb0*1e-3,'o-')
%     pause
    
        % find the dem tile file or download the data
    [status , cmdout ]=system(['find ',tiledir,' -name ',tilefile]); %status always 0, cmdout could be empty.
    if ~isempty(cmdout) && status ==0 % 
        tilefile=deblank(cmdout);
    else
        warning(['Tile file ',tilefile,' not found! Download it from website'])
        [dir,name,ext] =fileparts(tilefile);
    %     http://data.pgc.umn.edu/elev/dem/setsm/ArcticDEM/mosaic/v2.0/43_59/43_59_2_1_5m_v2.0.tar
        tarfile=[name(1:17),'.tar'];
%       webtile=deblank(['http://data.pgc.umn.edu/elev/dem/setsm/ArcticDEM/mosaic/v2.0/',name(1:5),'/',tarfile,'.gz']);
        webtile=deblank(['http://data.pgc.umn.edu/elev/dem/setsm/ArcticDEM/mosaic/v3.0/2m/',name(1:5),'/',tarfile,'.gz']);
        system(['wget  ',webtile,' -o downloadlog'])
	system(['gunzip -f ',tarfile,'.gz']);
        system(['tar -xvf ',tarfile]);
        %collect all downloaded mosaic dems to tiledirnew/ % to do
        %system(['mv *dem_meta.txt  *reg.txt *.tar *_reg_dem.tif *_reg_matchtag.tif ',tiledirnew]);
        system(['cp ',name(1:5),'* ', tiledirnew]);
        [status , cmdout ]=system(['find ./ -name ',tilefile]);
        %[status , cmdout ]=system(['find ',tiledirnew,' -name ',tilefile]);
	C = strsplit(strtrim(cmdout));
        %tilefile=deblank(cmdout);
        tilefile=deblank(C{1});
        if ~exist(tilefile,'file')
          warning(['Tile file ',tilefile,' not found! '])
          continue
        end
    end

    data=readGeotiff(tilefile,'map_subset',rang0);
    data.z(data.z == -9999) = NaN; % %convert nodata values to nan for interpolation
    tz =interp2(data.x,data.y,double(data.z),xout,yout','*linear',nan);%toc;%4s
    M=~isnan(tz);
    data0r.z(M)=tz(M);
%     dataref(icount)=data;
    end
end
% data0r.z(isnan(data0r.z))=-9999; %return to -9999

if flagplot==1
hills=hillshade(double(data0r.z),data0r.x,data0r.y,'plotit');
set(gcf,'Color','white');set(gca,'FontSize', 12);set(gcf, 'PaperPosition', [0.25 2.5 4 3]);
colorbar;%caxis([0 250])
hold on
title(['Reference DEM tile mosaic']);
axis equal
saveas(gcf,'RefDEM','fig')
end

if 0 % old way
data=readGeotiff(tilefile,'map_subset',rang0);
%   reduce resolution
ranget=[[min(data.x) max(data.x) min(data.y) max(data.y)]/resrc];
ranget=[ceil(ranget(1)) floor(ranget(2)) ceil(ranget(3)) floor(ranget(4))]*resrc;
tx=ranget(1):resrc:ranget(2);ty=ranget(4):-resrc:ranget(3);
res=data.info.map_info.dx;
nsr=resrc/res;
if nsr <1 
 data.z(data.z == -9999) = NaN; % %convert nodata values to nan for interpolation
 tz =interp2(data.x,data.y,double(data.z),tx,ty','*linear');%toc;%4s
 tz(isnan(tz))=-9999; %return to -9999
%the following get the exactly the same as from interp2
else
idrxs=find(abs(data.x-ranget(1))<1e-3);idrxe=find(abs(data.x-ranget(2))<1e-3);
idrys=find(abs(data.y-ranget(4))<1e-3);idrye=find(abs(data.y-ranget(3))<1e-3);
idrx=idrxs:nsr:idrxe;
idry=idrys:nsr:idrye;
dd=[idrx(2:end)-idrx(1:end-1),idry(2:end)-idry(1:end-1)];
if isempty(dd)||any(abs(dd-nsr)>1e-3);warning('Resized grid is not constant spacing.');end
if length(tx)~=length(idrx)||length(ty)~=length(idry);warning('Wrong downsizing of DEM.');end
tz=data.z(idry,idrx); %0.09s
end
datar= struct();
datar.x=tx;datar.y=ty;  datar.z=tz;
data0r=datar;

%from 54_06_2_2_5m_v2.0_reg_dem.tif get 54_06_5m_v2.0_reg.txt
[dir,name,ext] =fileparts(tilefile);
infile=[name(1:5),name(10:21),'.txt'];%tile reg file
infile=[dir,'/',infile];
if exist(infile)
c=textread(infile,'%s','delimiter','\n');
r=find(~cellfun(@isempty,strfind(c,'Translation')));
c1=c{r(1)};r1=strfind(c1,'=');c1(1:r1)='';dzxys=c1;
dzxyt=sscanf(dzxys,'%f, %f, %f');
dx4=dzxyt; %z, x,y%(2:3); %save global variable
dx4=zeros(3,1); % do not use the tile reg.txt data.
%Mean Vertical Residual
r=find(~cellfun(@isempty,strfind(c,'Mean Vertical Residual')));
c1=c{r(1)};r1=strfind(c1,'=');c1(1:r1)='';
rmstile=sscanf(c1, '%g', 1);
else
    Warning(['Tile Reg file do not exist:',infile])
end
end % if 0
dx4=zeros(3,1); % do not use the tile reg.txt data.

fprintf ('\n Step 1.1b: Coregister all strip DEMs to the reference DEM tile.')
tic
dzxyd=zeros(length(idregion),3);
% dxov=resr*10; %grid size for approximating overlapping area of polygons.
idreg=find(dzxyd(:,1)~=0);  %reg file
idregn=find(dzxyd(:,1)==0); %no reg file
pg=zeros(length(idregion),3);dX4Sg=pg;rmsreg2=zeros(length(idregion),1);
idd=[];idd2=[];
for j=1:length(idregion) % %hi
    %idreg' %idreg([2 5 6 7 10 12 14 15 18])'%[idregn(:)]'  %files that have no reg.txt file.
    i=idregion(j);
    Xb2=XYbg{i}(:,1);
    Yb2=XYbg{i}(:,2);
    
    demdir=fdir{i};
    infile= strrep([demdir,'/',f{i}],'meta.txt','dem.tif') 
    texttar=f{i}(1:13);

    data=readGeotiff(infile);
%get the mask
    infile= strrep([demdir,'/',f{i}],'meta.txt','bitmask.tif');
    if exist(infile,'file') 
    mask=readGeotiff(infile); %0 is good; 1 is edge; > 1 bad data.
    data.z(mask.z>0)=-9999;
    end

    % ratio of good DEMs for whole scene
    res=data.info.map_info.dx;
    XYbi=XYbg{i};
    Xb=XYbi(:,1);Yb=XYbi(:,2);
    if 0 %super slow
    [X,Y]=meshgrid(data.x,data.y);
    in = inpolygon(X,Y,Xb,Yb); % 1 inside the polygon
    nptsubrt=sum(sum(data.z~=-9999))/(sum(in(:))) 
	clear X Y in
    else %fast
	idx=round((Xb-data.x(1))/res)+1;
	idy=round(-(Yb-data.y(1))/res)+1;
	in= poly2mask(idx,idy, length(data.y),length(data.x)); % faster than inpolygon
    	nptsubrt=sum(sum(data.z~=-9999))/(sum(in(:)))
    end
	
    if nptsubrt<0.5 
        idd2=[idd2;j];
	continue
    end

%   reduce resolution
    ranget=[[min(data.x) max(data.x) min(data.y) max(data.y)]/resrc];
    ranget=[ceil(ranget(1)) floor(ranget(2)) ceil(ranget(3)) floor(ranget(4))]*resrc;
    tx=ranget(1):resrc:ranget(2);ty=ranget(4):-resrc:ranget(3);
    res=data.info.map_info.dx;
    nsr=resrc/res;
    if nsr <1 
     data.z(data.z == -9999) = NaN; % %convert nodata values to nan for interpolation
     tz =interp2(data.x,data.y,double(data.z),tx,ty','*linear');%toc;%4s
     tz(isnan(tz))=-9999; %return to -9999
    %the following get the exactly the same as from interp2
    else
    idrxs=find(abs(data.x-ranget(1))<1e-3);idrxe=find(abs(data.x-ranget(2))<1e-3);
    idrys=find(abs(data.y-ranget(4))<1e-3);idrye=find(abs(data.y-ranget(3))<1e-3);
    idrx=idrxs:nsr:idrxe;
    idry=idrys:nsr:idrye;
    dd=[idrx(2:end)-idrx(1:end-1),idry(2:end)-idry(1:end-1)];
    if isempty(dd)||any(abs(dd-nsr)>1e-3);warning('Resized grid is not constant spacing.');end
    if length(tx)~=length(idrx)||length(ty)~=length(idry);warning('Wrong downsizing of DEM.');end
    tz=data.z(idry,idrx); %0.09s
    end
    datar= struct();
    datar.x=tx;datar.y=ty;  datar.z=tz;
    
    rangtar=[min(datar.x) max(datar.x) min(datar.y) max(datar.y)];
    rangref=[min(data0r.x) max(data0r.x) min(data0r.y) max(data0r.y)];
    rangeov=[max(rangtar(1),rangref(1)),min(rangtar(2),rangref(2)), max(rangtar(3),rangref(3)),min(rangtar(4),rangref(4))];
    %crop the overlapping data
    refdem=[];tardem=[];
    idx=find(data0r.x>=rangeov(1) & data0r.x<=rangeov(2));
    idy=find(data0r.y>=rangeov(3) & data0r.y<=rangeov(4));
    refdem.x=data0r.x(idx);refdem.y=data0r.y(idy);
    refdem.z=data0r.z(idy,idx);
    idx=find(datar.x>=rangeov(1) & datar.x<=rangeov(2));
    idy=find(datar.y>=rangeov(3) & datar.y<=rangeov(4));
    tardem.x=datar.x(idx);tardem.y=datar.y(idy);
    tardem.z=datar.z(idy,idx);
    [n,m]=size(tardem.z);minm=min(n,m);
    if sum(size(tardem.z)~=size(refdem.z))||minm<3
        warning(['Wrong overlapping crop, check i:',num2str(idregion(j))]);
        idd=[idd;j];
        continue
    end

      datatarz=tardem.z;datatarz(datatarz== -9999) = NaN; 
      datarefz=refdem.z;datarefz(datarefz== -9999) = NaN;
      iter= 1;
    %  	[dx,dy,sigma]=vmimc(infile1,infile2) 
      [z2out,p,sigma] = coregisterdems(refdem.x,refdem.y,double(datarefz),tardem.x,tardem.y,double(datatarz));
        rmsreg2(j)=sigma;
      if sum(size(z2out)~=size(refdem.z)) || sigma>10 || isnan(sigma) || max(abs(p(2:3))) >= 10 %control parameter
          warning(['coregistration failure',infile]); p=zeros(3,1);
          idd=[idd;j];
      else
        dx=p(2);dy=p(3); %z, x, y
        pg(j,1:3)=p;
% 	dzxyt=dzxyd(idreg(iref),:); %dx in the reg.txt file
% 	if abs(dx)+abs(dy)~=0
%         dx3=p(2:3);
        dX4Sg(j,1:3)=p+dx4;%p(2:3)+dx4(2:3); %
    
%         df=p+dx4-dzxyd(j,:); %validation z, x, y
      
        %write to a reg2.txt file
        if  0
        i=idregion(j);
        infile= strrep([demdir,'/',f{i}],'meta.txt','reg2.txt');
        demfile= strrep([f{i}],'meta.txt','dem.tif');
        fid10 = fopen(infile);
        fprintf(fid10,'DEM Filename: %s \n',demfile);
        fprintf(fid10,'Registration Dataset 1 Name: %s \n',tilefile);
        fprintf(fid10,'Registration Software: coregisterdems  \n');
        fprintf(fid10,'Translation Vector (dz,dx,dy)(m)= %d \n',[dX4S(j,1:3)] );
        fclose(fid10)
        end

	if flagplot==1

% The case when coregistration is not applied.
z2n = interp2(tardem.x' ,tardem.y,double(datatarz) ,refdem.x',refdem.y,'*linear');

[X,Y]=meshgrid(refdem.x,refdem.y);
[LATa,LONa]=polarstereo_inv(X,Y,[],[],70,-45);

figure;set(gcf,'Color','white');set(gca,'FontSize', 12);set(gcf, 'PaperPosition', [0.25 2.5 6 2.5]);
% imagesc(data0r.x*1e-3,data0r.y*1e-3,z2n-data0r.z);caxis([-5 5]);colorbar
surf(LONa,LATa,z2n-double(datarefz)); shading interp;
colorbar;colormap jet;view(0,90)
hl=xlabel('Longitude ($^{\circ}$)');
set(hl, 'Interpreter', 'latex');
hl=ylabel('Latitude ($^{\circ}$)');
set(hl, 'Interpreter', 'latex');
caxis([-5 5])
title([texttar,' - RefDEM; NO coregistration'])
ofile=[texttar,'mrefNocoreg'];
saveas(gcf,ofile,'fig')

figure;set(gcf,'Color','white');set(gca,'FontSize', 12);set(gcf, 'PaperPosition', [0.25 2.5 6 2.5]);
surf(LONa,LATa,z2out-double(datarefz)); shading interp;
% imagesc(data0r.x*1e-3,data0r.y*1e-3,z2out-data0r.z);caxis([-5 5]);colorbar
% title('2011/10/08-2013/05/26 DEM (m); After coregistration')
title([texttar,' - RefDEM; After coregistration'])
colorbar;colormap jet;view(0,90)
hl=xlabel('Longitude ($^{\circ}$)');
set(hl, 'Interpreter', 'latex');
hl=ylabel('Latitude ($^{\circ}$)');
set(hl, 'Interpreter', 'latex');
caxis([-8 8])
ofile=[texttar,'mrefwcoreg'];
saveas(gcf,ofile,'fig')

	end

      end
close all
	
end % if j

% remove strips that are not coregistered.
if 1  %hi; still use the data even if the coregistration is bad, dx4=0.
%remove data that are cropped too much.
%idd=idd2;
idd=[idd(:);idd2(:)];
idregion(idd)=[];%XYb(idd)=[];dzxy(idd)=[];
dzxyd(idd,:)=[];%rmsreg(idd)=[];
dX4Sg(idd,:)=[];rmsreg2(idd)=[];
end
fprintf ('\n Done with coregistration.\n ')
toc
clear data datar 
%data0r

%plot validation of coregistration
if 0
df=pg+repmat(dx4',length(idregion),1)-dzxyd;
text1={'z','x','y'};
figure;
hold all
set(gcf,'Color','white')
set(gcf, 'PaperPosition', [0.25 2.5 6 4.5]);
nt=length(idreg);
for j=1:3 %z x y
    subplot(1,3,j);hold all
set(gca,'FontSize', 12);
    plot(1:nt,dx4(j)*ones(nt,1),'ro:','linewidth',4,'Markersize',8)
plot(1:nt,pg(idreg,j),'go:','linewidth',4,'Markersize',8)
plot(1:nt,pg(idreg,j)+dx4(j),'bo-','linewidth',4,'Markersize',8)
plot(1:nt,dzxyd(idreg,j),'ko-','linewidth',4,'Markersize',8)
plot(1:nt,df(idreg,j),'mo-','linewidth',4,'Markersize',8)
plot(1:nt,rmsreg2(idreg),'g>-','linewidth',4,'Markersize',8)
title(text1(j));box on;
legend('tile reg.txt','coregistration','sum','strip reg.txt','difference','coregistration rms')
plot([1,11],df(idd,j),'cs','linewidth',4,'Markersize',12)
plot(1:nt,-pg(idreg,j)+dx4(j),'yo-','linewidth',4,'Markersize',8)
meanx=[mean(dx4(j)),mean(pg(idreg,j)), mean(pg(idreg,j)+dx4(j)), mean(dzxyd(idreg,j)), mean(df(idreg,j))]
end %if j
end % if 0

return
end
