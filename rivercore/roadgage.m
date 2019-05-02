function groad=roadgage(rang0,odir)
constant

[status , cmdout ]=system(['find ',codedir,' -name groad.shp']);
shpname=deblank(cmdout);
x0=[rang0(1) rang0(2) rang0(2) rang0(1) rang0(1) ];y0=[rang0(4) rang0(4) rang0(3) rang0(3) rang0(4) ];
[lat0,lon0]=polarstereo_inv(x0,y0,[], [],70,-45);
bb = geoshape(lat0,lon0,'Geometry','Polygon');
tileshape=sprintf('%s/gageareabox.shp',odir);
tilecoastname=sprintf('%s/gagearearoads.shp',odir);
shapewrite(bb,tileshape);
if ~exist(tilecoastname,'file')
% system(['rm ',tilecoastname])
system(['time ogr2ogr -overwrite -clipsrc ',tileshape,' ',tilecoastname,' ',shpname]);
%ogr2ogr -overwrite -clipsrc tile.shp tilegshhs.shp GSHHS/GSHHS_f_L1.shp
end
S = shaperead(tilecoastname);
cnt=length(S); %figure;mapshow(S);

%Expand the road for 600 meters for control surfaces.
widthroad=100;
xw=rang0(1):resrc:rang0(2);yw=rang0(4):-resrc:rang0(3); % a priori water mask
nwx=length(xw);nwy=length(yw);
smg=false(nwy,nwx);%water mask from a priori coastline shapefiles
for j=1:cnt
    [sx,sy]=polarstereo_fwd(S(j).Y,S(j).X,[], [],70,-45);
    P=[sx(:),sy(:)];
    polyout = polybuffer(P,'lines',widthroad);
%     figure;plot(polyout)
    sx=polyout.Vertices(:,1);sy=polyout.Vertices(:,2);
    dx=resrc;
    idx=round((sx-rang0(1))/dx)+1;idy=round((sy-rang0(4))/(-dx))+1;
    %Oct 10, 2018:fix bug 14; separate polygons using the NaN;
    M=[idx(:),idy(:)];
    idx = any(isnan(M),2);
    idy = 1+cumsum(idx);
    idz = 1:size(M,1);
    C = accumarray(idy(~idx),idz(~idx),[],@(r){M(r,:)});
    for k=1:length(C)
        idx=C{k}(:,1);idy=C{k}(:,2);   
        sm=poly2mask(idx,idy,nwy,nwx); % fast, apply to each polygon one by one.
        smg=smg|sm;
    end
end
groad.x=xw;groad.y=yw;%initialize
groad.z=smg;clear smg; 
if(sum(groad.z(:))==0);fprintf('This area contain no road from groads database!');end % 
return
end