function [c,widave]=getcl(rang0,loneq,lateq)
%Given the region boundary [xmin xmax ymin ymax ] in polar stereographic coordicates
% To collect the river centerline from Elizabeth H. Altenau's database
% loneq,lateq is target gage location
% output: c.X c.Y: longitude and latitude of output rivercenterline (going
% uphill).
% widave: average width in m.

%input parameters
if 0
bbox=[-149.2 64.5; -149  64.6];
bbox=[-149.2-1 64.5-1; -149+1  64.6+1];
loneq=-149.1;lateq=64.55;
end
rang0in=rang0;
exb=15e3;%expand the box by 15km for the given centerline.
rang0=[rang0in(1)-exb rang0in(2)+exb rang0in(3)-exb rang0in(4)+exb];

%rang0=[xeq-exb xeq+exb yeq-exb yeq+exb ];
x0=[rang0(1) rang0(2) rang0(2) rang0(1) rang0(1) ];y0=[rang0(4) rang0(4) rang0(3) rang0(3) rang0(4) ];
[lat,lon]=polarstereo_inv(x0,y0,[],[],70,-45);
bbox=[min(lon) min(lat); max(lon) max(lat)];

range=[bbox(1,1) bbox(2,1) bbox(1,2) bbox(2,2)];
x0=[range(1) range(2) range(2) range(1) range(1) ];y0=[range(4) range(4) range(3) range(3) range(4) ];
% figure;hold all;plot(x0,y0,'-','linewidth',3)
% addpath(genpath(['/Users/chunlidai/surge/home/dai.56/arcticdemapp/river/rivergithub2v2/']))
% addpath(genpath(['/home/dai.56/arcticdemapp/river/rivergithub2v2/']))
S2=shaperead('NA_Merge.shp','BoundingBox',bbox);%10 minutes
% figure;mapshow(S2);
% figure;plot(S2(:).X,S2(:).Y,'k-');
ns=length(S2);

save t1.mat -v7.3 

if ns==0
    fprintf(['\n No centerline in the region.']);
c.X=[];c.Y=[];widave=[];c.widave=[];
return
end

sx=zeros(ns,2);
for j=1:ns
sx(j,:)=[S2(j).X,S2(j).Y];
end

idsc=zeros(ns,2);
hwf=zeros(ns,3);
for j=1:ns
idsc(j,:)=[S2(j).segID,S2(j).cl_id];
hwf(j,:)=[S2(j).p_height,S2(j).p_width,S2(j).flowacc]; %m m km^2
basin(j)=S2(j).basin_code;
end
% [~,idcl]=sort(idsx(:,2));

% figure;plot(sx(:,1),sx(:,2),'.-');
ds=sqrt(diff(sx(:,1)).^2+diff(sx(:,2)).^2);
ds=median(ds);%median point distance in degrees.
%Find the segment that's cloest to the gage
dist=sqrt((sx(:,1)-loneq).^2+(sx(:,2)-lateq).^2);
[mindist,idg]=min(dist);
mindistp=mindist*(pi/180)*6371e3;
if isempty(mindistp) || isnan(mindistp)
fprintf(['getcl.m error, mindistp=',num2str(mindistp),' m.']);
c.X=[];c.Y=[];widave=[];c.widave=[];
return
elseif mindistp > 1e3 
fprintf(['Gage is far away from the centerline:',num2str(mindistp),' m.']);
c.X=[];c.Y=[];widave=[];c.widave=[];
return
end

% segu=unique(idsc(:,1));
% consider seg i been divided into two segments by boundingbox;
id1=idsc(2:end,2)-idsc(1:end-1,2);
idd=find(abs(id1)>1);%potential risk: two segments, id differred by 1.
segui=[1;idd(:)+1];%start id of each segment.
id2=id1;id2(abs(id2)==1)=0;
id2(abs(id2)>1)=1;
idsc(:,3)=cumsum([1;id2(:)]);
% segu=idsc(segui,1);
% seg0=idsc(idg,1);
segu=1:length(segui);%rename the segments id segID
seg0=idsc(idg,3);
iseg0=find(segu==seg0);

%separate each segments;
datarsv(length(segu))=struct('X',[],'Y',[],'segID',[],'cl_id',[],'p_heigth',[],'p_width',[],'flowacc',[]);

%Find the end points of each segment;
pts=zeros(length(segu),2);pte=zeros(length(segu),2);
for j=1:length(segu)
    M=(idsc(:,3)==segu(j));

    datarsv(j).X=sx(M,1);
    datarsv(j).Y=sx(M,2);
    datarsv(j).segID=idsc(M,1);
    datarsv(j).cl_id=idsc(M,2);
    datarsv(j).p_height=hwf(M,1);
    datarsv(j).p_width=hwf(M,2);
    datarsv(j).flowacc=hwf(M,3);
    
    pts(j,:)=[datarsv(j).X(1),datarsv(j).Y(1)];
    pte(j,:)=[datarsv(j).X(end),datarsv(j).Y(end)];
end

%coordinates of start and end points;
figure
for j=1:length(segu)
    hold all;plot(datarsv(j).X,datarsv(j).Y,'.-','linewidth',1);
end
legend(['seg ',num2str(1)],num2str(2),num2str(3),num2str(4),num2str(5))
hold all;plot([sx(idg,1),loneq],[sx(idg,2),lateq],'g>-');
hold all;plot(pts(:,1),pts(:,2),'g*');
hold all;plot(pte(:,1),pte(:,2),'b*');
hold all;plot(x0,y0,'-','linewidth',3)

%Find the connect points from the START point of seg0 .
iter=0;isegcs1sv=[];ptcs1sv=[];
while true
iter=iter+1;
%given segment
if iter ==1
    seggvn=seg0;iseggvn=iseg0;
else
    iseggvn=isegcs1sv(iter-1);
    seggvn=segu(iseggvn);
end
ids=find(segu==seggvn);

%START point of the given segment.
if iter ==1
pt0=[datarsv(ids).X(1) datarsv(ids).Y(1)];
else
    id1=[1;length(datarsv(ids).X)];
    ptid=ptcs1sv(iter-1);%point id of previous connecting point.
    id1=id1(id1~=ptid);%exclude ptid.
    if length(id1)~=1;fprintf('warning: wrong start point.');end
    pt0=[datarsv(ids).X(id1) datarsv(ids).Y(id1)];
end

%start points of all segments, and the distance to pt0.
dists=sqrt((pt0(:,1)-pts(:,1)).^2+(pt0(:,2)-pts(:,2)).^2);

% find the segments that have start point close to pt0.
ids=find(dists<2*ds);
isegcs1=ids(ids~=iseggvn); %exclude the given segment itself.
idcs1=ones(size(isegcs1));

%find segments that have end point close to pt0.
dists=sqrt((pt0(:,1)-pte(:,1)).^2+(pt0(:,2)-pte(:,2)).^2);
ide=find(dists<2*ds);
%connecting segments candidates ;%excluding the given segment;
isegcs1e=ide(ide~=iseggvn);
idcs1e=zeros(size(isegcs1e));
for j=1:length(isegcs1e)
    js=isegcs1e(j);
idcs1e(j)=length(datarsv(js).cl_id);%end point id.
end
isegcs1=[isegcs1(:);isegcs1e(:)];%seg id
idcs1=[idcs1(:);idcs1e(:)]; %point id

if isempty(isegcs1)
    fprintf(['\n Searching the connecting point from START of seg0, total segments:',num2str(iter)])
    break
end

%compare the drainage area and pick the largest one.
flowacc_cs1=zeros(size(isegcs1));
for j=1:length(isegcs1)
    js=isegcs1(j);ptid=idcs1(j);
    flowacc_cs1(j)=datarsv(js).flowacc(ptid);
end
    
[~,idpick]=max(flowacc_cs1);
%selected segment id and point id
isegcs1sv(iter)=isegcs1(idpick);
ptcs1sv(iter)=idcs1(idpick);

if iter > 5 ;fprintf('Iteration of finding connecting segments larger than 5.');break ;end
%Find the connect points at the start point of segcs1.

end

%Find the connect points from the END point of seg0 .
iter=0;isegce1sv=[];ptce1sv=[];
while true
iter=iter+1;
%given segments
if iter ==1
    seggvn=seg0;iseggvn=iseg0;
else
    iseggvn=isegce1sv(iter-1);
    seggvn=segu(iseggvn);
end
ids=find(segu==seggvn);

%END points of the given segment
if iter ==1
pt0=[datarsv(ids).X(end) datarsv(ids).Y(end)];
else
    id1=[1;length(datarsv(ids).X)];
    ptid=ptce1sv(iter-1);%point id of previous connecting point.
    id1=id1(id1~=ptid);%exclude ptid.
    if length(id1)~=1;fprintf('warning: wrong end point.');end
    pt0=[datarsv(ids).X(id1) datarsv(ids).Y(id1)];
end

%start points of all segments, and the distance to pt0.
dists=sqrt((pt0(:,1)-pts(:,1)).^2+(pt0(:,2)-pts(:,2)).^2);

% find the segments that have start point close to pt0.
ids=find(dists<2*ds);
isegcs1=ids(ids~=iseggvn); %exclude the given segment itself.
idcs1=ones(size(isegcs1));

%find segments that have end point close to pt0.
dists=sqrt((pt0(:,1)-pte(:,1)).^2+(pt0(:,2)-pte(:,2)).^2);
ide=find(dists<2*ds);
%connecting segments candidates ;%excluding the given segment;
isegcs1e=ide(ide~=iseggvn);
idcs1e=zeros(size(isegcs1e));
for j=1:length(isegcs1e)
    js=isegcs1e(j);
idcs1e(j)=length(datarsv(js).cl_id);%end point id.
end
isegcs1=[isegcs1(:);isegcs1e(:)];%seg id
idcs1=[idcs1(:);idcs1e(:)]; %point id

if isempty(isegcs1)
    fprintf(['\n Searching the connecting point from END of seg0, total segments:',num2str(iter)])
    break
end

%compare the drainage area and pick the largest one.
flowacc_cs1=zeros(size(isegcs1));
for j=1:length(isegcs1)
    js=isegcs1(j);ptid=idcs1(j);
    flowacc_cs1(j)=datarsv(js).flowacc(ptid);
end
    
[~,idpick]=max(flowacc_cs1);
%selected segment id and point id
isegce1sv(iter)=isegcs1(idpick);
ptce1sv(iter)=idcs1(idpick);

if iter > 5 ;fprintf('Iteration of finding connecting segments larger than 5.');break ;end
%Find the connect points at the start point of segcs1.

end

%connect all segments.
ns=length(isegcs1sv);ne=length(isegce1sv);
nseg=ns+1+ne;
segc=zeros(nseg,1);
Xc=[];Yc=[];hwfc=[];%connected centerline coordinates
for iseg=1:ns
    iter=ns-iseg+1;
    ids=isegcs1sv(iter);
    segc(iseg)=ids;
    
    %start and end points
    id1=[1;length(datarsv(ids).X)];
    ptid=ptcs1sv(iter);%point id of the connecting point: End point of this segment when coonecting.
    id1=id1(id1~=ptid);%exclude ptid.
    dn=(ptid-id1)/abs(ptid-id1);%1 or -1
    id=id1:dn:ptid;
    pt0=[datarsv(ids).X(id) datarsv(ids).Y(id)];
    hwfcj=[datarsv(ids).p_height(id) datarsv(ids).p_width(id) datarsv(ids).flowacc(id)];
    
    Xc=[Xc;pt0(:,1);];Yc=[Yc;pt0(:,2);];
    hwfc=[hwfc;hwfcj];
end

%seg0
ids=iseg0;
pt0=[datarsv(ids).X  datarsv(ids).Y ];
Xc=[Xc;pt0(:,1);];Yc=[Yc;pt0(:,2);];
segc(ns+1)=iseg0;
hwfcj=[datarsv(ids).p_height(:) datarsv(ids).p_width(:) datarsv(ids).flowacc(:)];
hwfc=[hwfc;hwfcj];

%
for iseg=1:ne
    iter=iseg;
    ids=isegce1sv(iter);
    segc(ns+1+iseg)=ids;
    
    %start and end points
    id1=[1;length(datarsv(ids).X)];
    ptid=ptce1sv(iter);%point id of the connecting point: START point of this segment when coonecting.
    id1=id1(id1~=ptid);%exclude ptid.
    dn=(id1-ptid)/abs(ptid-id1);%1 or -1
    id=ptid:dn:id1;
    pt0=[datarsv(ids).X(id) datarsv(ids).Y(id)];
    hwfcj=[datarsv(ids).p_height(id) datarsv(ids).p_width(id) datarsv(ids).flowacc(id)];
    
    Xc=[Xc;pt0(:,1);];Yc=[Yc;pt0(:,2);];
    hwfc=[hwfc;hwfcj];

end
%
c.X=Xc;c.Y=Yc; %longitude X, latitude Y
widave=mean(hwfc(:,2));
c.widave=widave;

%adjust the centerline to go uphill for ProcessTananaFairbanks.m.
pth=[hwfc(1,1); hwfc(end,1)];
if pth(2)<pth(1) %downhill
    fprintf(['\n Flip centerline to go uphill. Heights of two ends:',num2str(pth(:)')]);
    c.X=flip(c.X);
    c.Y=flip(c.Y);
    hwfc=flip(hwfc);
end
save clsv2.mat c

hold on;plot(c.X,c.Y,'k-',c.X(1),c.Y(1),'k>','linewidth',3,'markersize',8);
box on
xlabel('Longitude');ylabel('Latitude')
saveas(gcf,'ElizabethCenterline','fig')

figure;
subplot(2,1,1);hold all;
plot(hwfc(:,1),'.-'); %height
plot(hwfc(:,2),'.-') %width
legend('height','width')
box on
subplot(2,1,2)
plot(hwfc(:,3),'.-') %flowacc
legend('flowacc (km^2)')
box on
% legend('height','width','flowacc')
set(gcf,'Color','white');set(gca,'FontSize', 12);set(gcf, 'PaperPosition', [0.3 0 3.8 3]); set(gcf, 'PaperSize', [4 3]);
xlabel('River ceterline node number')

 save t2.mat -v7.3

end





