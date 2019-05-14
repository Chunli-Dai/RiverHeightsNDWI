function [c]=mask2centerline(data)

if 0 %cons: 1\ sensitive to the input width parameters -> let it be small; 2\ output are river segments instead of one lines;
% run the centerline detection method
%https://www.mathworks.com/matlabcentral/fileexchange/63030-improved-feature-centerline-central-line-extraction-and-cleaning?focused=7687377&tab=function
Max_River_width=10;  %an estimated river width; small is good
nedgelist=Improved_CentralLine(data.z,Max_River_width); %run the core centerline method
n=length(nedgelist);

if 0
    len=zeros(n,1);
    for i=1:n;len(i)=length(nedgelist{i});end
    M=len<=10; %remove branches < 400 m;
    nedgelist(M)=[];
end

n=length(nedgelist);

if 1 %connecting small branches
%Find unique segments
sep=zeros(n,4);
for i=1:n
    sep(i,1:2)=nedgelist{i}(1,:); %start
    sep(i,3:4)=nedgelist{i}(end,:); %end
    %consider sort 1 and 3 in case the order of start and end is flipped.
end
[C,idxu,idxc]=unique(sep,'rows','stable');
%https://www.mathworks.com/matlabcentral/answers/175086-finding-non-unique-values-in-an-array
% count unique values (use histc in <=R2014b)
[count, ~, idxcount] = histcounts(idxc,numel(idxu));
% Where is greater than one occurence
idxrp = find(count(idxcount)>1);
idnn=1:n;
%remove repeat segments 
idu=idnn(~ismember(idnn,idxrp)); % B(~ismember(B,A)) excluding A from B 

% select each repeated segments
[C,idxu,idxc]=unique(sep(idxrp,:),'rows','stable');
idrk=[];
for i=1:length(idxu)
    row=C(i,:);
    mat1=repmat(row,n,1);
    df=sum(sep-mat1,2);
    idr=find(df==0);
    
    %pick one of the branches. 
    % choose the one with shorter length.
    %Choose the one with larger width. -> hard to implement
    len=zeros(length(idr),1);
    for k=1:length(idr)
        j=idr(k);
        x1=(nedgelist{j}(:,2));y1=(nedgelist{j}(:,1));
        cl=[x1,y1];
        dists = sqrt(diff(cl(:,1)).^2+diff(cl(:,2)).^2); % distances between each centerline node
        cumdists = [0; cumsum(dists)]; % cumulative distance along centerline
        len(k) = cumdists(end);

%         spacing = 2;
%         [Wm, SWm] = width_from_mask(data.z, cl, spacing);

    end
    [~,k]=min(len);
    idrk=[idrk(:);idr(k)];
end

id=[idu(:);idrk(:)];
x=[];y=[];
%connect the segments
for k=1:length(id)
    i=id(k);
    x1=data.x(nedgelist{i}(:,2));y1=data.y(nedgelist{i}(:,1));
x=[x(:);NaN;x1(:)];y=[y(:);NaN;y1(:)];
% hold on;plot(x,y,'.-')
end
[xm, ym] = polymerge(x, y);
hold on;plot(xm,ym,'go-')
else
    len=zeros(n,1);
    for i=1:n;len(i)=length(nedgelist{i});end
    [~,i]=max(len);
    xm=data.x(nedgelist{i}(:,2));ym=data.y(nedgelist{i}(:,1));
end

[lat,lon]=polarstereo_inv(xm,ym,[],[],70,-45);
c.X=lon;c.Y=lat; %.X (longitude), .Y (latitude).

else %Cons: 1\sometimes it is not working. imregionalmin deleted too many points. 2\ need es parameter.
Wn = 10; % nominal width
es = 'SW'; % exit sides
es = 'WE'; % START stream at W, End stream at E.
es ='EW';
es = getes(data);
plotornot=0; %1 plot;0 not plot
% 3a. Compute the centerline
[cl, Icl] = centerline_from_mask(data.z,es,Wn,plotornot);
x=data.x(cl(:,1));y=data.y(cl(:,2));
[lat,lon]=polarstereo_inv(x,y,[],[],70,-45);
c.X=lon;c.Y=lat; %.X (longitude), .Y (latitude).
end

hold on;plot(lon,lat,'k-')

return
end

