function wm=mask2river(data)
% Given a water mask, get the high-precision river mask, exclude the lakes mask.
% data: input water mask, data.x, data.y,data.z; data.z (int8): -1 edge, 0 land, 1 water
%       input water mask's small clusters are removed before this function.
% wm: output river mask, wm.x, wm.y,wm.z (int8)
% Refers to coastline/codec2/icluster.m

%
    width2=100; %river width m
    dataorg=data;

    %Plot water mask of 50% probability.
    [X,Y]=meshgrid(data.x,data.y);
    xmin=min(X(data.z==1));xmax=max(X(data.z==1));
    ymin=min(Y(data.z==1));ymax=max(Y(data.z==1));
    idx=data.x>=xmin&data.x<=xmax;idy=data.y>=ymin&data.y<=ymax;
    wm.x=data.x(idx);wm.y=data.y(idy);BW=data.z(idy,idx);
    Medge=(BW==-1);
    BW(BW==-1)=0;

    %OUTPUT 
    [m,n]=size(BW);
    if ~( min([m,n])>2 )
      %Attention: could be a dry out river.
      warning(['mask2river: water cluster all removed.'])
      wm=dataorg;
      wm.z=zeros(size(dataorg.z),'int8'); %could be a dry out river.
      wm.z(dataorg.z==-1)=int8(-1); %keep the edge for lowest.m
      %wm.z=-1*ones(size(dataorg.z),'int8');  %must be -1, no valid data;affect multispecstrip.m and lowest.m
      return
    end
    
    resr=abs(data.x(2)-data.x(1));
    narea=round(width2*width2/resr/resr);
%     nlb=round(width2/resr);

    %remove small clusters
    if 0 %input water mask already removed small clusters
    Modj= bwareaopen(BW, narea);
    BW=Modj;
    end
        
    % Test distinguishing rectangle from squares / circles.
    if 0
        %area/perimeter^2; L=nW;
        n=1:100;
        ap2=n./(n+1).^2/4;ap2c=1/4/pi;
        figure;plot(n,ap2,'.-');
        hold on;plot(n,ap2c*ones(size(n)),'r-')
        xlabel('Length/width Ratio')
        ylabel('Area/Perimeter Square Ratio')
        %n=5;ap2=0.035;
    end % if 0
%     ap2thres=0.035;%L=5W;

    %old method
    if ~exist('clsv2.mat','file') 
    fprintf(['\n River centerline file clsv2.mat not found; Using the old algorithm in mask2river.m. ']);
    
    CC = bwconncomp(BW);
    BW2=BW;BW2(:)=0;  
    
    data.x=wm.x;data.y=wm.y;
    [X,Y]=meshgrid(data.x,data.y);

    resrc=40.;dsr=resr/resrc;
    nareaf=round(width2*width2/resrc/resrc); %size to fill holes

    for k=1:CC.NumObjects
% k
            %check if this cluster is a long shape rather than square. 
            BW3=BW;BW3(:)=0;
            BW3(CC.PixelIdxList{k})=BW(CC.PixelIdxList{k});
            
            data.z=BW3;
            xmin=min(X(data.z==1));xmax=max(X(data.z==1));
            ymin=min(Y(data.z==1));ymax=max(Y(data.z==1));
            idx=data.x>=xmin&data.x<=xmax;idy=data.y>=ymin&data.y<=ymax;

            Modj=BW3(idy,idx);  
            
            %when checking the shape of the water body, use 40m resolution
            Modj=imresize(Modj,dsr);
            
            Modfil = bwareaopen(~Modj, nareaf); %fill small areas: 1e4*4m^2
            Modfil=~Modfil;
            
            %expand BW3 for geting the edges at boundary
            [n1,m1]=size(Modfil);
            BW3=zeros(n1+2,m1+2);BW3(2:end-1,2:end-1)=Modfil;
            
            Md1 = imerode(BW3, ones(3));
            M=logical(-Md1+BW3);

            area=sum(sum(BW3));peri=sum(M(:));
            ratio=area/peri^2; % area / perimeter^2
            
            a=4*ratio;b=8*ratio-1;c=4*ratio;
            rLW=(-b+sqrt(b^2-4*a*c))/(2*a); %estimated ratio: Length/width
	    %Not working well, due to the increased perimeter when points are scattered.
            if ~isreal(rLW);continue;end
            
            if 0 %use advanced method; time consuming
%             Max_River_width=sqrt(area/rLW);  %an estimated river width; small is good
%             nedgelist=Improved_CentralLine(BW3,Max_River_width); %run the core centerline method
            
            Wn = round(sqrt(area/rLW)); % nominal width
            es = 'SW'; % exit sides
            plotornot=1;
            % 3a. Compute the centerline
            try 
            [cl, Icl] = centerline_from_mask(BW3,es,Wn,plotornot);
		%Not working well if the centerline is not correctly retrieved.
            catch e
                %tif provided by pgc might not be compatable with readGeotiff.
                fprintf('There was an error! The message was:\n%s',e.message);
                %fprintf('Read tif file using importdata instead of readGeotiff.')
                continue
            end

            % Parameterizing centerline by along-stream distance
            if isempty(cl)
                continue
            end
            dists = sqrt(diff(cl(:,1)).^2+diff(cl(:,2)).^2); % distances between each centerline node
            cumdists = [0; cumsum(dists)]; % cumulative distance along centerline
            len = cumdists(end); % total centerline length
            width=area/len;rLW=len/width;
            title(['Length/Width Ratio:',num2str(rLW)]);axis equal
            end % if use advanced method
            
%             if ratio<ap2thres % squares /circles
            if rLW >40% 10 %30 
            %40 to remove the long shape water bodies near Tanana gage
                                    
            BW2(CC.PixelIdxList{k})=BW(CC.PixelIdxList{k});
            
            if 0
%                 tt.x=data.x(idx);tt.y=data.y(idy);
                tt.x=imresize(data.x(idx),dsr);tt.y=imresize(data.y(idy),dsr);
                figure;set(gcf,'Color','white');set(gca,'FontSize', 12);set(gcf, 'PaperPosition', [0.25 2.5 6 4]);hold all;
                imagesc(tt.x*1e-3,tt.y*1e-3,Modj);title(['Length/Width Ratio:',num2str(rLW)]);axis equal
                xlabel('x (km)');ylabel('y (km)')
                figure;set(gcf,'Color','white');set(gca,'FontSize', 12);set(gcf, 'PaperPosition', [0.25 2.5 6 4]);hold all;
                imagesc(BW3);title(['Length/Width Ratio:',num2str(rLW)]);axis equal
                [Xs,Ys]=meshgrid(1:(m1+2),1:(n1+2));
                hold on;plot(Xs(M),Ys(M),'r.')
                
                figure(1);imagesc(wm.x*1e-3,wm.y*1e-3,BW2);title(['Collected BW2;Length/Width Ratio:',num2str(rLW)]);axis equal
            end
            end
            
    end
    BW2(Medge)=-1;
    wm.z=BW2;

    else
%  Try SWOT algorithm.
%  Need the a priori river centerline.
    %given wm.x, wm.y, BW (1 water 0 non-water) 
        load('clsv2.mat','c')
%       fprintf(['\n River centerline file clsv2.mat found; Using SWOT algorithm in mask2river.m. ']);

        [clx,cly]=polarstereo_fwd(c.Y,c.X,[], [],70,-45);
        if isfield(c, 'widave')
        widave=c.widave;
        else
            % use a fixed max_distance
            widave=20;
        end
        
        %river buffer; refers to prepwm.m
        %get bufferzone along river centerline, to remove other tributaries.
        width80=widave;
        resx=mean(data.x(2:end)-data.x(1:end-1));resy=mean(data.y(2:end)-data.y(1:end-1));
        resr=mean([abs(resx),abs(resy)]);
        buf=zeros(size(BW));
        %polar stereographic coordinates to image coordinates.
        [ny,nx]=size(BW);
        clear cl
        cl(:,1)=round((clx-wm.x(1))/resx)+1;
        cl(:,2)=round((cly-wm.y(1))/resy)+1;
        M=cl(:,1)>=1&cl(:,1)<=nx&cl(:,2)>=1&cl(:,2)<=ny;
        cl(~M,:)=[];
%       buf(cl(:,2),cl(:,1))=1;
        for j=1:length(cl(:,1))
        buf(cl(j,2),cl(j,1))=1;
        end
        widpix=round(width80/resr);
        ncl=3; % expand along centerline by ncl times; try 10, 5 3
        widpix2=round(width80/resr*ncl);
        clbuf= imdilate(buf, ones(widpix2*2)); % width expansion
              
        if 0 %SWOT algorithm
        %step 1 : calculating distances of pixels to centerline nodes, and assign the closest node.
        %reduce centerline resolution to 200 m node interval
        S = [0; cumsum(sqrt(diff(clx(:)).^2+diff(cly(:)).^2))];
        rescl=nanmean(S(2:end)-S(1:end-1));
        nodeint=200;
        if rescl < nodeint %
        fprintf(['\n Reduce the centerline node interval to 200 m for SWOT algorithm!'])
        nsr=round(nodeint/rescl);
        clx=clx(1:nsr:end);cly=cly(1:nsr:end);
        end
        S = [0; cumsum(sqrt(diff(clx(:)).^2+diff(cly(:)).^2))];
        rescl=nanmean(S(2:end)-S(1:end-1));
        
        max_distance=max(widave,nodeint+1);
        [X,Y]=meshgrid(wm.x,wm.y);
        p1=[X(logical(BW)), Y(logical(BW))];
        p2=[clx(:),cly(:)];
        tic;Z=pdist2(p1,p2);toc %Error using pdist2mex Requested 6395498x15826 (754.1GB) array exceeds maximum array size preference. e.g. wprob2 nsubx=10001;
                                %Error fixed if reducing centerline node interval to 200m.
        [Zcl,idcl]=min(Z,[],2);
        %plot3(p1(:,1),p1(:,2),idcl)
        
        % step 2a: dist <= max_distance
        M=Zcl<=max_distance;
%         figure;plot(p1(M,1)*1e-3,p1(M,2)*1e-3,'.')
        %get the pixel matrix
        [ny,nx]=size(BW);
        BW2a=zeros(size(BW));
        x0=wm.x(1);dx=(wm.x(end)-wm.x(1))/(length(wm.x)-1);
        y0=wm.y(1);dy=(wm.y(end)-wm.y(1))/(length(wm.y)-1);
        idx=round((p1(M,1)-x0)/dx)+1;
        idy=round((p1(M,2)-y0)/dy)+1;
        BW2a(idy(:)+(idx(:)-1)*ny)=1;
% 	save t1.mat -v7.3
        
        %step 2b find the cluster that's contiguous with the dominant
        %segment.
        % find the dominant segment
        CC = bwconncomp(BW2a);
        segsize=zeros(CC.NumObjects,1);
        for k=1:CC.NumObjects
            segsize(k)=length(CC.PixelIdxList{k});
        end
        [~,k]=max(segsize);
        BW3=zeros(size(BW));
        BW3(CC.PixelIdxList{k})=BW2a(CC.PixelIdxList{k});
        BWdom=BW3;
        
        %find the contiguous cluster
        CC = bwconncomp(BW);
        BW4=false(size(BW));
        for k=1:CC.NumObjects
            BW3=zeros(size(BW));
            BW3(CC.PixelIdxList{k})=BW(CC.PixelIdxList{k});
            
            %check if this cluster is contiguous with the dominant segment.
            overlap=sum(sum(BW3&BWdom));
            
            if overlap >=1 
                BW4=BW4|BW3;
                ksv=k;
                break
            end
        end
        BW3=BW4;
        
        else %Modified SWOT algorithm
            %Keep any cluster that contains the given river centerline.
          %find the contiguous cluster
        CC = bwconncomp(BW);
        BW4=false(size(BW));
        for k=1:CC.NumObjects
            BW3=zeros(size(BW));
            BW3(CC.PixelIdxList{k})=BW(CC.PixelIdxList{k});
            
            %check if this cluster is contiguous with the dominant segment.
            overlap=sum(sum(BW3&buf));
            
            if overlap >=1 
                BW4=BW4|BW3;
                ksv=k;
%                 break
            end
        end
        BW3=BW4;            
        end
        

        BW2=BW3&clbuf;
        BW2(Medge)=-1;
        wm.z=BW2;
                
    end
	
    %OUTPUT 
    [m,n]=size(wm.z);
    if min([m,n])>2
    %recover the size of wm to be the same as dataorg
%   tz = interp2(wm.x,wm.y,double(wm.z),dataorg.x,dataorg.y','*nearest',-1);%recover its size.
    tz = interp2(wm.x,wm.y,int8(wm.z),dataorg.x,dataorg.y','*nearest',-1);%recover its size. %save space
    wm=dataorg;wm.z=tz; 
    else
      %Attention: could be a dry out river.
      warning(['mask2river: water cluster all removed.'])
      wm=dataorg;
      wm.z=zeros(size(dataorg.z),'int8'); %could be a dry out river.
      wm.z(dataorg.z==-1)=int8(-1); %keep the edge for lowest.m
      %wm.z=-1*ones(size(dataorg.z),'int8');  %must be -1, no valid data;affect multispecstrip.m and lowest.m
    end
    
return
end
