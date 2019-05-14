function wm=mask2river(data)
% Given a water mask, get the high-precision river mask, exclude the lakes mask.
% data: input water mask, data.x, data.y,data.z; data.z (int8): -1 edge, 0 land, 1 water
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
    
    resr=abs(data.x(2)-data.x(1));
    narea=round(width2*width2/resr/resr);
%     nlb=round(width2/resr);

    %remove small clusters
    if 0
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
