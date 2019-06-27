function [flagstats,id]=lowest(datarsv2,XYb2is)
% Find the lowest stage water masks
% Input: datarsv2 files of water masks;
%        XYb2is polygons of boundaries.
%Output flagstats: whether the lowest stage of all is found. 1 found. 0 not found.
%       id: the index of the lowest stage mask. 

        %Immune to Bug 2. Delete mask that are empty, all zeros, all -1 s.
        %save all zeros (could be dry river).
        n=length(datarsv2);
        idd=[];idkp=1:n;
        for i=1:n
            data=datarsv2(i);
            if isempty(data.z) || sum(sum(data.z~=-1))==0
                idd=[idd;i];
            end
        end
        XYb2is(idd)=[];datarsv2(idd)=[];
        idkp(idd)=[];
        
        n=length(datarsv2);
	if n==0
	    flagstats=0;id=[];
            return
	end
        range=zeros(n,4);
        for i=1:n
            Xb=XYb2is{i}(:,1); Yb=XYb2is{i}(:,2);
            range(i,1:4)=[min(Xb) max(Xb) min(Yb) max(Yb)];
        end
        
        resr=datarsv2(1).x(2)-datarsv2(1).x(1);
        
        rangeov=[max(range(:,1)),min(range(:,2)), max(range(:,3)),min(range(:,4))];
        if isempty(rangeov)
            flagstats=0;id=[];
            return
        
        else
            %stack all water masks to the overlapping area.
            ranget=round(rangeov/resr)*resr;
            tx=ranget(1):resr:ranget(2);ty=ranget(4):-resr:ranget(3);
            demg=-1*ones(length(ty),length(tx),n);
            for i=1:n
                data=datarsv2(i);
                tz = interp2(data.x,data.y,data.z,tx,ty','*nearest',-1);
                demg(:,:,i)=tz;
            end
            %Mkeep=any(demg~=-1,3);%matrix of pixels that are 0 or 1 (Non edges).
            Medge=any(demg==-1,3); Mkeep=~Medge; %For a pixel, if any of mask is -1, it is -1 void.
            if sum(Mkeep(:))==0;flagstats=0;id=[];return;end 
            
            area=zeros(n,1);
            for i=1:n
                tz=demg(:,:,i);
                tz(~Mkeep)=-1;
                area(i)=sum(sum(tz==1)); %water pixels;
            end
            [~,idsort]=sort(area);
            id=idsort(1);
            id=idkp(id);
            flagstats=1;
        end

end
