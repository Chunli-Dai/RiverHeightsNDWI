function [clx2,cly2,S2]=interpcl(clx,cly,dc)
%Interpolation of centerline to any given interval along centerline distance.
% clx, cly input centerline nodes, in meter.
% dc input, requested centerline interval in meter
% clx2, cly2 output centerline nodes in meter

clx=clx(:);cly=cly(:);% be a column rather than an array.

S = [0; cumsum(sqrt(diff(clx(:)).^2+diff(cly(:)).^2))];
lencl=max(S);% total centerline length

p.dx=dc*1e-3; %km
p.x=0:p.dx:lencl*1e-3; %km
p.x=p.x*1e3;%m

%For each new node, find the two original nodes that includes this new
%node.
kpre=1;clx3=[];cly3=[];idk3=[];
clx2=zeros(length(p.x),1);cly2=clx2;
for i=1:length(p.x)
    
    %find k that S(k)<= p.x(i) <S(k+1)
    [k, ~] = find(S <= p.x(i), 1, 'last'); % find index of middle centerline node

    if ~isempty(k)
        if k< length(S) && k>=1 
            pt1=[clx(k) cly(k)];pt2=[clx(k+1) cly(k+1)];

            ds12=sqrt((pt1(1)-pt2(1))^2+(pt1(2)-pt2(2))^2);
            ds12p=p.x(i)-S(k); %has to be >=0;
            if ds12p <0;printf(['Warning: interpcl.m: ds12p < 0, find the wrong nodes! \n']);end

            ratio=ds12p/ds12;
            dx=(clx(k+1)-clx(k))*ratio;
            dy=(cly(k+1)-cly(k))*ratio;
            clx2(i)=clx(k)+dx;cly2(i)=cly(k)+dy;
        else %k=length(S)
            if S(k)==p.x(i)
                clx2(i)=clx(k);cly2(i)=cly(k);
            else
                printf(['Warning: interpcl.m: k is not within a reasonable range! \n'])
            end
        end
    
    else
        printf(['Warning: interpcl.m: k is empty! \n'])
    end
    
    %for validation
    clx3=[clx3(:);clx(kpre+1:k);clx2(i)];
    cly3=[cly3(:);cly(kpre+1:k);cly2(i)];  
    kpre=k;
    idk3=[idk3;length(clx3)]; %id of the interpolated nodes;
end

S2b = [0; cumsum(sqrt(diff(clx2(:)).^2+diff(cly2(:)).^2))];
S3 = [0; cumsum(sqrt(diff(clx3(:)).^2+diff(cly3(:)).^2))];
S2=p.x(:); %cumulative centerline distance of new nodes along original centerline.
fprintf(['\n Check the new centerline node, maximum centerline distance difference S3-p.x:',num2str(max(abs(S3(idk3)-p.x(:)))),'\n']);


return
end
