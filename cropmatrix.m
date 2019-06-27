function [wm]=cropmatrix(data,M)
%crop matrix using the M matrix (1 valid data, 0 non valid data); 
%data.z and M should have the same size;
jump=M;
xout=data.x;yout=data.y;
[X,Y]=meshgrid(xout,yout);
xmin=min(X(jump==1));xmax=max(X(jump==1));
ymin=min(Y(jump==1));ymax=max(Y(jump==1));
idx=xout>=xmin&xout<=xmax;idy=yout>=ymin&yout<=ymax;
wm.x=xout(idx);wm.y=yout(idy);wm.z=data.z(idy,idx);

return
end