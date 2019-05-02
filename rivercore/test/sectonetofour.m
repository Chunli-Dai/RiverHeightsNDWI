addpath(genpath(['../rivergithub2/']))

sectionname='sagExtent.shp';
S = shaperead(sectionname);
cnt=length(S);
figure;plot(S.X,S.Y)

j=1;
[sx,sy]=polarstereo_fwd(S(j).Y,S(j).X,[], [],70,-45);
figure;plot(sx*1e-3,sy*1e-3)

p=[-2253, 554; -2235, 550; -2217, 545; -2200, 537; ]*1e3;
hold all;plot(p(:,1)*1e-3,p(:,2)*1e-3,'g>')

exb=10e3;
for j=1:4
    xeq=p(j,1);yeq=p(j,2);
rang0=[xeq-exb xeq+exb yeq-exb yeq+exb ];
% rang0=[-3467 -3455 110 124 ]*1e3;
% rang0=[-3288 -3281 352 358]*1e3;
x0=[rang0(1) rang0(2) rang0(2) rang0(1) rang0(1) ];y0=[rang0(4) rang0(4) rang0(3) rang0(3) rang0(4) ];
plot(x0*1e-3,y0*1e-3,'g-','linewidth',3)
end
xlabel('x (km)');ylabel('y (km)')

[lat,lon]=polarstereo_inv(p(:,1),p(:,2),[], [],70,-45);
[lon,lat]

ans =

 -148.8146   68.8173
 -148.8249   68.9819
 -148.8110   69.1488
 -148.7172   69.3135

