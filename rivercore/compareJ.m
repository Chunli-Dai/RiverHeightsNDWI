macdir='/Users/chunli/surge/';
infile=[macdir,'/data3/ArcticDEM/region_34_alaska_north/tif_results/2m/GE01_20130526_1050410002A6FC00_10504100029EAE00_053537405030_01_P001_053537402030_01_P001_2_dem.tif'];

[M,ratio]=maskentropy(infile);
frozen=0

M=(J<3);
M = bwareaopen(M, 1000*200);%20121011
figure;imagesc(M)
Modj=M;
Modfil = bwareaopen(~Modj, 1000*5);
Modfil=~Modfil;
figure;imagesc(Modfil)
Md1 = imdilate(Modfil, ones(3));
M=logical(Md1-Modfil);
% M=M&Mden&~(imdilate(or == 0,ones(round(9*2*8/resr))));
M=M&~(imdilate(or == 0,ones(round(9*2*8/resr))));
figure;imagesc(M)
or=readGeotiff(orFile);
[X,Y]=meshgrid(or.x,or.y);
p=zeros(3,1)
[LATa,LONa]=polarstereo_inv(X- p(2),Y- p(3),[],[],70,-45);
output=[LATa(M),LONa(M)];
save('2013MayCoastJ3.dat','output','-ascii')

figure;imagesc(or.x,or.y,J);colorbar;caxis([0 5])
hold on;plot(X(M),Y(M),'r.')

if 0
figure;
set(gcf,'Color','white')
set(gca,'FontSize', 12);
set(gcf, 'PaperPosition', [0.25 2.5 6 2.5]);
hold all;
surf(LONa,LATa,J);
shading interp;
colorbar;colormap jet;view(0,90)
caxis([0 6])
hold on;plot3(LONa(M),LATa(M),1e9*ones(size(LATa(M))),'r.','Markersize',8)
view(0,90)
box on
hl=xlabel('Longitude ($^{\circ}$)');
set(hl, 'Interpreter', 'latex');
hl=ylabel('Latitude ($^{\circ}$)');
set(hl, 'Interpreter', 'latex');
end


