
addpath(genpath('/home/dai.56/arcticdemapp/river/rivergithub2/'))

infile='/home/dai.56/data2/river/2018dec03/ortho/WV02_20160712063246_1030010059AE7E00_16JUL12063246-M1BS-500849898080_01_P002_u16ns3413.xml';
infile='/home/dai.56/data2/river/2018dec03/ortho/WV02_20110611224612_103001000B3D1A00_11JUN11224612-M1BS-052516089010_03_P009_u16ns3413.xml';
infile='/home/dai.56/data2/river/2018dec03/ortho/QB02_20100623213534_101001000BDD9E00_10JUN23213534-M1BS-500251602110_01_P013_u16ns3413.xml';

ndwisv=struct('x',[],'y',[],'z',[]);%Initialize
fprintf(['hi 1 ',num2str(int32(clock)),'.\n'])
whos

load wprob2.mat
clear Msum Msuma

        wm.z= jump;%water mask 
        wm.x=xout;wm.y=yout;   
        wm=cropmatrix(wm,wm.z);

fprintf(['hi 2 ',num2str(int32(clock)),'.\n'])
whos

%rangeov=[-2241000 -2220800 550200 539400];
rangeov=[-2242000 -2219000 532000 556000];
data=multispecmono(infile,wm,ndwisv,rangeov);
