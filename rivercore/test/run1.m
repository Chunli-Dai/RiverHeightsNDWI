addpath(genpath(['../rivergithub2/']))
load clsv1.mat
        c.X=flip(c.X);
        c.Y=flip(c.Y);
odir='gage32/';
odir='gage32/gage32sag/';
loneq=-148.8177778;lateq= 69.0158333;

M=isnan(c.X)|isnan(c.Y);
c.X(M)=[];c.Y(M)=[];
Co=ProcessTananaFairbanks(odir,c,lateq,loneq);

