% load('/Users/chunlidai/Downloads/WV03_20150331124938_104001000986FC00_15MAR31124938-M1BS-500358011050_01_P008PSDndwiwl.mat')
% a=ndwil(~isnan(ndwil));b=ndwilw(~isnan(ndwilw));
% filename='1';
function [ithreshold,badflag,sigma0hat]=adaptivethresNDWI(filename,a,b)
% a land; b water; 
constant
ofile=[filename,'ab.mat'];
save(ofile,'a','b','-v7.3')

% a=random('norm',1,2,[10000,1]);b=random('norm',5,1,[10000,1]); %for test
mean1=nanmean(a);std1=nanstd(a);cnt1=sum(~isnan(a));
mean2=nanmean(b);std2=nanstd(b);cnt2=sum(~isnan(b));
p1=0.5;p2=0.5; %Note: p1+p2 should be 1.
%flagplot=1;
display(['Initial mean1 std1 mean2 std2:',num2str([mean1 std1 mean2 std2])])
%here, mean1 is for land;

%Estimate mean and std using Gaussian function to fit the PDF histogram.
inv=0.1/10;
edges=[-10:inv:10];
x=edges(1:end-1)+inv/2;

%get the edges
[edgesa]=histr(a);
[edgesb]=histr(b);

figure
hold all
hb=histogram(b,edgesb,'Normalization','pdf');
yb=hb.Values(:);hbBin=hb.BinWidth;% 
ha=histogram(a,edgesa,'Normalization','pdf','Facecolor','red'); %,'Edgecolor','green'
ya=ha.Values(:);haBin=ha.BinWidth;
sigma0hat=ones(2,1);
for k=[1,2]
    if k==1
        y=ya;% hb.Values;
        ksai0=[mean1 std1]';
        x=ha.BinEdges(1:end-1)+ha.BinWidth/2;
        xa=x;
    elseif k==2
        y=yb;% hb.Values;
        ksai0=[mean2 std2]';
        x=hb.BinEdges(1:end-1)+hb.BinWidth/2;
        xb=x;
    end
    lenf=length(y);m=2;
count=0; dksai=1;
while  dksai > 1e-9
mu0=ksai0(1);s0=ksai0(2);
E=exp(-(x(:)-mu0).^2/(2*s0^2));
y0=1./(sqrt(2*pi)*s0)*E;
AM=[1./(sqrt(2*pi)*s0)*E*2.*(x(:)-mu0)/(2*s0^2), (-1)/(sqrt(2*pi)*s0^2)*E+1/(sqrt(2*pi)*s0)*E.*(x(:)-mu0).^2/(4*s0^3)];
est=AM\(y-y0); %equal weight of y
ksai0=ksai0+est;
dksai=sum(abs(est));
etilde=y-y0;
count=count+1;
if (count>99) 
    badflag=1;ithreshold =threshold;sigma0hat=[1 1];
    warning([filename,': Gaussian fitting failed.'])
    return
end
end
sigma0hat(k)=sqrt(etilde'*etilde/(lenf-m));
if k==1
    mean1=ksai0(1);std1=ksai0(2);
elseif k==2
    mean2=ksai0(1);std2=ksai0(2);
end
end
%close all
display(['Fitted mean1 std1 mean2 std2:',num2str([mean1 std1 mean2 std2])])

badflag=0;
%if mean2<=mean1
%    badflag=1;ithreshold =threshold;
if std1==0 && std2==0
    ithreshold=(mean1+mean2)/2;
elseif std1==0 
    ithreshold=mean1+0.1; %arbituary value
elseif std2==0 %water has constant ndwi
    ithreshold=mean2-0.1; %arbituary value
else
    A=std1^2-std2^2;B=2*(mean1*std2^2-mean2*std1^2);
    C=std1^2*mean2^2-std2^2*mean1^2+2*std1^2*std2^2*log((std2*p1)/(std1*p2));
    D=B^2-4*A*C;
    if D < 0
            badflag=1;ithreshold =threshold;
    else
        roots=[(-B+sqrt(D))/(2*A),(-B-sqrt(D))/(2*A)];
        ithreshold=roots(roots>=min(mean1,mean2)&roots<=max(mean1,mean2));
        if length(ithreshold)~=1
            warning([filename,':ithreshold does not have one root! ithreshold=',num2str(ithreshold)])
            badflag=1;ithreshold =threshold;%sigma0hat=[1 1];return %keep the following plot
        end
    end
end

if flagplot==1
	ofile=[filename,'ab.mat'];
	save(ofile,'a','b','-v7.3')
	
        x=min([hb.BinEdges(:);ha.BinEdges(:)]):hb.BinWidth:max([hb.BinEdges(:);ha.BinEdges(:)]);
fitbimodal=p1/(sqrt(2*pi)*std1)*exp(-(x-mean1).^2/(2*std1^2))+p2/(sqrt(2*pi)*std2)*exp(-(x-mean2).^2/(2*std2^2));
fitbimodal1=p1/(sqrt(2*pi)*std1)*exp(-(x-mean1).^2/(2*std1^2)); %land
fitbimodal2=p2/(sqrt(2*pi)*std2)*exp(-(x-mean2).^2/(2*std2^2)); %water

figure
set(gcf,'Color','white')
set(gca,'FontSize', 12);
set(gcf, 'PaperPosition', [0 0 6 3]); 
set(gcf, 'PaperSize', [ 6 3]); 
hold all;
hold on
ha=histogram(a,edgesa,'Normalization','pdf','Facecolor','red','Edgecolor','red'); %,'Edgecolor','green'
hb=histogram(b,edgesb,'Normalization','pdf','Facecolor','blue','Edgecolor','blue');
% h=p1*ha.Values+p2*hb.Values;
% plot(x,h*1/p1,'go-')
plot(x,fitbimodal1*1/p1,'mo-')
plot(x,fitbimodal2*1/p2,'co-')
plot([ithreshold ithreshold],[0 10],'k-','linewidth',2)
x1=ithreshold;
hh=p2/(sqrt(2*pi)*std2)*exp(-(x1-mean2).^2/(2*std2^2))*1/p2; 
x1=mean2-abs(ithreshold-mean2);x2= mean2+abs(ithreshold-mean2);
%hold on;plot([x1,x2],[hh,hh],'r-')
% plot(x,fitbimodal*1/p1,'ko-')
legend('Land','Ocean','Land fitting','Ocean fitting','Location','Northwest')
box on
% xlabel('Brightness') 
% axis([-1.5 2 0 10])
ofile=['pic/',filename,'pdf'];
print('-dpdf','-r300',ofile) 
saveas(gcf,ofile,'fig')

end

x1=mean2-abs(ithreshold-mean2);x2= mean2+abs(ithreshold-mean2);

%ithreshold=[mean2-abs(ithreshold-mean2), mean2+abs(ithreshold-mean2)];

%estimate commission errors: the ratio of land pixels fall within the threshold range over the total land pixels.
M=xa>ithreshold;
comerr=sum(ya(M))*haBin;%since pdf: vi=ci/N*wi, where ci is count within bin i; wi is bin width; N is total count.
%ommision error:
M=xb<=ithreshold;
omierr=sum(yb(M))*hbBin;

%Gonzalez and Woods, 1992; eq 10.3-9;Minimize the error in making decisions that a given pixel belongs to an object or to the background.
%err=0.5*comerr+0.5*omierr;
err=comerr+omierr;
if (err >= comerrthres) badflag=1;end
sigma0hat=[comerr omierr];

close all
return
end


