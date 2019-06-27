function [ithres]=histhres(name,a)
% Find the threshold of brightness for separating water using histogram method.
% input: a, a vector of brightness.
%	name, filename for save figures
% Output:ithres(1:2), the lower and upper threshold of water pixels
%Refers to:coastline/codec2/adaptivethresbp2.m
%	   rivergithub2/maskentropybp1.m
	constant
%save ndwil20090525.mat a
	ofile=[name,'ndwil.mat'];
	save(ofile,'a','-v7.3')

	% try change the BinWidth if there are zeros within the main lobe.
	mean1=nanmean(a);std1=nanstd(a);cnt1=sum(~isnan(a));
	hmin=exp(-9/2.)*0.5;%mu+-3*std range
	xamin=mean1-6*std1;xamax=mean1+6*std1;
	figure
	hold all
	ha=histogram(a,'Normalization','pdf','Facecolor','red'); %,'Edgecolor','green'
	inv=ha.BinWidth;edges=[xamin:inv:xamax];
	ha=histogram(a,edges,'Normalization','pdf','Facecolor','red'); %,'Edgecolor','green'
	ya=ha.Values(:);edgesa=ha.BinEdges;xa=ha.BinEdges(1:end-1)'+ha.BinWidth/2;
	% id=xa>(mean1-3*std1)&xa<(mean1+3*std1)&ya>=max(ya)*hmin;xa(id)=[];ya(id)=[]; 

	% adjust BinWdith
%id=xa>(xamin)&xa<(xamax)&ya<=max(ya)*hmin;
	id=xa>(xamin)&xa<(xamax)&ya<=max(ya)*hmin*0;
	%ratio=sum(id)/(sum(xa>(xamin)&xa<(xamax))-sum(id));
	nall=sum(xa>(xamin)&xa<(xamax));nbad=sum(id);
% 	ng=nall-nbad;
% 	ratio=nall/ng;
    % Reduce Binwidth
	if 1
	nbins=4000; %larger the better, the void ones will be adjusted again.
	    if nall<nbins
	       scale=nall/nbins;
	    else
	       scale=1;
	    end
    else
        scale=1/10; %somehow more fluctuation
    end
	inv=ha.BinWidth*scale;edges=[xamin:inv:xamax];
	ha=histogram(a,edges,'Normalization','pdf','Facecolor','yellow');
	ya=ha.Values(:);edgesa=ha.BinEdges;xa=ha.BinEdges(1:end-1)'+ha.BinWidth/2;

	% adjust BinWdith Again
	id=xa>(xamin)&xa<(xamax)&ya<=max(ya)*hmin*0;
	nall=sum(xa>(xamin)&xa<(xamax));nbad=sum(id);
	ng=nall-nbad;
	ratio=nall/ng;%>=1
	if ratio ==1 
	    %good
	    scale=1;
	else %increase the binwidth
% 	    scale=ratio*3;%reduce the spiky signals in dense distribution.; but the interval can be too big 
	    scale=ratio;
	end
	inv=ha.BinWidth*scale;edges=[xamin:inv:xamax];
	ha=histogram(a,edges,'Normalization','pdf','Facecolor','green');
	ya=ha.Values(:);edgesa=ha.BinEdges;xa=ha.BinEdges(1:end-1)'+ha.BinWidth/2;
	
	% check the historgram 
	ndwil=a;
        [hpdf,edges]=histcounts(ndwil,edges,'Normalization','pdf');
        %[hpdf,edges]=histcounts(ndwil,'Normalization','pdf');
        %smooth to reduce flunctuation
        if 1 %better; no shift
        hpdfs=zeros(size(hpdf));
        hpdfs(end)=hpdf(end)/2;
        hpdfs(1:end-1)=(hpdf(2:end)+hpdf(1:end-1))/2;
        else % shift one bin to right
            hpdfs=movmean(hpdf,2); %can only be even number to reduce the fluctuation
        end
        hpdf=hpdfs;
        [hmax,idm]=max(hpdf);
        meanor=(edges(idm)+edges(idm+1))/2;
        hh=hmax*exp(-1/2.); %1 sigma,About 68% of values drawn from a normal distribution are within one standard deviation σ away from the mean;
        hh=hmax*exp(-4/2.); %2 sigma,95% 
        id2=idm+1:length(hpdf);
        id2i=find((hpdf(id2)-hpdf(id2-1))>0&hpdf(id2)<hmax/2.);%find the curve turning point
        if ~isempty(id2i); id2=id2(1):id2(id2i(1));end
        [~,id2i]=min(abs(hpdf(id2)-hh));
        hw2=id2(id2i(1));

        id1=1:min(idm-1,length(hpdf));
        id1i=find((hpdf(id1+1)-hpdf(id1))<0&hpdf(id1)<hmax/2.);%find the curve turning point
        if ~isempty(id1i); id1=id1(id1i(end)):id1(end);end
        [~,id1i]=min(abs(hpdf(id1)-hh));
        hw1=id1(id1i(end));
   %    hw=edges(hw2)-edges(hw1);
    	if flagplot==1;figure;plot(edges(2:end),hpdf,'b+-')
        	hold on;plot([edges(hw1+1),edges(hw2+1)],[hh,hh],'r-')
		ofile=[name,'pdf'];
		saveas(gcf,ofile,'fig')
    	end
	
	ithres(1)=edges(hw1+1);ithres(2)=edges(hw2+1);
	close all

return
end
