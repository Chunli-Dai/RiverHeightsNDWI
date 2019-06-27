function [edges]=histr(a)
% input: a, a vector of brightness.
% Output:edge, for histogram with better intervals.
%Refers to:coastline/codec2/adaptivethresbp2.m
%	   rivergithub2/maskentropybp1.m
	constant

	% try change the BinWidth if there are zeros within the main lobe.
	mean1=nanmean(a);std1=nanstd(a);cnt1=sum(~isnan(a));
	hmin=exp(-9/2.)*0.5;%mu+-3*std range
	xamin=mean1-6*std1;xamax=mean1+6*std1;
	figure
	hold all
	ha=histogram(a,'Normalization','pdf','Facecolor','red'); %,'Edgecolor','green'
	inv=ha.BinWidth;
	edgesa=ha.BinEdges;xamin=min([edgesa(:);xamin]);xamax=max([edgesa(:);xamax]);
	edges=[xamin:inv:xamax];
	ha=histogram(a,edges,'Normalization','pdf','Facecolor','red'); %,'Edgecolor','green'
% 	ya=ha.Values(:);edgesa=ha.BinEdges;xa=ha.BinEdges(1:end-1)'+ha.BinWidth/2;
% 	%ha.BinWidth could be 'nonuniform'
	ya=ha.Values(:);edgesa=ha.BinEdges;xa=ha.BinEdges(1:end-1)'+inv/2;
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
	inv=inv*scale;edges=[xamin:inv:xamax];
	ha=histogram(a,edges,'Normalization','pdf','Facecolor','yellow');
	ya=ha.Values(:);edgesa=ha.BinEdges;xa=ha.BinEdges(1:end-1)'+inv/2;

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
	inv=inv*scale;edges=[xamin:inv:xamax];
	ha=histogram(a,edges,'Normalization','pdf','Facecolor','green');
	ya=ha.Values(:);edgesa=ha.BinEdges;xa=ha.BinEdges(1:end-1)'+inv/2;
	
	close all

return
end
