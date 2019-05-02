function [dX4Sg]=pxpymono(dX4Sg2,f2is,fdir2is,XYb2is,fis,fdiris,XYbis)
%To collect the translational parameters of reference mono images that are components of strip files.
%And to coregsiter all mono image to the reference image
%dX4Sg2,f2is,fdir2is: Translational parameters, filenames, and directories of strip files.
% fis: filenames of mono images.
%output: Mrefmo : flag for fis that whether it is reference mono image or not.
%	psrefmo: the pzxy for the reference mono image.

n=length(fis);n2=length(f2is);
Mrefmo=false(n,1);
psrefmo=zeros(n,3);

for j=1:n2
stripmetafile= [fdir2is{j},'/',f2is{j}];
[~,filename,~]=fileparts(stripmetafile);
satname=filename(1:4);
p3=dX4Sg2(j,:);
end % for j


%coregister all mono image to the reference mono images.
idtar=1:n;%find(Mrefmo~=1);
idtar=79:80; %test two mono images of strip WV02_20160827_103001005A735B00_103001005CAE7900_seg1_2m_dem.tif
fprintf(['Number of targe images:',num2str(length(idtar))])
for j=1:length(idtar)
	is=idtar(j);
	tarimage=[fdiris{is},'/',fis{is}];

        Xb=XYbis{is}(:,1); Yb=XYbis{is}(:,2);
        rang0=[min(Xb) max(Xb) min(Yb) max(Yb)];
	dx=2;dy=dx;
	nwx=round((rang0(2)-rang0(1))/dx);nwy=round((rang0(4)-rang0(3))/dy);
	sx=Xb;sy=Yb;
        idx=round((sx-rang0(1))/dx)+1;idy=round((sy-rang0(3))/(dy))+1;
        sm0=poly2mask(idx,idy,nwy,nwx); % fast, apply to each polygon one by one.

	%find the most overlapped reference image.
	nptall=zeros(n2,1);
        for jj=1:n2
	  k=jj;
          sx=XYb2is{k}(:,1); sy=XYb2is{k}(:,2);
          idx=round((sx-rang0(1))/dx)+1;idy=round((sy-rang0(3))/(dy))+1;
          sm=poly2mask(idx,idy,nwy,nwx); % fast, apply to each polygon one by one.
          nptall(jj)=sum(sum(sm.*sm0))/sum(sum(sm0))*100;
        end %k
	[~,jjref]=max(nptall);
	isref=jjref;
	refimage=[fdir2is{isref},'/',f2is{isref}];
    	fprintf(['\n j tar ref: ',num2str([j]),' ', tarimage,' ',refimage])

	%prepare images for coregistration
	[refimagep,tarimagep,dataref,datatar]=prepareMJ(refimage,tarimage);

	odir='./outcoreg/'
	if ~exist(odir,'dir')
	  mkdir(odir)
	end
	str=['time setsm -Coreg 1 -image ',refimagep,' -image ', tarimagep, ' -outpath ', odir];
	fprintf([str])
	[status, cmdout]=system(str);
	
	%plot the animation of images and control points before and after coregistration
	%refers to /home/dai.56/arcticdemapp/river/riverwork/coregtest1/plotcontrolpts.m
	cpts=load([odir,'/txt/GCPs_Image_ID_1_level_0.txt']);
	%Before Coregistration
    rangtar=[min(datatar.x) max(datatar.x) min(datatar.y) max(datatar.y)];
    rangref=[min(dataref.x) max(dataref.x) min(dataref.y) max(dataref.y)];
    rangeov=[max(rangtar(1),rangref(1)),min(rangtar(2),rangref(2)), max(rangtar(3),rangref(3)),min(rangtar(4),rangref(4))];
    range0=rangeov;
	range0=[-2231000 -2230800 549200 549400]; %%sag river %hi

	[co]=testgif(dataref,datatar,range0,cpts,refimagep,tarimagep,1);
	
	%get coregistration parameters
	%txy=[-1.86 5.04 ];
	coregfile=[odir,'coreg_result.txt'];
	c=textread(coregfile,'%s','delimiter','\n');	
	r=find(~cellfun(@isempty,strfind(c,tarimagep)));
	%Ty[meter]       Tx[meter]       avg_roh(average correlation)
	c2=c{r};
 	r1=strfind(c2,'tif');c2([1:r1(1)+2])='';
	[tmp]=sscanf(c2, '%f',[1,5]);

	txy=[tmp(4), tmp(3)]; rho=tmp(5);
	
	datatarc=datatar;datatarc.x=datatar.x-txy(1);datatarc.y=datatar.y-txy(2);
	cpts(:,3)=cpts(:,3)-txy(1); cpts(:,4)=cpts(:,4)-txy(2);
	[co]=testgif(dataref,datatarc,range0,cpts,refimagep,tarimagep,2);

	system(['rm ',refimagep, ' ',tarimagep])
	system(['rm ',odir, '/tmp/*'])
	
	ps=[0 txy]; %zxy
%	psref=psrefmo(isref,1:3);
	psref=[dX4Sg2(isref,1:3)];psref(1)=0;
	dX4Sg(is,1:3)=ps+psref;
    	fprintf(['\n j pzxy pzxytile: ',num2str([j, ps,ps+psref]),' ',tarimage,' ',refimage])

end

return
end
