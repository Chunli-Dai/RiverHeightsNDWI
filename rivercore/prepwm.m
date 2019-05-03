function [wmi]=prepwm(xout,yout,jump,prob,stepi)
%Prepare water mask for a priori information of water classification
%	Output:
%       wmi.z;%river mask (50% (step 2) or 90% (step 3))
%       wmi.w50; %water mask 50% probability, for getting land area
%       wmi.buf; % river buffer, twice the width of river mask with 50% probability.
	constant

%   wm1.x=[];wm1.y=[];wm1.z=[];
    wm.z= jump;%water mask  (50% (step 2) or 90% (step 3))
    wm.x=xout;wm.y=yout;   
	wm=cropmatrix(wm,wm.z);
	wm1=wm; %For the river mask!;
	
	%river mask only (no lakes), for a priori water mask.
    % % Assumption: the river mask is very reliable with 50% or 90% water probability.
	if 0
            BW=wm.z;
            Medge=(BW==-1);
            BW(BW==-1)=0;
            Modj= bwareaopen(BW, round(lakearea/2/resr/resr)); %remove clusters of lakes 500m by 500m
            Modj=int8(Modj);
            Modj(Medge)=-1;
            wm1=wm;wm1.z=Modj;
            wm1=mask2river(wm1); %remove lakes, no fill
        %   tz = interp2(wm1.x,wm1.y,wm1.z,double(wm.x),wm.y','*nearest',-1);
        %   wm1=wm;wm1.z=tz; %keep wm1(river mask without lakes) unchanged;still could be lakes undetected, or small tributaries.
	end % if 0

	if stepi ==3 
	    %get the 50% mask for buffering.
	    iprobthre=probthre;
            [nsuby,nsubx]=size(prob);
	    jump=-1*ones(nsuby,nsubx,'int8'); %%-1,non value;1 water; 0 non-water
	    jump(prob>=iprobthre&prob~=255)=1;
	    jump(prob<iprobthre)=0;
            tz = interp2(xout,yout,jump,wm.x,wm.y','*nearest',-1);
	    wm.z=tz; %keep wm (water mask with 50%) unchanged;
	    clear prob jump
	end
	
	%Get river centerline from the 50% water mask, wm.
        %Get river buffer zone;
	tic
        BW=wm.z;
        Medge=(BW==-1);
        BW(BW==-1)=0;
        Modj= bwareaopen(BW, round(lakearea/2/resr/resr)); %remove clusters of lakes 500m by 500m
        Modj=int8(Modj);
        Modj(Medge)=-1;
        data=wm;data.z=Modj;
        data=mask2river(data); %remove lakes, no fill
    %   tz = interp2(data.x,data.y,data.z,wm.x,wm.y','*nearest',-1);
    %   data=wm;data.z=tz; %50% water probability without lakes.
	fprintf('\nGetting a priori river mask (stacked), excluding lakes.')
	toc

	%river mask only (no lakes), for a priori water mask.
	wm1.z=(wm1.z==1)&(data.z==1); %keep wm1(river mask without lakes) unchanged;still could be lakes undetected, or small tributaries.
    
% 	data=wm;
	resx=mean(data.x(2:end)-data.x(1:end-1));resy=mean(data.y(2:end)-data.y(1:end-1));
	resr=mean([abs(resx),abs(resy)]);
    
    if exist('clsv2.mat','file')
         load clsv2.mat %use the saved one %hi
	     %Notice: centerline cannot have NaNs; center line monotonic
	     M=isnan(c.X)|isnan(c.Y);
	     c.X(M)=[];c.Y(M)=[];
    else
        res1=resr;

        tz=data.z;
        tz(tz==-1)=0;
        Modj= bwareaopen(tz, round(lakearea/res1/res1)); %remove small clusters
        Modfil = bwareaopen(~Modj, round(cloudarea/res1/res1)); %fill small areas: 1e4*4m^2
        Modfil=~Modfil;
        data.z=Modfil;

        npt=sum(sum(data.z==1));
        if npt>0

            try     
                tic
                [c]=mask2centerline(data);
                save cl1.mat c
                fprintf('\nGetting river centerline from river mask.')
                toc
            catch e
                fprintf('There was an error! The message was:\n%s',e.message);
                fprintf('\n')
            end
            %adjust the centerline to go uphill for ProcessTananaFairbanks.m.
            width=500; %m 
            pt=[c.X(1) c.Y(1);c.X(end) c.Y(end)];
            %use data0r to get the height
            if 0
            c.X=flip(c.X);
            c.Y=flip(c.Y);
            end
        else %npt<=0
           fprintf(['No water pixels in the water mask. ']);
        end % if npt
    end % if exist clsv2.mat

	%get river width
	tic
 	infile='watermaskstacked';
	[gagewidth,widthp]=getwidth(data,infile,c);
	%get buffer width: 80 percentile of river width times 2.
	width80 = prctile(widthp.y(widthp.y>0),80);	
	fprintf(['\n80 percentile river width is:',num2str(width80),'m.']);
	if width80<10||isnan(width80)||isinf(width80)
	   width80=200;
	   fprintf(['\n Unable to get a reasonable width, use the fixed width:',num2str(width80),'m.']);
	end
	toc

	%get bufferzone along river centerline, to remove some tributaries.
	buf=zeros(size(data.z));
	[clx,cly]=polarstereo_fwd(c.Y,c.X,[], [],70,-45);
	%polar stereographic coordinates to image coordinates.
	[ny,nx]=size(data.z);
	clear cl
	cl(:,1)=round((clx-data.x(1))/resr)+1;
	cl(:,2)=round(-(cly-data.y(1))/resr)+1;
	M=cl(:,1)>=1&cl(:,1)<=nx&cl(:,2)>=1&cl(:,2)<=ny;
	cl(~M,:)=[];
%	buf(cl(:,2),cl(:,1))=1;
	for j=1:length(cl(:,1))
	buf(cl(j,2),cl(j,1))=1;
	end
	widpix=round(width80/resr);
    ncl=3; % expand along centerline by ncl times; try 10, 5 3
	widpix2=round(width80/resr*ncl);
	clbuf= imdilate(buf, ones(widpix2*2)); % width expansion

    
    %Get buffer zone along river water body
	Modj=(data.z==1)&clbuf; %remove undetected lakes or unwanted tributaries in the far field.
    %futher remove narrow rivers <80m
    if 1 %very tricky; can remove main stream
    widsm=40;%m;any river tributar less than widsm m will be removed for the generation of buff zone.
    widsmpix=round(widsm/resr);
    Modj=imerode(Modj, ones(widsmpix));
    Modj=bwareaopen(Modj, round(lakearea/10/resr/resr));
    Modj=imdilate(Modj, ones(widsmpix));
    end
	Modj= bwareaopen(Modj, round(lakearea/resr/resr)); %remove small clusters
	rivbuf=imdilate(Modj, ones(widpix*2)); % width expansion; 

	%get river mask without lake (50% (step 2) or 90% (step 3))
	Modj=(wm1.z==1)&rivbuf;%remove undetected lakes or unwanted tributaries.
	Modj= bwareaopen(Modj, round(lakearea/resr/resr)); %remove small clusters
	rivm=Modj;

	wmi=wm1;
	wmi.z=rivm;%wmi.z=wm1.z;%river mask (50% (step 2) or 90% (step 3)). 1 water, 0 land
	wmi.w50=wm.z; %water mask 50% probability, for getting land area. 1 water, 0 land, -1 edge
	wmi.buf=rivbuf; % river buffer, expand the river mask of 50% probability by twice the width. 1,within buffer, 0 otherwise.

projstr='polar stereo north';
ofile=['riverbuff.tif'];
writeGeotiff(ofile,wmi.x,wmi.y,uint8(wmi.buf),1,255,projstr)

figure
set(gcf,'Color','white')
set(gca,'FontSize', 12);
set(gcf, 'PaperPosition', [0 0 6 3]);
set(gcf, 'PaperSize', [ 6 3]);
hold all;
imagesc(wmi.x*1e-3,wmi.y*1e-3,double(wmi.z))
colormap jet;colorbar
title(['Stacked river mask']);
axis equal; box on
caxis([-1 1])
hold on
ofile=['rivermaskstackedStep',num2str(stepi)];
print('-dpng','-r400',ofile)
%saveas(gcf,ofile,'fig')
close all

figure
set(gcf,'Color','white')
set(gca,'FontSize', 12);
set(gcf, 'PaperPosition', [0 0 6 3]);
set(gcf, 'PaperSize', [ 6 3]);
hold all;
imagesc(wmi.x*1e-3,wmi.y*1e-3,double(wmi.w50))
colormap jet;colorbar
title(['Water mask based on 50% probability']);
axis equal; box on
caxis([-1 1])
hold on
ofile=['watermaskstacked50pStep',num2str(stepi)];
print('-dpng','-r400',ofile)
close all

figure
set(gcf,'Color','white')
set(gca,'FontSize', 12);
set(gcf, 'PaperPosition', [0 0 6 3]);
set(gcf, 'PaperSize', [ 6 3]);
hold all;
imagesc(wmi.x*1e-3,wmi.y*1e-3,double(wmi.w50),'alphadata',double(wmi.buf))
colormap jet;colorbar
title(['Water mask based on 50% probability within river buffer']);
axis equal; box on
caxis([-1 1])
hold on
ofile=['riverbuffStep',num2str(stepi)];
print('-dpng','-r400',ofile)

close all
