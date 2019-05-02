%this document aims to find the centerline using Meansift method.
%  oringinal the meanshift is designed for road detection:
% Xiangyun Hu; Yijing Li; Jie Shan; Jianqing Zhang; Yongjun Zhang, "Road Centerline Extraction in Complex Urban Scenes From LiDAR Data Based on Multiple Features," Geoscience and Remote Sensing, IEEE Transactions on , vol.52, no.11, pp.7448,7456, Nov. 2014
% doi: 10.1109/TGRS.2014.2312793
% http://ieeexplore.ieee.org/xpl/articleDetails.jsp?arnumber=6803040
%
%here it is implemented for river centerline extraction
%
%%---------------the first step is to find the river centre line -------------
clear all;
clc;

%%load a sample data
%river_mask=imread('E:\Radarsat_CSA\River_connection\Bow_river_mask_after_SAR.tif');
[data,info]=enviread('.\Test_Data\Bow_river_refined_mask_after_SAR_complete_envi.dat','.\Test_Data\Bow_river_refined_mask_after_SAR_complete_envi.hdr');

% experiment on a subset:  process a subset of this data for experimental purpose
river_mask_subset=data(250:1515, 900:2480);
river_mask_subset=river_mask_subset(676:1252, 326:644);
%river_mask_subset=river_mask_subset(284:473, 2:312);
%river_mask_subset=river_mask;

%the output root
outputpath='BowRiver_mask_subet_skeleton_';

%the mean shift method to find the central line of rivers
bankwidth = 20;

%tile processing
% [size_col,size_row]=size(river_mask);
% 
% Tile_x_width=250;  %size_col/(tile_col*0.9-0.1); %given 10% of overlap on both side  
% Tile_y_width=250;  %size_row/(tile_row*0.9-0.1); %given 10% of overlap on both side
% 
% tile_col=10;  %the number of tiles in column and row directions
% tile_row=5;
% 
% for tileX=1:tile_col
%     for tileY=1:tile_row

        %obtain a subset for processing
        %river_mask_subset=river_mask(((tileX-1)*Tile_x_width+1):(tileX*Tile_x_width), ((tileY-1)*Tile_y_width+1):(tileY*Tile_y_width));   
        
        %x=[rand(1,10000)*10; rand(1,10000)*0.5];
        [idy idx]=find(river_mask_subset>-50);
        River_mask_coords=([idx';idy']);
        
%         idx=idx.+Tile_y_width;
%         idy=idy.+Tile_x_width;
        
        %tic
        %[clustCent,point2cluster,clustMembsCell] = MeanShiftCluster(River_mask_coords,bankwidth);
        %[clustCent,point2cluster,clustMembsCell] = MeanShiftCentralLine_backup(River_mask_coords,bankwidth);
        [Meanshift_Pts Window_Size_Pts]= MeanShiftCentralLine(River_mask_coords,bankwidth);
%         for i=1:1
%             bankwidth = bankwidth/2;
%             Meanshift_Pts = MeanShiftCentralLine(Meanshift_Pts,bankwidth);
%         end
        %toc
%     end
%     display(tileY)
% end

% con = struct('xc',[0,2*pi]);
% pp1 = splinefit(Meanshift_Pts(1,:),Meanshift_Pts(2,:),8,con,0.25); % Robust fitting
% xx = linspace(min(Meanshift_Pts(1,:)),max(Meanshift_Pts(1,:)),1000);
% spline_fit = ppval(pp1,xx);

%%save the result
%%read the hdr file to find the file offset at X and Y , as well as the
%%point space
% Meanshift_Pts(1,:)=Meanshift_Pts(1,:)*2+704344;
% Meanshift_Pts(2,:)=Meanshift_Pts(2,:)*-2+5660810;
mkdir('.\Result\')
fileID = fopen('.\Result\Bow_river_pts_after_Meanshift.csv','w');
fprintf(fileID,'X,Y\n');
fprintf(fileID,'%10.4f, %10.4f\n',[Meanshift_Pts(1,:)*2+704344; Meanshift_Pts(2,:)*-2+5660810]);  %export the point but keep the UTM coords
fclose(fileID);
%%this point result file ,as a CSV, can be import and overlapped directly
%%with other layers in ARCGIS, by importing --> adding XY. menu.

% plot the result
if (1)
	figure;
	imshow(river_mask_subset>-50)
	% scatter(River_mask_coords(1,:),River_mask_coords(2,:),'.k');

	cVec = 'bgrcmykbgrcmykbgrcmykbgrcmyk';
	hold on;
	gscatter(Meanshift_Pts(1,:),Meanshift_Pts(2,:),Window_Size_Pts,cVec(unique(fix(Window_Size_Pts/bankwidth))), '.');
	%scatter(Meanshift_Pts(1,:),Meanshift_Pts(2,:),'.r','MarkerFaceColor',cVec(fix(Window_Size_Pts/bankwidth)));

	%plot(xx,spline_fit,'r');

	hold off;
end
% numClust = length(clustMembsCell);
% 
% %to show the result
% figure;
% scatter(River_mask_coords(1,:),River_mask_coords(2,:));
% hold on;
% % cVec = 'bgrcmykbgrcmykbgrcmykbgrcmyk';%, cVec = [cVec cVec];
% for k = 1: numClust %min(numClust,length(cVec))
%     myClustCen = clustCent(:,k);
%     plot(myClustCen(1),myClustCen(2),'s','MarkerEdgeColor','k');  %,'MarkerFaceColor',cVec(k)
% end 
% % [n,idx]=sort(clustCent(1,:)); %sort by the X
% % clustCent=clustCent(:,idx);
% % plot(clustCent(1,:),clustCent(2,:),'-r')
% hold off;

return



%the original image
%imshow(river_mask_subset>-50)
imwrite(uint8(river_mask_subset>-50)*255,strcat(outputpath,'original.tif'),'TIFF');

%show the default methods from matlab to find the outline and centreline
river_mask_subset_outline=bwmorph(river_mask_subset>-50,'remove');
%imshow(river_mask_subset_outline);
imwrite(uint8(river_mask_subset_outline)*255,strcat(outputpath,'original_outline.tif'),'TIFF');
%default skeleton method from Matlab
river_mask_subset_centreline=bwmorph(river_mask_subset>-50,'skel',Inf);
%figure; imshow(river_mask_subset_centreline)
imwrite(uint8(river_mask_subset_centreline)*255,strcat(outputpath,'original_Ske.tif'),'TIFF');


%using the DIP processing package

river_mask_subet_skeleton=dip_EuclideanSkeleton(river_mask_subset>-50,'natural',0); 
%figure; imshow(uint8(dip_array(river_mask_subet_skeleton))>0)
% print -deps2c river_mask_subet_skeleton_natural.eps  % -r0  means screen resolution; //print as EPS file
% print -dpdf  river_mask_subet_skeleton_natural.pdf %print as PDF file
% print -dtiff river_mask_subet_skeleton_natural.tif %print as TIFF file
imwrite(uint8(dip_array(river_mask_subet_skeleton))*255,strcat(outputpath,'natural.tif'),'TIFF');

river_mask_subet_skeleton=dip_EuclideanSkeleton(river_mask_subset>-50,'1neighbor',0); 
%figure; imshow(uint8(dip_array(river_mask_subet_skeleton))>0)
imwrite(uint8(dip_array(river_mask_subet_skeleton))*255,strcat(outputpath,'1neighbor.tif'),'TIFF');

river_mask_subet_skeleton=dip_EuclideanSkeleton(river_mask_subset>-50,'2neighbors',0); 
%figure; imshow(uint8(dip_array(river_mask_subet_skeleton))>0)
imwrite(uint8(dip_array(river_mask_subet_skeleton))*255,strcat(outputpath,'2neighbors.tif'),'TIFF');

river_mask_subet_skeleton=dip_EuclideanSkeleton(river_mask_subset>-50,'3neighbors',0); 
%figure; imshow(uint8(dip_array(river_mask_subet_skeleton))>0)
imwrite(uint8(dip_array(river_mask_subet_skeleton))*255,strcat(outputpath,'3neighbors.tif'),'TIFF');