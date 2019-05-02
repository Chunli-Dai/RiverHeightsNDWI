%this document aims to find the centerline using an improved morphology method and the Meansift method.
%
% you need to include the ENVI header into the MATLAB PATH
%
%  oringinal the meanshift is designed for road detection:
% Xiangyun Hu; Yijing Li; Jie Shan; Jianqing Zhang; Yongjun Zhang, "Road Centerline Extraction in Complex Urban Scenes From LiDAR Data Based on Multiple Features," Geoscience and Remote Sensing, IEEE Transactions on , vol.52, no.11, pp.7448,7456, Nov. 2014
% doi: 10.1109/TGRS.2014.2312793
% http://ieeexplore.ieee.org/xpl/articleDetails.jsp?arnumber=6803040
%
clear all;
clc;

%%load a sample data
[data,info]=enviread('./Test_Data/River_final_mask.hdf','./Test_Data/River_final_mask.hdr');


%
%  --------step 1-------------here it is implemented for river centerline
%  extraction using an improved centerline extraction method
%  -------------------------------------
%

% run the centerline detection method
Max_River_width=170;  %an estimated river width 
nedgelist=Improved_CentralLine(data,Max_River_width); %run the core centerline method

%export the result
str='.\Results\';
mkdir(str);
ExportEdgelist2ShpFile(nedgelist,704344,5660810,2,strcat(str,'river_complete_nedgelist.shp'));  %save the result using the UTM coords

%
%  --------step 2-------------here it is implemented for river centerline
%  extraction using an Mean Shift method
%  -------------------------------------
%


% experiment on a subset:  process a subset of this data for experimental purpose
river_mask_subset=data(250:1515, 900:2480);
river_mask_subset=river_mask_subset(676:1252, 326:644);

%the output root
outputpath='BowRiver_mask_subet_skeleton_';

%the mean shift method to find the central line of rivers
bandwidth = 20;

 [idy idx]=find(river_mask_subset>-50);
 River_mask_coords=([idx';idy']);
        
 [Meanshift_Pts Window_Size_Pts]= MeanShiftCentralLine(River_mask_coords,bandwidth);


%%save the result

% mkdir('.\Results\')
fileID = fopen('.\Results\Bow_river_pts_after_Meanshift.csv','w');
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
	gscatter(Meanshift_Pts(1,:),Meanshift_Pts(2,:),Window_Size_Pts,cVec(unique(fix(Window_Size_Pts/bandwidth))), '.');
	hold off;
end


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