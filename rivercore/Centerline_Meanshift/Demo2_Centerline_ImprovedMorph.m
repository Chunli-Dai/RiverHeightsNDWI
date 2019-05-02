%
%here it is implemented for river centerline extraction using an improved centerline extraction method
%

%read the water mask after river connection
[data,info]=enviread('.\Test_Data\Bow_river_refined_mask_after_SAR_complete_envi.dat','.\Test_Data\Bow_river_refined_mask_after_SAR_complete_envi.hdr');

% run the centerline detection method
nedgelist=Improved_CentralLine(Mask_Path_filter,Max_River_width)

%export the result
str='.\Results\';
mkdir(str);
ExportEdgelist2ShpFile(nedgelist,704344,5660810,2,strcat(str,'river_complete_nedgelist'));  %save the result using the UTM coords
