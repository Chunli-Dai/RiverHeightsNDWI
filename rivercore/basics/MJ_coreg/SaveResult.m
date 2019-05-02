function [out] = SaveResult(varargin)
ref_select_nor      = varargin{1};
tar_select_nor_z    = varargin{2};
ref_select_H        = varargin{3};
tar_H_array         = varargin{4};
RefImg_gridspace    = varargin{5};
TarImg_gridspace    = varargin{6};
Boundary            = varargin{7};
Tarboundary         = varargin{8};
Mean                = varargin{9};
Scale               = varargin{10};
ref_select_size     = varargin{11};
X                   = varargin{12};
Method              = varargin{13};
Tar_row             = varargin{14};
Tar_col             = varargin{15};
thresh_h            = varargin{16};
H_accuracy          = varargin{17};


null_flag           = (tar_select_nor_z > -1 & tar_select_nor_z < 1);

%calculate all height and distance on location in ref img
coord{1}            = ref_select_nor{1};
coord{2}            = ref_select_nor{2};
coord{3}            = tar_select_nor_z;
clear Diff_array;

transformed_coord   = ConformalTransform(coord,X);
flag                = transformed_coord{1} > -1 & transformed_coord{1} < 1 & transformed_coord{2} > -1 & transformed_coord{2} < 1;

Diff_array_ori      = ref_select_nor{3} - tar_select_nor_z;
Diff_array          = ref_select_nor{3}(null_flag & flag) - tar_select_nor_z(null_flag & flag);
Mean_s              = mean(Diff_array);
Std                 = std(Diff_array);

index_1             = (null_flag & flag);
t_size              = size(transformed_coord{1});
iter_tar_height     = ones(t_size(1),1)*(-9999);

diff_distance_ortho = ones(t_size(1),1)*(-9999);
clear temp_coord t_image flag;

if strcmp(Method,'3D')
    diff_distance_surface = ones(t_size(1),1)*(-9999);
    [iter_tar_height_t G{1}  G{2}  G{3}  D  garbage garbage tar_flag] = interpolationTIN_Modify(tar_H_array,TarImg_gridspace,transformed_coord{1}, transformed_coord{2},Mean, Scale,Tarboundary,X,[Tar_row Tar_col],0);
    clear D tar_hroh Tar_adj_h;

    diff_distance_surface(index_1)  = (ref_select_nor{3}(index_1) - iter_tar_height_t(index_1))*Scale(3);
    diff_distance_ortho             = diff_distance_surface;
    iter_tar_height(index_1 & tar_flag)       = iter_tar_height_t(index_1 & tar_flag)*Scale(3) + Mean(3);
    
    t_coord{1}                      = ref_select_nor{1};
    t_coord{2}                      = ref_select_nor{2}; 
    t_coord{3}                      = iter_tar_height_t;
    [ref_X ref_Y ref_Z select_flag] = SufaceDistance(ref_select_H,G,t_coord,RefImg_gridspace,Boundary,[ref_select_size(1) ref_select_size(2)], Mean, Scale, ref_select_nor{3});
    
    dx                              = (t_coord{1} - ref_X)*Scale(1);
    dy                              = (t_coord{2} - ref_Y)*Scale(2);
    dz                              = (t_coord{3} - ref_Z)*Scale(3);
    distance                        = sqrt(dx.*dx + dy.*dy + dz.*dz);

    tt_index                        = find(dz < 0);
    distance(tt_index)              = -distance(tt_index);
            
    coord_flag                      =  ref_X > -1 & ref_X < 1 & ref_Y > -1 & ref_Y < 1 & ref_Z > -1 & ref_Z < 1;        
    index_1                         = index_1 & tar_flag & select_flag & coord_flag;

    diff_distance_surface(index_1)  = -distance(index_1);

    null_index                      = find(diff_distance_surface < -9999 | diff_distance_surface > 5000);
    t_size                          = size(null_index);
    diff_distance_surface(null_index)=ones(t_size(1),1)*(-9999);

    clear RefImg ref_select_H ref_Z Gr ref_flag G D tar_flag t_coord ref_X ref_Y ref_Z select_flag flag tar_hroh Tar_adj_h dx dy dz distance;  
else
    R                               = RotationMatrix([X(2,1) X(3,1) X(4,1)]);
    t_coord                         = Denormalize_coord(transformed_coord,Mean,Scale,'XY');
    t_coord                         = GeoToImage(t_coord,Tarboundary,TarImg_gridspace);
    
    [tar_pt_height_1 flag]          = interpolation(tar_H_array,t_coord{1},t_coord{2},Tar_row, Tar_col);
    tar_pt_height                   = (tar_pt_height_1-Mean(3))/Scale(3);
    iter_tar_height                 = (X(1,1)*(R(1,3)*transformed_coord{1} + R(2,3)*transformed_coord{2} + R(3,3)*tar_pt_height) + X(7,1));
    index_1                         = index_1 & flag;
    diff_distance_ortho(index_1)    = (ref_select_nor{3}(index_1) - iter_tar_height(index_1))*Scale(3);
    iter_tar_height                 = iter_tar_height*Scale(3) + Mean(3);
    clear tar_pt_height_1 tar_pt_height;
end
clear coord tar_H_array index_1;

clear transformed_coord tar_select_nor_h_ori flag DiffArray;

%save co-registration img
index                               = [1:ref_select_size(1)*ref_select_size(2)]'; 

value                               = iter_tar_height;  
out_array                           = savefromindex(ref_select_size,index,value,-9999);
out                                 = out_array;
%geotiffwrite(save_img_file_name, bbox, out_array,option.bit_depth,option);
%geotiffwrite(save_name         , bbox, out_array,option.bit_depth,option);

clear index value out_array diff_h iter_tar_height ref_select_nor_ori diff_distance diff_distance_ortho;