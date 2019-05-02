function [diff_distance index_1]=SaveBeforeData(varargin)
tar_flag            = varargin{1};
ref_select_nor      = varargin{2};
tar_select_nor_h    = varargin{3};
ref_select_H        = varargin{4};
G                   = varargin{5};
RefImg_gridspace    = varargin{6};
TarImg_gridspace    = varargin{7};
Boundary            = varargin{8};
Tarboundary         = varargin{9};
Mean                = varargin{10};
Scale               = varargin{11};
ref_select_size     = varargin{12};
Method              = varargin{13};
Tar_row             = varargin{14};
Tar_col             = varargin{15};
thresh_h            = varargin{16};
H_accuracy          = varargin{17};


temp_coord_z        = tar_select_nor_h*Scale(3) + Mean(3);
null_flag           = (ref_select_nor{3} > -1 & ref_select_nor{3} < 1 & tar_select_nor_h > -1 & tar_select_nor_h < 1);

index_1             = null_flag;
t_size              = size(ref_select_nor{1});

diff_distance       = ones(t_size(1),1)*(-9999);

if strcmp(Method,'3D')
    diff_distance(index_1)          = (ref_select_nor{3}(index_1) - tar_select_nor_h(index_1))*Scale(3);
    
    t_coord{1}                      = ref_select_nor{1};
    t_coord{2}                      = ref_select_nor{2}; 
    t_coord{3}                      = tar_select_nor_h;
    [ref_X ref_Y ref_Z select_flag] = SufaceDistance(ref_select_H,G,t_coord,RefImg_gridspace,Boundary,[ref_select_size(1) ref_select_size(2)], Mean, Scale, ref_select_nor{3});
    
    dx                              = (t_coord{1} - ref_X)*Scale(1);
    dy                              = (t_coord{2} - ref_Y)*Scale(2);
    dz                              = (t_coord{3} - ref_Z)*Scale(3);
    distance                        = sqrt(dx.*dx + dy.*dy + dz.*dz);

    tt_index                        = find(dz < 0);
    distance(tt_index)              = -distance(tt_index);
    coord_flag                      =  ref_X > -1 & ref_X < 1 & ref_Y > -1 & ref_Y < 1 & ref_Z > -1 & ref_Z < 1;        
    index_1                         = index_1 & tar_flag & select_flag & coord_flag;
    diff_distance(index_1)          = -distance(index_1);
    
    clear ref_select_H G t_coord ref_X ref_Y ref_Z select_flag dx dy dz distance;  
else
    index_1                         = index_1 & tar_flag;
    diff_distance(index_1)          = (ref_select_nor{3}(index_1) - tar_select_nor_h(index_1))*Scale(3);
end

