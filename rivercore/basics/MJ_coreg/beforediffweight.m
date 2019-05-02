function [diff_distance weight] = beforediffweight(varargin)

ref_select_nor      = varargin{1};
ref_select_nor_h    = varargin{2};
ref_select_H        = varargin{3};
tar_H_array         = varargin{4};
RefImg_gridspace    = varargin{5};
TarImg_gridspace    = varargin{6};
Boundary            = varargin{7};
Tarboundary         = varargin{8};
Mean                = varargin{9};
Scale               = varargin{10};
ref_select_size     = varargin{11};
Method              = varargin{12};
Tar_row             = varargin{13};
Tar_col             = varargin{14};
X                   = zeros(7,1);
X(1,1)              = 1.0;

[ref_Z_ori         Gr{1} Gr{2} Gr{3} Dr ref_hroh ref_adj_h ref_flag] = interpolationTIN_ref(ref_select_H,RefImg_gridspace,ref_select_nor{1}, ref_select_nor{2},Mean, Scale,Boundary,[ref_select_size(1) ref_select_size(2)],1);
[iter_tar_height_t G{1}  G{2}  G{3}  D  tar_hroh Tar_adj_h tar_flag] = interpolationTIN_Modify(tar_H_array,TarImg_gridspace,ref_select_nor{1}, ref_select_nor{2},Mean, Scale,Tarboundary,X,[Tar_row Tar_col],1);
[ref_slope ref_aspect]=SlopeAspect(Gr{1}, Gr{2}, Gr{3}, Scale);
[tar_slope tar_aspect]=SlopeAspect(G{1} , G{2} , G{3} , Scale);

pt_index            = ref_flag & tar_flag;

%similarity measurement
Rho                 = Correlation(ref_adj_h, Tar_adj_h, 3);
Ration_slope        = 1- ((abs(ref_slope - tar_slope))./(90));
Ration_aspect       = 1- ((abs(ref_aspect - tar_aspect))./360);

Wa                  = 20;
Ws                  = 40;            
Wr                  = 40;    
weight              = Wa*Ration_aspect + Ws*Ration_slope + Wr*Rho;


t_size              = size(ref_select_nor{1});

if t_size(1) > 0
    if strcmp(Method,'3D')
        t_coord{1}      = ref_select_nor{1}; t_coord{2} = ref_select_nor{2}; t_coord{3} = iter_tar_height_t;
        [ref_X ref_Y ref_Z select_flag] = SufaceDistance(ref_select_H,G,t_coord,RefImg_gridspace,Boundary,[ref_select_size(1) ref_select_size(2)], Mean, Scale, ref_Z_ori);
        diff_distance   = (ref_Z_ori - iter_tar_height_t)*Scale(3);
        dx              = (t_coord{1} - ref_X)*Scale(1);
        dy              = (t_coord{2} - ref_Y)*Scale(2);
        dz              = (t_coord{3} - ref_Z)*Scale(3);
        distance        = sqrt(dx.*dx + dy.*dy + dz.*dz);

        tt_index        = find(dz < 0);
        distance(tt_index)          = -distance(tt_index);

        diff_distance(select_flag)  = -distance(select_flag);


        null_pt         = abs(dx) < 200 & abs(dy) < 200;
        diff_distance   = diff_distance(null_pt);
        clear RefImg ref_select_H ref_Z Gr ref_flag G D tar_flag t_coord ref_X ref_Y ref_Z select_flag flag tar_hroh Tar_adj_h dx dy dz distance;  
    else
        t_coord         = Denormalize_coord(ref_select_nor,Mean,Scale,'XY');
        t_coord         = GeoToImage(t_coord,Tarboundary,TarImg_gridspace);

        [tar_pt_height_1 flag]= interpolation(tar_H_array,t_coord{1},t_coord{2},Tar_row, Tar_col);
        iter_tar_height = (tar_pt_height_1-Mean(3))/Scale(3);
        diff_distance   = (ref_select_nor_h - iter_tar_height)*Scale(3);

        diff_distance   = diff_distance(flag);

        clear tar_pt_height_1 tar_pt_height;
    end
else
    diff_distance       = zeros(0,0);
end

% null_height     = abs(diff_distance) < 1000;
% diff_distance   = diff_distance(null_height);
