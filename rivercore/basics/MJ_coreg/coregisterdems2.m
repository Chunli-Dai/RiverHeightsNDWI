function [outs] = coregisterdems2(varargin)

%% outs = [Zc, Zcs, Zref, stats],  varargin = [Zt, Zr, gridspace, params]
%% Zc, Zt, Zr structs
%  Zc : coregistered DEM
%  Zc.x : x coordinate vector (1D matrix)
%  Zc.y : y coordinate vector (1D matrix)
%  Zc.z : z coordinate vector (2D matrix)
%  Zcs  : Control surfaces used in weighted least squares adjustment
%  Zref : refernce DEM overlapped with target DEM
%  Zt   : target DEM
%  Zr   : reference DEM
%  Zcs, Zt, Zr sturcts with same format Zc

%% stats structs
%  1. P : Result on sigma0 and 7 parameters for scale, rotates(w,p,k) , shifts(tx,ty,tz)
%   stats.P[So X(7) M_scale(3) S_scale(3)]
%   where, 
%   So : sigma0 of adjustment
%   X[s w p k tx ty tz]  : 7 parameters
%     where, 
%     s  : scale
%     w  : omega angle along the x-axis (unit = deg)
%     p  : phi angle along the y-axis   (unit = deg)
%     k  : kappa angle along the z-axis (unit = deg)
%     tx : x-axis shift (unit = m)
%     ty : y-axis shift (unit = m)
%     tz : z-axis shift (unit = m)
%   M_scale[X_mean Y_mean Z_mean] : mean value for normalizing xyz coordinates between -1 and 1
%   S_scale[X_scale Y_scale Z_scale] : scale value for normalizing xyz coordinates between -1 and 1
%  2. S : statistics 
%   stats.S[iteration ncs gridspace before(5) after(5)]
%   where,
%   iteration    : the number of iteration
%   ncs          : the number of control surfaces   
%   gridspace    : gridspace of DEM
%   before[dh mean std max min] : statistics before coregistration
%   after[dh mean std max min] : statistics after coregistration
%     where,
%     dh    : height differencing value of control surfaces
%     mean  : mean of dh
%     std   : standard deviation of dh
%     max   : max of dh
%     min   : min of dh
%  3. WH : weight histogram on surface similarities
%   stats.WH[before(1D matrix) after(1D matrix)]
%   where,
%   before     : weight of control surfaces before coregistration
%   after      : weight of control surfaces after coregistration
%  4. E  : histogram on height differencing according to along and across track
%   stats.E[diff diff_std interval]
%   where,
%   diff[along(1D matrix) across(1D matrix) Z(1D matrix)]  : height difference on each direction
%   diff_std[along(1D matrix) across(1D matrix) Z(1D matrix)] : standard deviation of height differencs on each direction
%   interval[along(1D matrix) across(1D matrix) Z(1D matrix)] : histogram interval on each direction
%% gridspace : gridspace of target and reference DSM respectively
%   gridspace[target_DSM_gridspace reference_DSM_gridspace]
%% params structs : if params dose not input, all value set up default values as following.
%  1. I : Initial 7 parameters for adjustment
%   params.I[s w p k tx ty tz]
%  If these values don't be known, use as following values.
%   params.I[1 0 0 0 0 0 0] (default)
%  2. C : Stop criteria (criteria for stopping the iterative least squares adjustment)
%   params.C[MI Mds Mdr Mds MdS]
%   where, (default values)
%   MI = 100 (default): the number of maximum iteration
%   Mds = 0.001 (default): the maximum correction value of scale
%   Mdr = 0.001 (default): the maximum correction value of rotations (unit = deg)
%   Mds = 0.05  (default): the maximum correction value of shifts (unit = m)
%   MdS = 0.0001 (default): the maximum changing-rate of sigma0
%  3. G : General information
%   params.G[MH AH]
%   where,
%   MH = 3000 (default) : Maximum Height value for removing null heights such as clouds or false matching value (unit = m)
%   AH = 20   (default) : DEM height accuracy (unit = m)  
%  4. M : Method of surface representation. Select between '2.5D' and '3D'
%   params.M = 3D (default)
%  5. V : minimum similarity value for selecting control surfaces
%   params.V[weight minimum_angle minimum_undulation]
%   where,
%   weight = 10 (default): maximum percentage of weights accumulated from right end on weight histogram
%   minumum_angle = 20 (default) (unit = deg) : minimum angle of surface slope for selecting control surface
%   minumum_undulation = 10 (default) (unit = m): minimum undulation of height differences within windows mask (3X3) for selecting control surface
%%

%% loading information 
Zt                                  = varargin{1};
Zr                                  = varargin{2};
IsResultOpened                      = 'No';
RockfilterOnly                      = 'No';
if nargin == 3
    %default value
    Gridspace                       = varargin{3};
    General_info.max_height         = 3000;
    General_info.height_accuracy    = 20;
    General_info.boundary_buffer    = 0;
    General_info.point_buffer       = 0;
    General_info.height_diff_th     = 20;

    Method                          = '3D';

    Sim_param(1,1)      = 1.0;  Sim_param(2,1)      = 0.0;  Sim_param(3,1)      = 0.0;  Sim_param(4,1)      = 0.0;
    Sim_param(5,1)      = 0.0;  Sim_param(6,1)      = 0.0;  Sim_param(7,1)      = 0.0;

    Ini_param                       = Sim_param;

    Stop_condition.max_iteration    = 100;
    Stop_condition.delta_scale      = 0.001;
    Stop_condition.delta_angle      = 0.001;
    Stop_condition.delta_shift      = 0.05;
    Stop_condition.delta_sigma      = 0.0001;

    weight_percen                   = 10;    
    adjust_method                   = 'step_by_step';   
    non_featrue_value               = [20 10];
    
elseif nargin == 4
  
    Gridspace                       = varargin{3};
    params                          = varargin{4};
    
    General_info.max_height         = params.G(1);
    General_info.height_accuracy    = params.G(2);
    General_info.boundary_buffer    = 0;
    General_info.point_buffer       = 0;
    General_info.height_diff_th     = params.G(2);

    Method                          = params.M;

    Sim_param(1,1)      = 1.0;  Sim_param(2,1)      = 0.0;  Sim_param(3,1)      = 0.0;  Sim_param(4,1)      = 0.0;
    Sim_param(5,1)      = 0.0;  Sim_param(6,1)      = 0.0;  Sim_param(7,1)      = 0.0;

    for i=1:7
        Ini_param(i,1)              = params.I(i);
    end

    Stop_condition.max_iteration    = params.C(1);
    Stop_condition.delta_scale      = params.C(2);
    Stop_condition.delta_angle      = params.C(3);
    Stop_condition.delta_shift      = params.C(4);
    Stop_condition.delta_sigma      = params.C(5);


    weight_percen                   = params.V(1);    
    adjust_method                   = 'step_by_step';   
    non_featrue_value(1)            = params.V(2);
    non_featrue_value(2)            = params.V(3);
end


interpolation_method= 'TIN';

clear varargin params;
%% set up target and reference DSM and rock fildter data 
RefImg                              = Zr;
TarImg                              = Zt;

null_tar_value                      = find(TarImg.z < 1);
TarImg.z(null_tar_value)            = -9999;

% figure(1);
% imshow(RefImg.z);
% figure(2);
% imshow(TarImg.z);
clear Zr Zt;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Boundary Detect
RefImg_gridspace                    = Gridspace(2);
TarImg_gridspace                    = Gridspace(1);
gridspace                           = RefImg_gridspace;

Null_height                         = General_info.max_height;
[Boundary Refboundary Tarboundary RefHeight_min_max TarHeight_min_max] = BoundarySetting(RefImg, TarImg, Null_height, General_info.boundary_buffer,gridspace);
row_num                             = floor((Boundary(4)-Boundary(2))/gridspace)+1;     
col_num                             = floor((Boundary(3)-Boundary(1))/gridspace)+1;

%rock filter data
if 1 % chunli
tempRockz=0;
else
[tempRockz tempRockx tempRocky]     = subsetGimpIceMask(Boundary(1),Boundary(3),Boundary(2),Boundary(4));
end

if tempRockz == 0
    IsRockOpened                    = 'No';
else
    IsRockOpened                    = 'Yes';
    RockImg.rock.x                  = tempRockx;
    RockImg.rock.y                  = tempRocky;
    RockImg.rock.z                  = tempRockz;
    RockImg_gridspace               = 15; %default
end

clear tempRockz tempRockx tempRocky;

%% Main procedure
if row_num > 10 && col_num > 10
   
    if strcmp(IsRockOpened,'Yes')
        Rockboundary(1)             = min(RockImg.rock.x); 
        Rockboundary(2)             = min(RockImg.rock.y);
        Rockboundary(3)             = max(RockImg.rock.x); 
        Rockboundary(4)             = max(RockImg.rock.y);
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Normalized Mean X,Y,Z
    Mean_X                          = double(double(col_num-1)/2.0)*gridspace + Boundary(1);
    Mean_Y                          = Boundary(4) - double(double(row_num-1)/2.0)*gridspace;
    Mean_Z                          = double(max((RefHeight_min_max(2)-RefHeight_min_max(1)),(TarHeight_min_max(2)-TarHeight_min_max(1))))/2.0 + double(min(RefHeight_min_max(1),TarHeight_min_max(1)));
    Mean                            = [Mean_X Mean_Y Mean_Z];
    %Normalized Scale
    X_scale                         = col_num*gridspace/2.0;
    Y_scale                         = row_num*gridspace/2.0;
    Z_scale                         = double(max((RefHeight_min_max(2)-RefHeight_min_max(1)),(TarHeight_min_max(2)-TarHeight_min_max(1))))/2.0;
    Scale                           = [X_scale Y_scale Z_scale];
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %stop condition and Initial Value
    delta_angle                     = Stop_condition.delta_angle;
    delta_shift                     = Stop_condition.delta_shift;  
    delta_scale                     = Stop_condition.delta_scale;   
    delta_sigma                     = Stop_condition.delta_sigma;

    PI                              = 3.141592;
    DegtoRad                        = PI/180.0;
    RadtoDeg                        = 180.0/PI;
    
    img_size                        = size(TarImg.z);

    %%find top left coordinate
    stop                            = 0;
    count                           = 1;
    Top_left.row                    = 1;
    Top_left.column                 = 1;
    while stop == 0 && count < img_size(1)
         index                      = find(TarImg.z(count,:) > -1 & TarImg.z(count,:) < 5000);
         t_size                     = size(index);
         if t_size(2) > 0
             Top_left.row           = count;
             Top_left.column        = index(1);
             stop                   = 1;
         end
         count                      = count + 1;
    end
    %%find bottom left coordinate
    stop                            = 0;
    count                           = 1;
    Bottom_left.row                 = 1;
    Bottom_left.column              = 1;
    while stop == 0 && count < img_size(2)
         index                      = find(TarImg.z(:,count) > -1 & TarImg.z(:,count) < 5000);
         t_size                     = size(index);
         if t_size(1) > 0
             Bottom_left.row        = index(1);
             Bottom_left.column     = count;
             stop                   = 1;
         end
         count                      = count + 1;
    end

    bb                              = Top_left.column-Bottom_left.column;
    
    if bb == 0
        rotate_angle                = 0;
    else
        rotate_angle                = atan(abs(Top_left.row-Bottom_left.row)/abs(Top_left.column-Bottom_left.column));
    end

    % dscale, domega, dphi, dkappa, tx, ty, tz
    X                               = Ini_param;
    X(2,1)                          = X(2,1)*DegtoRad;
    X(3,1)                          = X(3,1)*DegtoRad; 
    X(4,1)                          = X(4,1)*DegtoRad;   
    X(5,1)                          = X(5,1)/Scale(1);  
    X(6,1)                          = X(6,1)/Scale(2);  
    X(7,1)                          = X(7,1)/Scale(3);
    TestX                           = Sim_param;     
    TestX(2,1)                      = TestX(2,1)*DegtoRad;
    TestX(3,1)                      = TestX(3,1)*DegtoRad;
    TestX(4,1)                      = TestX(4,1)*DegtoRad; 
    TestX(5,1)                      = TestX(5,1)/Scale(1); 
    TestX(6,1)                      = TestX(6,1)/Scale(2);
    TestX(7,1)                      = TestX(7,1)/Scale(3); 
    iter_angle                      = 100000;
    iter_shift                      = 100000;    
    iter_sigma                      = 100000;
    Sigma_pre                       = 100000;

    Max_iteration                   = Stop_condition.max_iteration;     
    thresh_h                        = General_info.height_diff_th; 
    H_accuracy                      = General_info.height_accuracy;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %conjugate point setting
    [ConjugateInfo]                 = ConjugatePointSetting(Boundary, RefImg_gridspace, Refboundary, RefImg, [X_scale Y_scale Z_scale], [Mean_X Mean_Y Mean_Z], TestX, Method);
    [Ref_row Ref_col]               = size(RefImg.z);

    ref_select                      = ConjugateInfo{1};    
    ref_select_size                 = ConjugateInfo{2};     
    ref_H_array                     = ConjugateInfo{3};  
    ref_image_index                 = ConjugateInfo{5};
    
    minx = min(ref_select{1});
    miny = min(ref_select{2});
    maxx = max(ref_select{1});
    maxt = max(ref_select{2});
    %ref image save
    Zref.x                          = ref_select{1};
    Zref.y                          = ref_select{2};
    Zref.z                          = ConjugateInfo{4};
    Zcs.x                           = ref_select{1};
    Zcs.y                           = ref_select{2};
    Zc.x                            = ref_select{1};
    Zc.y                            = ref_select{2};
        
    clear ConjugateInfo ;

    ref_select_nor_ori              = Normalize_coord(ref_select,[Mean_X Mean_Y Mean_Z], [X_scale Y_scale Z_scale],'XYZ');

    tar_select                      = ref_select;     
    image_coord_float               = GeoToImage(tar_select,Tarboundary,TarImg_gridspace);
    [Tar_row Tar_col]               = size(TarImg.z);
    tar_H_array                     = TarImg.z(:);
    null_tar                        = find(tar_H_array < 0 | tar_H_array > Null_height);
    t_size                          = size(null_tar);
    tar_H_array(null_tar)           = ones(t_size(1),1)*(-9999);
    t_size                          = size(ref_select{3});
    tar_select{3}                   = ones(t_size(1),1)*(-9999);

    ref_select_H                    = ref_select{3};
    null_ref                        = find(ref_select_H < 0 | ref_select_H > Null_height);
    t_size                          = size(null_ref);
    ref_select_H(null_ref)          = ones(t_size(1),1)*(-9999);
    clear null_tar null_ref;

    %compute heights of target surface on location of reference surface 
    if strcmp(interpolation_method,'TIN')
        [tar_height G_tar{1} G_tar{2} G_tar{3} garbage garbage garbage tar_h_flag] = interpolationTIN_ref(tar_H_array,TarImg_gridspace,ref_select_nor_ori{1}, ref_select_nor_ori{2},Mean, Scale,Tarboundary,[Tar_row Tar_col],0);
        tar_select{3}(tar_h_flag)   = tar_height(tar_h_flag);
        tar_select{3}(tar_h_flag)   = tar_select{3}(tar_h_flag) * Scale(3) + Mean(3);
        clear garbage tar_height;
    else
        G_tar                       = 0;
        [tar_select{3} tar_h_flag]  = interpolation(tar_H_array,image_coord_float{1},image_coord_float{2},Tar_row, Tar_col);
    end
    tar_select_nor_h_ori            = Normalize_coord(tar_select,Mean,Scale,'Z');

    clear image_coord_float ref_H_array;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %main process
    window_size = 3;
    if strcmp(IsResultOpened, 'No')
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %null height remove
        disp('saving before co-registration');
        [ori_diff null_flag]                = SaveBeforeData(tar_h_flag, ref_select_nor_ori, tar_select_nor_h_ori, ref_select_H, G_tar, RefImg_gridspace, TarImg_gridspace, Boundary, Tarboundary, Mean, Scale, ref_select_size, Method,Tar_row, Tar_col,thresh_h,H_accuracy);
        disp('end saving before co-registration');
        clear G_tar tar_h_flag;

        histogram_save.null_height_remove   = ori_diff(null_flag);


        %Rock and Lidar conjuate point setting
        if strcmp(IsRockOpened,'Yes')
            rock_select_ori         = FindImageFromCoord(RockImg.rock.z,Rockboundary,RockImg_gridspace,ref_select);
            clear save_overlap_array RockImg;
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %initial difference and statistic for selecting candidate corresponding surfaces
        if strcmp(IsRockOpened,'Yes')
            size_t1                 = size(ori_diff(null_flag));
            size_t2                 = size(ori_diff(rock_select_ori == 0 & null_flag));
    
            %rock area covered-ratio
            percen                  = size_t2(1) / size_t1(1) * 100;
        end

        DiffArray                   = ori_diff(null_flag);
        Std_ori                     = std(DiffArray);
        Mean_ori                    = mean(DiffArray);    

        if strcmp(IsRockOpened,'Yes')
            if percen > 1.0 && strcmp(RockfilterOnly,'Yes');
                sigma_flag          = rock_select_ori == 0;
            else
                sigma_flag          = (ori_diff > (Mean_ori - Std_ori) & ori_diff < (Mean_ori + Std_ori) | rock_select_ori == 0);
            end
        else
            sigma_flag              = (ori_diff > (Mean_ori - Std_ori) & ori_diff < (Mean_ori + Std_ori));
        end

        flag                        = null_flag & sigma_flag;

        ref_select_nor{1}           = ref_select_nor_ori{1}(flag);   
        ref_select_nor{2}           = ref_select_nor_ori{2}(flag); 
        ref_select_nor{3}           = ref_select_nor_ori{3}(flag);
        tar_select_nor_z            = tar_select_nor_h_ori(flag);

        if strcmp(IsRockOpened,'Yes')
            rock_select             = rock_select_ori(flag);
        end    
        histogram_save.initial_surfaces = ori_diff(flag);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        image                       = Denormalize_coord(ref_select_nor,Mean,Scale,'XY');
        image                       = GeoToImage(image,Boundary,RefImg_gridspace);
        index                       = (floor(image{1})-1) * ref_select_size(1) + floor(image{2});
        value                       = ori_diff(flag);
        out_array                   = savefromindex(ref_select_size,index,value,-9999);
        %save_diff_ori_img_file_name = sprintf('%s/height_difference_before_coregistration_sigma_remove.tif',save_folder);
        %geotiffwrite(save_diff_ori_img_file_name, bbox, out_array,option.bit_depth,option);

        clear tar_select ori_diff  index  out_array value RefImg TarImg ref_select flag null_flag sigma_flag; 
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


        %selecting initial corresponding surfaces
        disp('computing reference surface information');
        [Z_ref Gx_ref Gy_ref Gz_ref D_ref ref_hroh Ref_adj_h index_flag] = interpolationTIN_ref(ref_select_H,RefImg_gridspace,ref_select_nor{1}, ref_select_nor{2}, Mean, Scale, Boundary,[ref_select_size(1) ref_select_size(2)],1);
        [ref_slope ref_aspect]      = SlopeAspect(Gx_ref, Gy_ref, Gz_ref, Scale);
        clear Gx_ref Gy_ref Gz_ref D_ref;
        disp('end computing reference surface information');

        %null pt remove
        flag                        = index_flag & ref_slope > 0 & ref_hroh > 0 & index_flag;

        for i=1:3
            ref_select_nor{i}       = ref_select_nor{i}(flag);
        end
        tar_select_nor_z            = tar_select_nor_z(flag);
        ref_slope                   = ref_slope(flag);
        ref_aspect                  = ref_aspect(flag);
        ref_hroh                    = ref_hroh(flag);
        Z_ref                       = Z_ref(flag);
        for i=1:window_size*window_size
            Ref_adj_h{i}            = Ref_adj_h{i}(flag);
        end
        if strcmp(IsRockOpened,'Yes')
            rock_select             = rock_select(flag);
        end

        clear index_flag flag;


        candidate_coord{1}          = ref_select_nor{1};     
        candidate_coord{2}          = ref_select_nor{2};     
        candidate_coord{3}          = tar_select_nor_z;
        candidate_ref_coord_z       = ref_select_nor{3};
        clear ref_select_nor tar_select_nor_z;

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %adjustment
        total_size_ini              = size(candidate_coord{1});
        count = 1;     
        if strcmp(adjust_method,'step_by_step')
            nnnn                    = 2;
        else
            nnnn                    = 1;
        end

        for nn=1:nnnn

            pts_index               = true(total_size_ini(1),1);
            while ( count < Max_iteration && (iter_shift > delta_shift || iter_angle > delta_angle) && iter_sigma > delta_sigma)
                pts_index           = true(total_size_ini(1),1);
                diff_distance       = zeros(total_size_ini(1),1);
                distance            = zeros(total_size_ini(1),1);
                weight              = zeros(total_size_ini(1),1);

                %Transformation
                transformed_coord   = ConformalTransform(candidate_coord,X);

                %TarImg Gx, Gy, Gz, Coeff
                disp('computing target surface information');
                [iter_height G{1} G{2} G{3} D tar_hroh Tar_adj_h index_flag] = interpolationTIN_Modify(tar_H_array,TarImg_gridspace,transformed_coord{1}, transformed_coord{2},Mean, Scale,Tarboundary,X,[Tar_row Tar_col],1);
                disp('end  computing target surface information');
                [tar_slope tar_aspect]  = SlopeAspect(G{1}, G{2}, G{3}, Scale);
                clear D;


                %null pt remove
                t_index                 = find(~index_flag | (G{1} == 0 & G{2} == 0 & G{3} == 0) | iter_height > 1 |iter_height < -1 | transformed_coord{3} < -1 | transformed_coord{3} > 1 | transformed_coord{1} < -1 | transformed_coord{1} > 1 | transformed_coord{2} < -1 | transformed_coord{2} > 1);
                t_size                  = size(t_index);
                pts_index(t_index) = false(t_size(1),1);
                clear t_index index_flag;

                %similarity measurement
                disp('computing height correlation');
                Rho             = Correlation(Ref_adj_h, Tar_adj_h, window_size);
                Ration_slope    = 1- ((abs(ref_slope - tar_slope))./(90));
                Ration_aspect   = 1- ((abs(ref_aspect - tar_aspect))./360);
                disp('end computing height correlation');

                %conjugate pt selection
                Ta              = non_featrue_value(1);   
                Troh            = non_featrue_value(2);
                t_index         = find(tar_slope < Ta | tar_hroh < Troh);
                t_size          = size(t_index);
                pts_index(t_index) = false(t_size(1),1);
                clear t_index;

                Wa              = 20;
                Ws              = 40;
                Wr              = 40;    
                weight          = Wa*Ration_aspect + Ws*Ration_slope + Wr*Rho;

                [t_n t_xout]    = hist(weight(pts_index),[0:2:100]);
                t_sum           = sum(t_n);
                t_sum           = t_sum*(weight_percen)*0.01;
                t_size          = size(t_n);       
                tt_sum          = 0;
                ii              = 0;

                t_select        = 0;
                while tt_sum < t_sum 
                    tt_sum      = tt_sum + t_n(t_size(2) - ii);
                    t_select    = t_xout(t_size(2) - ii);
                    ii = ii + 1; 
                end


                if t_select > 100 - weight_percen
                    t_index     = find(weight < 100 - weight_percen);
                else
                    t_index     = find(weight < t_select);
                end

                t_size          = size(t_index);
                pts_index(t_index) = false(t_size(1),1);
                clear t_index;

                clear Tar_adj_h Ration_slope Ration_aspect  Rho;

                pts_pt_count    = size(pts_index);
                
                if pts_pt_count(1) > 10
                    
                    disp('computing surface distance');
                    if strcmp(Method,'3D')
                        t_coord{1}  = candidate_coord{1}; 
                        t_coord{2}  = candidate_coord{2};  
                        t_coord{3}  = iter_height;   

                        [ref_X ref_Y ref_Z select_flag] = SufaceDistance(ref_select_H,G,t_coord,RefImg_gridspace,Boundary,ref_select_size, Mean, Scale, Z_ref);
                        pts_index   = pts_index & select_flag;
                        D           = -(G{1}.*candidate_coord{1} + G{2}.*candidate_coord{2} + G{3}.*iter_height);
                        diff_distance = -(G{1}.*ref_X + G{2}.*ref_Y + G{3}.*ref_Z + D);

                        dx = (t_coord{1} - ref_X)*Scale(1);
                        dy = (t_coord{2} - ref_Y)*Scale(2);
                        dz = (t_coord{3} - ref_Z)*Scale(3);
                        distance    = sqrt(dx.*dx + dy.*dy + dz.*dz);

                        tt_index            = find(dz < 0);
                        distance(tt_index)  = -distance(tt_index);
                        clear D ref_X ref_Y ref_Z t_coord select_flag t_coord dx dy dz;
                    else
                        clear G;
                        if strcmp(interpolation_method,'TIN')
                            diff_distance   = (Z_ref - iter_height); 
                            %Gx
                            tar_col         = transformed_coord{1} + TarImg_gridspace/Scale(1);
                            tar_row         = transformed_coord{2};
                            [Gx_height garbage garbage garbage garbage garbage garbage col_index] = interpolationTIN_Modify(tar_H_array,TarImg_gridspace,tar_col,tar_row,Mean, Scale,Tarboundary,X,[Tar_row Tar_col],0);
                            G{1}            = (Gx_height - iter_height)./(TarImg_gridspace/Scale(1));
                            %Gy     
                            tar_col         = transformed_coord{1};
                            tar_row         = transformed_coord{2} + TarImg_gridspace/Scale(2);
                            [Gy_height garbage garbage garbage garbage garbage garbage row_index] = interpolationTIN_Modify(tar_H_array,TarImg_gridspace,tar_col,tar_row,Mean, Scale,Tarboundary,X,[Tar_row Tar_col],0);
                            G{2}            = (Gy_height - iter_height)./(TarImg_gridspace/Scale(2));
                            %Gz
                            t_size          = size(G{1});
                            G{3}            = ones(t_size,1)*(-1);

                            pts_index       = pts_index & col_index & row_index;
                            clear tar_col tar_row Gx_height col_index row_index garbage;
                        else
                            clear iter_height;

                            R               = RotationMatrix([X(2,1) X(3,1) X(4,1)]);
                            coord           = Denormalize_coord(transformed_coord,Mean,Scale,'XY');
                            coord           = GeoToImage(coord,Tarboundary,TarImg_gridspace);
                            [tar_pt_height_1 ori_flag]= interpolation(tar_H_array,coord{1},coord{2},Tar_row, Tar_col);
                            tar_pt_height   = (tar_pt_height_1-Mean_Z)/Z_scale;
                            iter_height     = X(1,1)*(R(1,3)*transformed_coord{1} + R(2,3)*transformed_coord{2} + R(3,3)*tar_pt_height) + X(7,1);

                            diff_distance   = (candidate_ref_coord_z - iter_height);         
                            ori_h           = tar_pt_height > -1 & tar_pt_height < 1;
                            clear tar_pt_height_1 tar_pt_height;

                            %Gx
                            tar_col         = coord{1} +1;
                            tar_row         = coord{2};
                            [temp_tar_col_height_1 gx_flag] = interpolation(tar_H_array,tar_col,tar_row,Tar_row, Tar_col);
                            temp_tar_col_height =  (temp_tar_col_height_1-Mean_Z)/Z_scale;
                            temp_tar_col_height = X(1,1)*(R(1,3)*transformed_coord{1} + R(2,3)*transformed_coord{2} + R(3,3)*temp_tar_col_height) + X(7,1);
                            G{1}            = double(temp_tar_col_height - iter_height)/(TarImg_gridspace/X_scale);
                            gx_h            = temp_tar_col_height > -1 & temp_tar_col_height < 1;
                            clear tar_col tar_row temp_tar_col_height_1 temp_tar_col_height;

                            %Gy
                            tar_col         = coord{1};
                            tar_row         = coord{2} -1;
                            [temp_tar_row_height_1 gy_flag]= interpolation(tar_H_array,tar_col,tar_row,Tar_row, Tar_col);
                            temp_tar_row_height = (temp_tar_row_height_1-Mean_Z)/Z_scale;
                            temp_tar_row_height = X(1,1)*(R(1,3)*transformed_coord{1} + R(2,3)*transformed_coord{2} + R(3,3)*temp_tar_row_height) + X(7,1);
                            G{2}            = double(temp_tar_row_height - iter_height)/(TarImg_gridspace/Y_scale);
                            gy_h            = temp_tar_row_height > -1 & temp_tar_row_height < 1;
                            clear tar_col tar_row temp_tar_row_height_1 temp_tar_row_height;

                            t_size          = size(G{1});
                            G{3}            = ones(t_size,1)*(-1);

                            pts_index       = pts_index & ori_flag & gx_flag & gy_flag & gx_h & gy_h & ori_h;
                            clear ori_flag gx_flag gy_flag coord gx_h gy_h ori_h;
                        end
                    end
                    disp('end computing surface distance');

                    %weighted-least square adjustment
                    Wb              = zeros(7,7);            
                    if (strcmp(adjust_method,'shift_only') || strcmp(adjust_method,'step_by_step')) && nn == 1
                        Wb(1,1)     = 10000000000000000000000000;          
                        Wb(2,2)     = 10000000000000000000000000;  
                    else
                        Wb(1,1)     = 1;
                        Wb(2,2)     = 1;
                    end

                    Wb(3,3)         = Wb(2,2);  
                    Wb(4,4)         = Wb(2,2);
                    Wb(5,5)         = 1;
                    Wb(6,6)         = 1;   
                    Wb(7,7)         = 1;
                    Lb              = ones(7,1)*0.000001;

                    [A AW L]        = CoeffMatrix_25D(transformed_coord{1}(pts_index),transformed_coord{2}(pts_index),transformed_coord{3}(pts_index),G{1}(pts_index), G{2}(pts_index), -G{3}(pts_index),[X(2,1) X(3,1) X(4,1)],X(1,1),diff_distance(pts_index),weight(pts_index));
                    Qxx             = inv(AW.'*A + Wb);     
                    X0              = Qxx*(AW.'*L + Wb*Lb);      
                    X0              = real(X0);
                    V               = A*X0 - L;   
                    Vb              = X0 - Lb;           
                    VW              = V.*weight(pts_index);
                    VWV             = V.'*(V.*weight(pts_index));  
                    VbWbVb          = Vb.'*Wb*Vb;       
                    total_size      = size(V); 
                    sigma0          = sqrt((V.'*V)/(total_size(1) - 7));

                    clear AW A L Vb;

                    %remove blunder
                    find_index = find(pts_index);
                    if strcmp(Method,'3D')
                        sigma_th_index  = find(abs(V) > 3.29*sigma0 | abs(distance(pts_index)) > 50);
                    else
                        sigma_th_index  = find(abs(V) > 3.29*sigma0 | abs(diff_distance(pts_index)*Scale(3)) > 50);
                    end
                    t_size          = size(sigma_th_index);
                    t_index         = find_index(sigma_th_index);
                    pts_index(t_index) = false(t_size(1),1);

                    [A AW L]        = CoeffMatrix_25D(transformed_coord{1}(pts_index),transformed_coord{2}(pts_index),transformed_coord{3}(pts_index),G{1}(pts_index), G{2}(pts_index), -G{3}(pts_index),[X(2,1) X(3,1) X(4,1)],X(1,1),diff_distance(pts_index),weight(pts_index));
                    Qxx             = inv(AW.'*A + Wb);    
                    X0              = Qxx*(AW.'*L + Wb*Lb);    
                    X0              = real(X0);
                    V               = A*X0 - L;          
                    Vb              = X0 - Lb;           
                    VWV             = V.'*(V.*weight(pts_index));  
                    VbWbVb          = Vb.'*Wb*Vb;
                    total_size      = size(V); 
                    sigma0          = sqrt((V.'*V)/(total_size(1) - 7));
                    clear G;

                    if strcmp(Method,'3D')
                       diff_distance= -distance./Scale(3);
                    end

                    clear AW A L V Vb;

                    clear tar_slope tar_aspect tar_hroh;
                    

                    sigmaX          = sigma0*sqrt(Qxx);
                    sigmaX(2,2)     = sigmaX(2,2)*RadtoDeg;   
                    sigmaX(3,3)     = sigmaX(3,3)*RadtoDeg;     
                    sigmaX(4,4)     = sigmaX(4,4)*RadtoDeg;
                    sigmaX(5,5) 	= sigmaX(5,5)*X_scale;      
                    sigmaX(6,6)     = sigmaX(6,6)*Y_scale;    
                    sigmaX(7,7)     = sigmaX(7,7)*Z_scale;

                    iter_sigma      = abs(Sigma_pre-(sigma0))/Sigma_pre;        
                    Sigma_pre       = sigma0;

                    temp_matrix(1)  = abs(X0(2,1));     
                    temp_matrix(2)  = abs(X0(3,1));     
                    temp_matrix(3)  = abs(X0(4,1));
                    iter_angle      = max(temp_matrix);      
                    iter_angle      = iter_angle*RadtoDeg;

                    temp_matrix(1)  = abs(X0(5,1))*X_scale;            
                    temp_matrix(2)  = abs(X0(6,1))*Y_scale;         
                    temp_matrix(3)  = abs(X0(7,1))*Z_scale;
                    iter_shift      = max(temp_matrix);

                    if iter_shift < 0.001;
                        iter_sigma  = 0.000001;
                    end

                   
                    X = X + X0;

                    if strcmp(Method,'3D')
                       sigma0       = sqrt((diff_distance(pts_index).'*diff_distance(pts_index))/(total_size(1)-7));
                    end
                    str             = sprintf('iter=%d, S0=%f, ds=%f, dw=%f, dp=%f, dk=%f, dtx=%f, dty=%f, dtz=%f', count, sigma0*Z_scale, X0(1,1),X0(2,1)*RadtoDeg,X0(3,1)*RadtoDeg,X0(4,1)*RadtoDeg,X0(5,1)*X_scale,X0(6,1)*Y_scale,X0(7,1)*Z_scale);
                    disp(str);
                    str             = sprintf('iter=%d, S0=%f, s=%f, w=%f, p=%f, k=%f, tx=%f, ty=%f, tz=%f', count, sigma0*Z_scale, X(1,1),X(2,1)*RadtoDeg,X(3,1)*RadtoDeg,X(4,1)*RadtoDeg,X(5,1)*X_scale,X(6,1)*Y_scale,X(7,1)*Z_scale);
                    disp(str);

                    count = count +1;

                    if nn==2
                        ref_slope       = ref_slope(pts_index);
                        ref_aspect      = ref_aspect(pts_index);
                        ref_hroh        = ref_hroh(pts_index);
                        Z_ref           = Z_ref(pts_index);
                        for i=1:window_size*window_size
                            Ref_adj_h{i}= Ref_adj_h{i}(pts_index);
                        end
                        if strcmp(IsRockOpened,'Yes')
                            rock_select = rock_select(pts_index);
                        end
                        candidate_coord{1}  = candidate_coord{1}(pts_index);     
                        candidate_coord{2}  = candidate_coord{2}(pts_index);   
                        candidate_coord{3}  = candidate_coord{3}(pts_index);
                        candidate_ref_coord_z = candidate_ref_coord_z(pts_index);

                        total_size_ini = size(candidate_coord{1});
                    end
                else
                    count       = Max_iteration;
                    sigma0      = 0;
                    sigmaX      = 0;
                    Qxx         = 0;
                end
                
            end

            if nn == 1 && strcmp(adjust_method,'step_by_step')
                iter_angle      = 100000;     
                iter_shift      = 100000;    
                iter_sigma      = 100000;     
                Sigma_pre       = 100000;

                ref_slope       = ref_slope(pts_index);
                ref_aspect      = ref_aspect(pts_index);
                ref_hroh        = ref_hroh(pts_index);
                Z_ref           = Z_ref(pts_index);
                for i=1:window_size*window_size
                    Ref_adj_h{i}= Ref_adj_h{i}(pts_index);
                end
                if strcmp(IsRockOpened,'Yes')
                    rock_select = rock_select(pts_index);
                end
            end

            if nn == 1
                candidate_coord{1}  = candidate_coord{1}(pts_index);     
                candidate_coord{2}  = candidate_coord{2}(pts_index);   
                candidate_coord{3}  = candidate_coord{3}(pts_index);
                candidate_ref_coord_z = candidate_ref_coord_z(pts_index);

                total_size_ini = size(candidate_coord{1});
            end

        end
        clear Ref_adj_h ref_slope ref_aspect Z_ref ref_hroh ref_select_nor tar_select_nor_z ;
        disp('end computing parameters');

        [diff_distance_first weight_first] = beforediffweight(candidate_coord, candidate_ref_coord_z, ref_select_H, tar_H_array, RefImg_gridspace, TarImg_gridspace, Boundary, Tarboundary, Mean, Scale, ref_select_size, Method,Tar_row, Tar_col);
        histogram_save.weight_before = weight_first;

        parameters_save.sigma0      = sigma0*Scale(3);
        parameters_save.X(1,1)      = X(1,1);
        parameters_save.X(2,1)      = X(2,1)*RadtoDeg;
        parameters_save.X(3,1)      = X(3,1)*RadtoDeg;
        parameters_save.X(4,1)      = X(4,1)*RadtoDeg;
        parameters_save.X(5,1)      = X(5,1)*Scale(1);
        parameters_save.X(6,1)      = X(6,1)*Scale(2);
        parameters_save.X(7,1)      = X(7,1)*Scale(3);
        parameters_save.Qxx         = Qxx;
        parameters_save.sigmaX      = sigmaX;
        parameters_save.Scaling.Mean= Mean;
        parameters_save.Scaling.Scale = Scale;
        parameters_save.boundary    = Boundary;

        histogram_save.weight_after = weight(pts_index);

        t_size              = size(transformed_coord{1}(pts_index));
        if t_size(1) > 10
            X_space             = 5*RefImg_gridspace/Scale(1);
            Y_space             = 5*RefImg_gridspace/Scale(2);
            Z_space             = H_accuracy/Scale(3);
            %inclined axis
            rotated_coord_x     =  transformed_coord{1}(pts_index)*cos(rotate_angle) + transformed_coord{2}(pts_index)*sin(rotate_angle);
            rotated_coord_y     = -transformed_coord{1}(pts_index)*sin(rotate_angle) + transformed_coord{2}(pts_index)*cos(rotate_angle);

            rotated_coord_z     =  transformed_coord{3}(pts_index);
            dh                  = diff_distance(pts_index)*Z_scale;
            x_min               = min(rotated_coord_x);
            x_max               = max(rotated_coord_x);
            y_min               = min(rotated_coord_y);
            y_max               = max(rotated_coord_y);
            z_min               = min(rotated_coord_z);
            z_max               = max(rotated_coord_z);


            x_index             = fix((rotated_coord_x - x_min)./X_space) + 1;
            y_index             = fix((rotated_coord_y - y_min)./Y_space) + 1;
            z_index             = fix((rotated_coord_z - z_min)./Z_space) + 1;


            % '2' is a range from -1 and 1
            tx_range            = abs(x_max - x_min);
            ty_range            = abs(y_max - y_min);
            tz_range            = abs(z_max - z_min);
            x_size              = floor(tx_range/X_space)+1;
            y_size              = floor(ty_range/Y_space)+1;
            z_size              = floor(tz_range/Z_space)+1;

            % coordinate ranges from 0 and 2 => x, y axis
            X_diff              = zeros(x_size,1); 
            Y_diff              = zeros(y_size,1);
            Z_diff              = zeros(z_size,1);
            X_diff_count        = zeros(x_size,1); 
            Y_diff_count        = zeros(y_size,1);
            Z_diff_count        = zeros(z_size,1);
            X_diff_std          = zeros(x_size,1); 
            Y_diff_std          = zeros(y_size,1);
            Z_diff_std          = zeros(z_size,1);

            pt_cell             = cell(1,4);
            pt_cell{1,1}        = x_index;
            pt_cell{1,2}        = y_index;
            pt_cell{1,3}        = z_index;
            pt_cell{1,4}        = dh;

            for ii=1:max(x_index)
                temp            = find(pt_cell{1,1} == ii);
                X_diff(ii)      = mean(pt_cell{1,4}(temp));
                X_diff_std(ii)  = std(pt_cell{1,4}(temp));
            end

            for ii=1:max(y_index)
                temp            = find(pt_cell{1,2} == ii);
                Y_diff(ii)      = mean(pt_cell{1,4}(temp));
                Y_diff_std(ii)  = std(pt_cell{1,4}(temp));
            end

            for ii=1:max(z_index)
                temp            = find(pt_cell{1,3} == ii);
                Z_diff(ii)      = mean(pt_cell{1,4}(temp));
                Z_diff_std(ii)  = std(pt_cell{1,4}(temp));
            end

            X_interval          = (([1:x_size] - 1)*X_space )*Scale(1)/1000;
            Y_interval          = (([1:y_size] - 1)*Y_space )*Scale(2)/1000;
            Z_interval          = (([1:z_size] - 1)*Z_space + z_min)*Scale(3) + Mean(3);
            errorbar_save.X_diff= X_diff;
            errorbar_save.Y_diff= Y_diff;
            errorbar_save.Z_diff= Z_diff;
            errorbar_save.X_diff_std= X_diff_std;
            errorbar_save.Y_diff_std= Y_diff_std;
            errorbar_save.Z_diff_std= Z_diff_std;
            errorbar_save.X_interval=X_interval;
            errorbar_save.Y_interval=Y_interval;
            errorbar_save.Z_interval=Z_interval;
        else
            errorbar_save.X_diff= 0;
            errorbar_save.Y_diff= 0;
            errorbar_save.Z_diff= 0;
            errorbar_save.X_diff_std= 0;
            errorbar_save.Y_diff_std= 0;
            errorbar_save.Z_diff_std= 0;
            errorbar_save.X_interval= 0;
            errorbar_save.Y_interval= 0;
            errorbar_save.Z_interval= 0;
        end
        
        clear X_direction_diff Y_direction_diff Z_direction_diff x_index y_index z_index X_diff_count Y_diff_count Z_diff_count X_diff Y_diff Z_diff;
        %save selected rotated-pts in target img
        denormal            = Denormalize_coord(candidate_coord,Mean,Scale,'XY');
        image_coord         = GeoToImage(denormal,Boundary,RefImg_gridspace);
        index               = (floor(image_coord{1})-1)*ref_select_size(1) + floor(image_coord{2});
        value               = 0;
        out_array           = savefromindex(ref_select_size,index,value,1);
        Zcs.z               = out_array;
        %geotiffwrite(save_overlap_pt_img_file_name, bbox, out_array,option.bit_depth,option);
        clear index value out_array tar_select_h image_coord result_coord denormal;

        DiffArray           = diff_distance(pts_index)*Z_scale;
        Mean_diff_H         = mean(DiffArray);
        std_diff_H          = std(DiffArray);
        Max_diff            = max(DiffArray);
        Min_diff            = min(DiffArray);
        statistic_save.after.height_diff    = DiffArray;
        statistic_save.after.height_mean    = Mean_diff_H;
        statistic_save.after.height_std     = std_diff_H;
        statistic_save.after.height_max     = Max_diff;
        statistic_save.after.height_min     = Min_diff;
        statistic_save.after.pt_counts      = size(DiffArray);
        statistic_save.after.weight_ave     = mean(histogram_save.weight_after);
        statistic_save.after.weight_min     = min(histogram_save.weight_after);
        statistic_save.after.weight_max     = max(histogram_save.weight_after);

        clear DiffArray;

        DiffArray           = diff_distance_first;
        Ori_Mean_diff_H     = mean(DiffArray);
        Ori_Std_diff_H      = std(DiffArray); 
        Ori_max             = max(DiffArray);
        Ori_min             = min(DiffArray);
        statistic_save.before.height_diff   = DiffArray;
        statistic_save.before.height_mean   = Ori_Mean_diff_H;
        statistic_save.before.height_std    = Ori_Std_diff_H;
        statistic_save.before.height_max    = Ori_max;
        statistic_save.before.height_min    = Ori_min;
        statistic_save.before.pt_counts     = size(DiffArray);
        statistic_save.before.weight_ave    = mean(histogram_save.weight_before);
        statistic_save.before.weight_min    = min(histogram_save.weight_before);
        statistic_save.before.weight_max    = max(histogram_save.weight_before);
        clear DiffArray diff_distance diff_distance_first candidate_ref_coord_z candidate_coord iter_height rock_select rock_select_ori transformed_coord pts_index;

        disp('saving result');
        Zc.z            = SaveResult(ref_select_nor_ori, tar_select_nor_h_ori, ref_select_H, tar_H_array, RefImg_gridspace, TarImg_gridspace, Boundary, Tarboundary, Mean, Scale, ref_select_size, X, Method,Tar_row, Tar_col,thresh_h,H_accuracy);
        
        stats.P{1}      = parameters_save.sigma0;
        stats.P{2}      = parameters_save.X;
        stats.P{3}      = [parameters_save.Scaling.Mean(1) parameters_save.Scaling.Mean(2) parameters_save.Scaling.Mean(3)];        
        stats.P{4}      = [parameters_save.Scaling.Scale(1) parameters_save.Scaling.Scale(2) parameters_save.Scaling.Scale(3)];
        stats.S{1}      = count-1;
        stats.S{2}      = statistic_save.before.pt_counts;
        stats.S{3}      = gridspace;
        stats.S{4}      = statistic_save.before;
        stats.S{5}      = statistic_save.after;
        stats.WH{1}     = histogram_save.weight_before;
        stats.WH{2}     = histogram_save.weight_after;
        
        diff{1}         = errorbar_save.X_diff;
        diff{2}         = errorbar_save.Y_diff;
        diff{3}         = errorbar_save.Z_diff;
        std_e{1}        = errorbar_save.X_diff_std;
        std_e{2}        = errorbar_save.Y_diff_std;
        std_e{3}        = errorbar_save.Z_diff_std;
        interval{1}     = errorbar_save.X_interval;
        interval{2}     = errorbar_save.Y_interval;
        interval{3}     = errorbar_save.Z_interval;
        stats.E         = [diff std_e interval];     
        
        Zref.x          = Zref.x - gridspace;
        Zref.y          = Zref.y + gridspace;
        Zcs.x           = Zref.x;
        Zcs.y           = Zref.y;
        Zc.x            = Zref.x;
        Zc.y            = Zref.y;
    
        outs{1}         = Zc;
        outs{2}         = Zcs;
        outs{3}         = Zref;
        outs{4}         = stats;
        clear diff std_e interval;
    else
        X(1,1) = parameters_save.X(1,1);                 X(2,1) = parameters_save.X(2,1)*DegtoRad;                 X(3,1) = parameters_save.X(3,1)*DegtoRad;
        X(4,1) = parameters_save.X(4,1)*DegtoRad;        X(5,1) = parameters_save.X(5,1)/Scale(1);        X(6,1) = parameters_save.X(6,1)/Scale(2);
        X(7,1) = parameters_save.X(7,1)/Scale(3);                                                   sigma0 = parameters_save.sigma0/Scale(3);

        Zc.z            = SaveResult(ref_select_nor_ori, tar_select_nor_h_ori, ref_select_H, tar_H_array, RefImg_gridspace, TarImg_gridspace, Boundary, Tarboundary, Mean, Scale, ref_select_size, X, Method,Tar_row, Tar_col,thresh_h,H_accuracy);
        outs{1} = 0;
        outs{2} = 0;
        outs{3} = 0;
        outs{4} = 0;
    end
    disp('end processing');
else
    outs{1} = 0;
    outs{2} = 0;
    outs{3} = 0;
    outs{4} = 0;
end
