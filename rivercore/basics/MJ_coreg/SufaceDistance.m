function [X Y Z select_flag]=SufaceDistance(varargin)

refer_Z         = varargin{1};
target_vector   = varargin{2};
target_point    = varargin{3};
ref_grid_space  = varargin{4};
ref_boundary    = varargin{5};
ref_img_size    = varargin{6};
Mean            = varargin{7};
Scale           = varargin{8};
p_ref_z         = varargin{9};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%initializing
total_pts           = size(target_point{1});
pts_flag            = true(total_pts(1),1);
select_flag         = false(total_pts(1),1);
X                   = ones(total_pts(1),1)*(-9999);
Y                   = ones(total_pts(1),1)*(-9999);
Z                   = ones(total_pts(1),1)*(-9999);
Z_ref               = ones(total_pts(1),1)*(-9999);
diff                = ones(total_pts(1),1)*(9999);
line_node_pt_after{1} = ones(total_pts(1),1)*(-9999);
line_node_pt_after{2} = ones(total_pts(1),1)*(-9999);
line_node_pt_after{3} = ones(total_pts(1),1)*(-9999);
mag                 = zeros(total_pts(1),1); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%denomalizing normal vector
target_vector{1}    = target_vector{1}/Scale(1);
target_vector{2}    = target_vector{2}/Scale(2);
target_vector{3}    = target_vector{3}/Scale(3);

mag(:)              = sqrt(target_vector{1}(:).*target_vector{1}(:) + target_vector{2}(:).*target_vector{2}(:) + target_vector{3}(:).*target_vector{3}(:));
target_vector{1}(:) = target_vector{1}(:)./mag(:);
target_vector{2}(:) = target_vector{2}(:)./mag(:);
target_vector{3}(:) = target_vector{3}(:)./mag(:);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%vector slope about x and y directions
vector_slope_x      = atan2(abs(target_vector{3}),abs(target_vector{1}));
vector_slope_y      = atan2(abs(target_vector{3}),abs(target_vector{2}));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%set up for start points
line_node_pt        = Denormalize_coord(target_point,Mean,Scale,'XYZ');
pre_pt_on_surface   = line_node_pt;
pre_pt_on_surface{3}= p_ref_z*Scale(3) + Mean(3);

after_pt_on_surface = pre_pt_on_surface; 
max_iter            = 50;
stop_number         = total_pts(1);
count = 1;    


while count < max_iter && stop_number > 0
    line_node_pt_after{1}(pts_flag) = line_node_pt{1}(pts_flag) + (pre_pt_on_surface{3}(pts_flag) - line_node_pt{3}(pts_flag)).*(target_vector{1}(pts_flag)./target_vector{3}(pts_flag));
    line_node_pt_after{2}(pts_flag) = line_node_pt{2}(pts_flag) + (pre_pt_on_surface{3}(pts_flag) - line_node_pt{3}(pts_flag)).*(target_vector{2}(pts_flag)./target_vector{3}(pts_flag));
    line_node_pt_after{3}(pts_flag) = pre_pt_on_surface{3}(pts_flag);
    
    find_pt{1}                      = line_node_pt_after{1}(pts_flag);
    find_pt{2}                      = line_node_pt_after{2}(pts_flag);
    find_pt                         = Normalize_coord(find_pt,Mean,Scale,'XY');
    [Z_ref tt tt tt tt tt tt index_flag] = interpolationTIN_ref(refer_Z,ref_grid_space,find_pt{1}, find_pt{2}, Mean, Scale, ref_boundary,[ref_img_size(1) ref_img_size(2)],0);
    clear tt find_pt;
    
    Z_ref                           = Z_ref*Scale(3) + Mean(3);

    selected_index                  = find(pts_flag);
    after_pt_on_surface{1}(selected_index(index_flag))  = line_node_pt_after{1}(selected_index(index_flag));
    after_pt_on_surface{2}(selected_index(index_flag))  = line_node_pt_after{2}(selected_index(index_flag));
    after_pt_on_surface{3}(selected_index(index_flag))  = Z_ref(index_flag);
    clear index_flag;
    
    dx                              = pre_pt_on_surface{1}(pts_flag) - after_pt_on_surface{1}(pts_flag);
    dy                              = pre_pt_on_surface{2}(pts_flag) - after_pt_on_surface{2}(pts_flag);
    dz                              = pre_pt_on_surface{3}(pts_flag) - after_pt_on_surface{3}(pts_flag);
    diff(selected_index)            = sqrt(dx.*dx + dy.*dy + dz.*dz);
    
    stop_pt                         = find(abs(diff) < 0.1 & pts_flag);
    X(stop_pt)                      = (after_pt_on_surface{1}(stop_pt)-Mean(1))./Scale(1);
    Y(stop_pt)                      = (after_pt_on_surface{2}(stop_pt)-Mean(2))./Scale(2);
    Z(stop_pt)                      = (after_pt_on_surface{3}(stop_pt)-Mean(3))./Scale(3);
    
    pre_angle_x                     = atan2(abs(dz),abs(dx));
    pre_angle_y                     = atan2(abs(dz),abs(dy));
    divergent_index_x               = (pre_angle_x >= vector_slope_x(pts_flag) & target_vector{1}(pts_flag) ~= 0 & dx ~= 0);
    divergent_index_y               = (pre_angle_y >= vector_slope_y(pts_flag) & target_vector{2}(pts_flag) ~= 0 & dy ~= 0);
    divergent_index                 = divergent_index_x | divergent_index_y;
    selected_index                  = selected_index(divergent_index);
    t_size                          = size(selected_index);
    
    if t_size(1) > 0
        temp_pt{1}                  = pre_pt_on_surface{1}(selected_index) + (after_pt_on_surface{1}(selected_index) - pre_pt_on_surface{1}(selected_index))./2;                
        temp_pt{2}                  = pre_pt_on_surface{2}(selected_index) + (after_pt_on_surface{2}(selected_index) - pre_pt_on_surface{2}(selected_index))./2;
        
        find_pt                     = Normalize_coord(temp_pt,Mean,Scale,'XY');
        [temp_Z tt tt tt tt tt tt temp_flag] = interpolationTIN_ref(refer_Z,ref_grid_space,find_pt{1}, find_pt{2}, Mean, Scale, ref_boundary,[ref_img_size(1) ref_img_size(2)],0);
        clear tt;
        temp_Z                                  = temp_Z*Scale(3) + Mean(3);
        after_pt_on_surface{1}(selected_index)  = temp_pt{1};
        after_pt_on_surface{2}(selected_index)  = temp_pt{2};
        after_pt_on_surface{3}(selected_index)  = temp_Z;
        clear temp_pt temp_Z find_pt temp_flag;
    end
    
    pre_pt_on_surface               = after_pt_on_surface;
    
    t_size                          = size(stop_pt);
    pts_flag(stop_pt)               = false(t_size(1),1);
    select_flag(stop_pt)            = true(t_size(1),1);
    stop_index                      = find(pts_flag);
    t_size                          = size(stop_index);
    stop_number                     = t_size(1);
    
    
    count                           = count + 1;
    
    clear dx dy dz selected_index divergent_index_x divergent_index_y divergent_index pre_angle_x pre_angle_y;
end


    
    