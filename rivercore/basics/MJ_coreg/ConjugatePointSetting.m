function [out] = ConjugatePointSetting(varargin)

Boundary            = varargin{1};
RefImg_gridspace    = varargin{2};
Refboundary         = varargin{3};
RefImg              = varargin{4};
Scale               = varargin{5};
Mean                = varargin{6};
test_X              = varargin{7};
method              = varargin{8};

ref_start_f(1)      = (Boundary(1)   - Refboundary(1) )/RefImg_gridspace +1;
ref_start_f(2)      = (Refboundary(4)- Boundary(4))    /RefImg_gridspace +1;
ref_end_f(1)        = (Boundary(3)   - Refboundary(1)) /RefImg_gridspace +1;
ref_end_f(2)        = (Refboundary(4)- Boundary(2))    /RefImg_gridspace +1;

ref_start(1)        = floor(ref_start_f(1));
ref_start(2)        = floor(ref_start_f(2));
ref_end(1)          = floor(ref_end_f(1));
ref_end(2)          = floor(ref_end_f(2));

RefArray            = RefImg.z(ref_start(2):ref_end(2),ref_start(1):ref_end(1));
[Ref_row Ref_col]   = size(RefImg.z);
ref_H_array         = RefImg.z(:);

[select_row select_col] = find(RefArray>-99999999.);
ref_select_size     = size(RefArray);

temp_height         = RefArray(:);
ref_select{3}       = double(temp_height);

ref_select_col      = (select_col - 1) + ref_start(1);
ref_select_row      = (select_row - 1) + ref_start(2);

Image_index         = (ref_select_col-1)*Ref_row + ref_select_row;

image{1}            = ref_select_col;
image{2}            = ref_select_row;

image               = ImageToGeo(image,Refboundary,RefImg_gridspace);
ref_select{1}       = image{1};
ref_select{2}       = image{2};
clear image;

if test_X(1,1) ~= 1 || test_X(2,1) ~= 0 || test_X(3,1) ~= 0 || test_X(4,1) ~= 0 || test_X(5,1) ~= 0 || test_X(6,1) ~= 0 || test_X(7,1) ~= 0
    height          = Normalize_coord(ref_select,Mean,Scale,'XYZ');

    R               = RotationMatrix([test_X(2,1) test_X(3,1) test_X(4,1)]);
    nor{1}          = test_X(1,1)*(R(1,1)*height{1}+R(2,1)*height{2} + R(3,1)*height{3}) + test_X(5,1);
    nor{2}          = test_X(1,1)*(R(1,2)*height{1}+R(2,2)*height{2} + R(3,2)*height{3}) + test_X(6,1);
    nor{3}          = test_X(1,1)*(R(1,3)*height{1}+R(2,3)*height{2} + R(3,3)*height{3}) + test_X(7,1);
    select_index    = nor{1} >= -1 & nor{1} <= 1 & nor{2} >= -1 & nor{2} <= 1 & nor{3} >= -1; 
    if strcmp(method,'TIN')
        [iter_height G{1} G{2} G{3} D tar_hroh Tar_adj_h index_flag] = interpolationTIN_Modify(ref_H_array,RefImg_gridspace,nor{1}, nor{2},Mean, Scale,Refboundary,test_X,[Ref_row Ref_col],0);
        t_size = size(nor{1});
        ref_select{3}                           = ones(t_size(1),1)*(-9999);
        
        ref_select{3}(select_index & index_flag)= iter_height(select_index & index_flag);
        ref_select{3}(select_index & index_flag)= ref_select{3}(select_index & index_flag) * Scale(3) + Mean(3);
    else
        geo         = Denormalize_coord(nor,Mean,Scale,'XYZ');
        image       = GeoToImage(geo,Refboundary,RefImg_gridspace);

        height_col  = image{1};
        height_row  = image{2};

        t_size      = size(height_col);
        ref_select{3}= ones(t_size(1),1)*(-9999);
        
        [ref_select_t select_flag]              = interpolation(ref_H_array,height_col,height_row,Ref_row, Ref_col);
        ref_select_t= (ref_select_t-Mean(3))/Scale(3);
        ref_select_t= test_X(1,1)*(R(1,3)*nor{1}+R(2,3)*nor{2} + R(3,3)*ref_select_t) + test_X(7,1);
        ref_select_t= ref_select_t*Scale(3) + Mean(3);
        ref_select{3}(select_flag & select_index)= ref_select_t(select_flag  & select_index);
    end
end

% t_size = size(ref_select{3});
% noise = random('norm',0,20,1,t_size(1));
% ref_select{3} = ref_select{3} + noise';
% AA = mean(noise);
% BB = std(noise);

save_ref_ori_array_t= ones(ref_select_size(1)*ref_select_size(2),1)*(-9999);
index               = (floor(select_col-1))*ref_select_size(1) + floor(select_row);
save_ref_ori_array_t(index) = ref_select{3};
save_ref_ori_array  = reshape(save_ref_ori_array_t,ref_select_size(1),ref_select_size(2));

out{1}              = ref_select;
out{2}              = ref_select_size;
out{3}              = ref_H_array;
out{4}              = save_ref_ori_array;
out{5}              = Image_index;

