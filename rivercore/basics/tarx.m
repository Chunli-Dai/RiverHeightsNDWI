function [outx]=tarx(out,TarImg,RefImg,Gridspace,params)
% given conformal transformation parameters to get transformed DEM 
% out is the output of coregisterdems2.m
% outx is the output of this function, 
%   outx{2} and outx{4} are the same with that in out.

RefImg_gridspace                    = Gridspace(2);
TarImg_gridspace                    = Gridspace(1);
gridspace                           = RefImg_gridspace;
% gridspace                           = out{1,4}.S{1,3};
% TarImg_gridspace=gridspace;RefImg_gridspace=gridspace;
X=out{1,4}.P{1,2};
Method=params.M ;
Sim_param(1,1)      = 1.0;  Sim_param(2,1)      = 0.0;  Sim_param(3,1)      = 0.0;  Sim_param(4,1)      = 0.0;
Sim_param(5,1)      = 0.0;  Sim_param(6,1)      = 0.0;  Sim_param(7,1)      = 0.0;
General_info.height_accuracy    = params.G(2);
General_info.height_diff_th     = params.G(2);
General_info.max_height         = params.G(1);
%Mean=out{1,4}.P{1,3};Scale=out{1,4}.P{1,4};

PI                              = 3.141592;
DegtoRad                        = PI/180.0;
RadtoDeg                        = 180.0/PI;
null_tar_value                      = find(TarImg.z < 1);
TarImg.z(null_tar_value)            = -9999;

General_info.boundary_buffer    = 0;
Null_height                         = General_info.max_height;

[Boundary Refboundary Tarboundary RefHeight_min_max TarHeight_min_max] = BoundarySetting(RefImg, TarImg, Null_height, General_info.boundary_buffer,gridspace);

row_num                             = floor((Boundary(4)-Boundary(2))/gridspace)+1;     
col_num                             = floor((Boundary(3)-Boundary(1))/gridspace)+1;

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

X(2,1)      = X(2,1)/RadtoDeg;
X(3,1)      = X(3,1)/RadtoDeg;
X(4,1)      = X(4,1)/RadtoDeg;
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
%conjugate point setting
[ConjugateInfo]                 = ConjugatePointSetting(Boundary, RefImg_gridspace, Refboundary, RefImg, [X_scale Y_scale Z_scale], [Mean_X Mean_Y Mean_Z], TestX, Method);
[Ref_row Ref_col]               = size(RefImg.z);

ref_select                      = ConjugateInfo{1};    
ref_select_size                 = ConjugateInfo{2}; 

Zref.x                          = ref_select{1};
Zref.y                          = ref_select{2};
Zref.z                          = ConjugateInfo{4};

% ref_select={out{1,3}.x,out{1,3}.y,out{1,3}.z(:)}; 
% % corresponding to line 964 in coregisterdems2.m
% ref_select{1}=ref_select{1}+ gridspace;ref_select{2}=ref_select{2}- gridspace;
% ref_select_size=size(out{1,3}.z);
ref_select_nor_ori = Normalize_coord(ref_select,[Mean_X Mean_Y Mean_Z], [X_scale Y_scale Z_scale],'XYZ');

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
interpolation_method= 'TIN';
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



thresh_h                        = General_info.height_diff_th; 
H_accuracy                      = General_info.height_accuracy;
   
% disp('saving result');
Zc.z  =SaveResult(ref_select_nor_ori, tar_select_nor_h_ori, ref_select_H, tar_H_array, RefImg_gridspace, TarImg_gridspace ...
, Boundary, Tarboundary, Mean, Scale, ref_select_size, X, Method,Tar_row, Tar_col,thresh_h,H_accuracy);

if 0
df=Zc.z-out{1,1}.z;
figure;imagesc(df);colorbar;caxis([-10 10])
figure;plot(df(:))
end

if 0 %for checking
tz=reshape(ref_select_nor_ori{3},ref_select_size);
ztar=reshape(tar_select_nor_h_ori,ref_select_size);
figure;imagesc(-tz+ztar);colorbar;caxis([0 0.001]);figure;imagesc(ztar);colorbar
figure;imagesc(reshape(ref_select_H,ref_select_size)-reshape(tar_H_array,([Tar_row,Tar_col]))); colorbar
[Boundary, Tarboundary, Mean, Scale]
load('tzs');load('ztars');
df=tz-tzs;figure;imagesc(df);colorbar;caxis([0 0.001]);
df=ztar-ztars;figure;imagesc(df);colorbar;caxis([0 0.001]);
load('tar_selects');
tz=reshape(tar_selects(:,3),ref_select_size);
tmp=reshape(tar_select{3},ref_select_size);
df=tmp-tz;figure;imagesc(df);colorbar;
load('tar_H_arrays'); %0
df=tar_H_array-tar_H_arrays;figure;plot(df)
end

Zref.x          = Zref.x - gridspace;
Zref.y          = Zref.y + gridspace;
Zc.x            = Zref.x;
Zc.y            = Zref.y;

outx{2} = out{2};%control points
outx{4} = out{4};%statistics
outx{1} = Zc;
outx{3} = Zref;
end
