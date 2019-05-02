function [rho] = Correlation(varargin)

value1 = varargin{1};
value2 = varargin{2};
window = varargin{3};

t_size = size(value1{1});
SUM_LR = zeros(t_size(1),1);
SUM_L = zeros(t_size(1),1);
SUM_R = zeros(t_size(1),1);
SUM_L2 = zeros(t_size(1),1);
SUM_R2 = zeros(t_size(1),1);
rho = zeros(t_size(1),1);
for i=1:window*window
    SUM_LR = SUM_LR + value1{i}.*value2{i};
    SUM_L  = SUM_L + value1{i};
    SUM_R  = SUM_R + value2{i};
    SUM_L2 = SUM_L2 + value1{i}.*value1{i};
    SUM_R2 = SUM_R2 + value2{i}.*value2{i};
end
val1 = SUM_L2 - (SUM_L .*SUM_L)./(window*window);
val2 = SUM_R2 - (SUM_R .*SUM_R)./(window*window);


de1 = sqrt(val1.*val2);
t_index_1 = find(de1 == 0);
t_size_1 = size(t_index_1);
t_index_2 = find(de1 ~= 0);
t_size_2 = size(t_index_2);

temp = (SUM_L(t_index_2).*SUM_R(t_index_2))./(window*window);
de2 = (SUM_LR(t_index_2) - temp);
rho(t_index_2,1) = (de2./de1(t_index_2) + ones(t_size_2(1),1))./2;
rho(t_index_1,1) = zeros(t_size_1(1),1);

