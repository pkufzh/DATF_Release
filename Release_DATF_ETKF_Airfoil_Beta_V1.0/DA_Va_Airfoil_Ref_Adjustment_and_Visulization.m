
% 更新日期：2021/03/18
% 函数功能：% 对计算压力系数曲线节点编号进行梳理，并实时输出当前同化轮次压力分布曲线
            % 集合生成过程、内迭代过程中全部成员完成计算后保存曲线包络线
% 输入输出变量说明
% 输入变量：shape_upper —— 翼型上翼面原始外形离散点信息
%           shape_lower —— 翼型下翼面原始外形离散点信息
%           com_data ——  计算数据路径
% 输出变量：shape_upper_arr_sor —— 经过排序后的计算数据（上翼面）
%           shape_lower_arr_sor —— 经过排序后的计算数据（下翼面）

function [shape_upper_num,shape_lower_num,shape_upper_arr_sor,shape_lower_arr_sor,n_com]=DA_Va_Airfoil_Ref_Adjustment_and_Visulization(com_data,fid_com_adj_path,fid_tar_new_path)

global shape_upper shape_lower shape_combined_arr_init;

% 计算数据操作
% 导入翼型原始外形数据，并进行上下翼型插值
xx_upper_com=(0:0.001:1);
yy_upper_com=interp1(shape_upper(:,1),shape_upper(:,2),xx_upper_com,'linear');
xx_lower_com=(0:0.001:1);
yy_lower_com=interp1(shape_lower(:,1),shape_lower(:,2),xx_lower_com,'linear');

% 导入计算数据
com_x=com_data(:,1);
com_y=com_data(:,2);
com_Cp=com_data(:,3);
com_Cf=com_data(:,4);
com_y_plus=com_data(:,5);
% 记录计算节点个数
n_com=length(com_data);
% 判断计算数据中边界点所在翼型位置（上或下）
% 0 代表上翼面；1 代表下翼面
flag_shape_pos=zeros(1,n_com);
for i=1:n_com
    y_itep_upper=interp1(shape_upper(:,1),shape_upper(:,2),com_x(i),'linear');
    y_itep_lower=interp1(shape_lower(:,1),shape_lower(:,2),com_x(i),'linear');
    y_itep_upper_err=abs(y_itep_upper-com_y(i));
    y_itep_lower_err=abs(y_itep_lower-com_y(i));
    if (y_itep_lower_err<y_itep_upper_err)
        flag_shape_pos(i)=1;
    end
end
% 分开储存上下翼面数据点及其计算编号
% 计算上下翼面计算数据点个数
shape_lower_num=sum(flag_shape_pos);
shape_upper_num=n_com-shape_lower_num;
% 创建存储上翼翼面计算数据点的数组
shape_upper_arr=zeros(shape_upper_num,7);
shape_lower_arr=zeros(shape_lower_num,7);
% arr数组中第一列存储原始计算数据点的编号
%          第二列存储计算数据点的x坐标
%          第三列存储计算数据点的y坐标
%          第四列存储计算数据点对应的实验值（插值方法）
cnt_upper=0;
cnt_lower=0;
for i=1:n_com
    % 若该点为上翼型数据点
    if (flag_shape_pos(i)==0)
        cnt_upper=cnt_upper+1;
        shape_upper_arr(cnt_upper,1)=i;          % 保存原始计算数据点编号
        shape_upper_arr(cnt_upper,2)=com_x(i);   % 保存计算数据点的x坐标
        shape_upper_arr(cnt_upper,3)=com_y(i);   % 保存计算数据点的y坐标
        shape_upper_arr(cnt_upper,4)=com_Cp(i);   % 保存计算数据点的y坐标
        shape_upper_arr(cnt_upper,5)=com_Cf(i);   % 保存计算数据点的y坐标
        shape_upper_arr(cnt_upper,6)=com_y_plus(i);   % 保存计算数据点的y坐标
    else
        cnt_lower=cnt_lower+1;
        shape_lower_arr(cnt_lower,1)=i;          % 保存原始计算数据点编号
        shape_lower_arr(cnt_lower,2)=com_x(i);   % 保存计算数据点的x坐标
        shape_lower_arr(cnt_lower,3)=com_y(i);   % 保存计算数据点的y坐标
        shape_lower_arr(cnt_lower,4)=com_Cp(i);   % 保存计算数据点的y坐标
        shape_lower_arr(cnt_lower,5)=com_Cf(i);   % 保存计算数据点的y坐标
        shape_lower_arr(cnt_lower,6)=com_y_plus(i);   % 保存计算数据点的y坐标
    end
end
% 将上下翼型数据数组根据x坐标进行排序（指定[0,1]区间），上翼型倒序，下翼型升序
shape_upper_arr_sor=sortrows(shape_upper_arr,2,'descend');
shape_lower_arr_sor=sortrows(shape_lower_arr,2,'ascend');

% 合并上下翼面计算与实验值（插值）点数据
shape_combined_arr_vis=[shape_lower_arr_sor;shape_upper_arr_sor];
%     % 按计算数据点原始编号进行升序排序
%     shape_combined_arr_sor=sortrows(shape_combined_arr_vis,1);

% 保存整理后的计算数据至新文件
fid_com_adj_da=fopen(fid_com_adj_path,'w+');
% fprintf(fid_com_adj_da,'VARIABLES="X","Y","CP","Cf","yplus"\n');
% fprintf(fid_com_adj_da,'ZONE T="ZONE CAL_ADJ for DA"\n');
for i=1:length(shape_combined_arr_vis)
    fprintf(fid_com_adj_da,'%20.8f %20.8f %20.8f %20.8f %20.8f\n',shape_combined_arr_vis(i,2),shape_combined_arr_vis(i,3),shape_combined_arr_vis(i,4),shape_combined_arr_vis(i,5),shape_combined_arr_vis(i,6));
end
fclose(fid_com_adj_da);

% 输出排序后的目标同化变量分布（Cp）
fid_tar_new_da=fopen(fid_tar_new_path,'w+');
fprintf(fid_tar_new_da,'VARIABLES="X","CP"\n');
fprintf(fid_tar_new_da,'ZONE T="ZONE CP_NEW"\n');
for i=1:length(shape_combined_arr_vis)
    fprintf(fid_tar_new_da,'%20.8f %20.8f\n',shape_combined_arr_vis(i,2),shape_combined_arr_vis(i,4));
end
fclose(fid_tar_new_da);
    
end
