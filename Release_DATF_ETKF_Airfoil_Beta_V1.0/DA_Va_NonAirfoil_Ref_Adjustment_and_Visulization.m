
% 更新日期：2021/05/03
% 函数功能：% 对计算压力系数曲线节点编号进行梳理，并实时输出当前同化轮次压力分布曲线
            % 集合生成过程、内迭代过程中全部成员完成计算后保存曲线包络线
% 输入输出变量说明
% 输入变量：shape_upper —— 翼型上翼面原始外形离散点信息
%           shape_lower —— 翼型下翼面原始外形离散点信息
%           com_data ——  计算数据路径
% 输出变量：shape_upper_arr_sor —— 经过排序后的计算数据（上翼面）
%           shape_lower_arr_sor —— 经过排序后的计算数据（下翼面）

function [shape_num,shape_arr_sor,n_com]=DA_Va_NonAirfoil_Ref_Adjustment_and_Visulization(com_data,fid_com_adj_path,fid_tar_new_path)

% % 定义全局变量
% global savename new_output_folder fid_all fid_info PPNS initPath CFL_Setting Resi_target Resi_target_lev JX_NS STD_Comp JX_Turb L2NORM_Setting flag_Va_init;
% global S_Va_turb_Name S_Va_turb_Type S_Va_turb_Num S_Va_turb_FilePath S_Va_turb_ModiLine S_Va_turb_ModiFormat S_Va_turb_Std;
% global S_Va_target_Name S_Va_target_Num S_Va_target_TypeTurb S_Va_target_TypeDA S_Va_turb_target S_Va_Flow_length PT PS;
% global S_Va_ref_Name S_Va_ref_Num;
% global shape_upper shape_lower shape_outline shape_combined_arr_init Va_Beta Va_NNP_Beta Va_NNP_Beta_new;
% global flag_da_relative flag_sigma_obs_adapt flag_sigma_obs_coeff D_obs_e_ori D_obs_e;
% global Resi_flag Resi_flag_label;

% 导入计算数据
com_x=com_data(:,1);
com_y=com_data(:,2);
com_Cp=com_data(:,3);
com_Cf=com_data(:,4);
com_y_plus=com_data(:,5);
% 记录计算节点个数
n_com=length(com_data);
shape_num=n_com;
shape_arr=zeros(shape_num,7);
% arr数组中第一列存储原始计算数据点的编号
%          第二列存储计算数据点的x坐标
%          第三列存储计算数据点的y坐标
%          第四列存储计算数据点对应的实验值（插值方法）
for i=1:n_com
    shape_arr(i,1)=i;          % 保存原始计算数据点编号
    shape_arr(i,2)=com_x(i);   % 保存计算数据点的x坐标
    shape_arr(i,3)=com_y(i);   % 保存计算数据点的y坐标
    shape_arr(i,4)=com_Cp(i);   % 保存计算数据点的y坐标
    shape_arr(i,5)=com_Cf(i);   % 保存计算数据点的y坐标
    shape_arr(i,6)=com_y_plus(i);   % 保存计算数据点的y坐标
end
% 将外形数据数组根据x坐标进行升序排序（指定[0,1]区间）
shape_arr_sor=sortrows(shape_arr,2,'ascend');

% 合并上下翼面计算与实验值（插值）点数据
shape_combined_arr_vis=shape_arr_sor;
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
