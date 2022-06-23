
% �������ڣ�2021/05/03
% �������ܣ�% �Լ���ѹ��ϵ�����߽ڵ��Ž���������ʵʱ�����ǰͬ���ִ�ѹ���ֲ�����
            % �������ɹ��̡��ڵ���������ȫ����Ա��ɼ���󱣴����߰�����
% �����������˵��
% ���������shape_upper ���� ����������ԭʼ������ɢ����Ϣ
%           shape_lower ���� ����������ԭʼ������ɢ����Ϣ
%           com_data ����  ��������·��
% ���������shape_upper_arr_sor ���� ���������ļ������ݣ������棩
%           shape_lower_arr_sor ���� ���������ļ������ݣ������棩

function [shape_num,shape_arr_sor,n_com]=DA_Va_NonAirfoil_Ref_Adjustment_and_Visulization(com_data,fid_com_adj_path,fid_tar_new_path)

% % ����ȫ�ֱ���
% global savename new_output_folder fid_all fid_info PPNS initPath CFL_Setting Resi_target Resi_target_lev JX_NS STD_Comp JX_Turb L2NORM_Setting flag_Va_init;
% global S_Va_turb_Name S_Va_turb_Type S_Va_turb_Num S_Va_turb_FilePath S_Va_turb_ModiLine S_Va_turb_ModiFormat S_Va_turb_Std;
% global S_Va_target_Name S_Va_target_Num S_Va_target_TypeTurb S_Va_target_TypeDA S_Va_turb_target S_Va_Flow_length PT PS;
% global S_Va_ref_Name S_Va_ref_Num;
% global shape_upper shape_lower shape_outline shape_combined_arr_init Va_Beta Va_NNP_Beta Va_NNP_Beta_new;
% global flag_da_relative flag_sigma_obs_adapt flag_sigma_obs_coeff D_obs_e_ori D_obs_e;
% global Resi_flag Resi_flag_label;

% �����������
com_x=com_data(:,1);
com_y=com_data(:,2);
com_Cp=com_data(:,3);
com_Cf=com_data(:,4);
com_y_plus=com_data(:,5);
% ��¼����ڵ����
n_com=length(com_data);
shape_num=n_com;
shape_arr=zeros(shape_num,7);
% arr�����е�һ�д洢ԭʼ�������ݵ�ı��
%          �ڶ��д洢�������ݵ��x����
%          �����д洢�������ݵ��y����
%          �����д洢�������ݵ��Ӧ��ʵ��ֵ����ֵ������
for i=1:n_com
    shape_arr(i,1)=i;          % ����ԭʼ�������ݵ���
    shape_arr(i,2)=com_x(i);   % ����������ݵ��x����
    shape_arr(i,3)=com_y(i);   % ����������ݵ��y����
    shape_arr(i,4)=com_Cp(i);   % ����������ݵ��y����
    shape_arr(i,5)=com_Cf(i);   % ����������ݵ��y����
    shape_arr(i,6)=com_y_plus(i);   % ����������ݵ��y����
end
% �����������������x���������������ָ��[0,1]���䣩
shape_arr_sor=sortrows(shape_arr,2,'ascend');

% �ϲ��������������ʵ��ֵ����ֵ��������
shape_combined_arr_vis=shape_arr_sor;
%     % ���������ݵ�ԭʼ��Ž�����������
%     shape_combined_arr_sor=sortrows(shape_combined_arr_vis,1);

% ���������ļ������������ļ�
fid_com_adj_da=fopen(fid_com_adj_path,'w+');
% fprintf(fid_com_adj_da,'VARIABLES="X","Y","CP","Cf","yplus"\n');
% fprintf(fid_com_adj_da,'ZONE T="ZONE CAL_ADJ for DA"\n');
for i=1:length(shape_combined_arr_vis)
    fprintf(fid_com_adj_da,'%20.8f %20.8f %20.8f %20.8f %20.8f\n',shape_combined_arr_vis(i,2),shape_combined_arr_vis(i,3),shape_combined_arr_vis(i,4),shape_combined_arr_vis(i,5),shape_combined_arr_vis(i,6));
end
fclose(fid_com_adj_da);

% ���������Ŀ��ͬ�������ֲ���Cp��
fid_tar_new_da=fopen(fid_tar_new_path,'w+');
fprintf(fid_tar_new_da,'VARIABLES="X","CP"\n');
fprintf(fid_tar_new_da,'ZONE T="ZONE CP_NEW"\n');
for i=1:length(shape_combined_arr_vis)
    fprintf(fid_tar_new_da,'%20.8f %20.8f\n',shape_combined_arr_vis(i,2),shape_combined_arr_vis(i,4));
end
fclose(fid_tar_new_da);
    
end
