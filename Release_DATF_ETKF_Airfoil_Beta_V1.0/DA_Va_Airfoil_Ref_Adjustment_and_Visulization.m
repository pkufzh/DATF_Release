
% �������ڣ�2021/03/18
% �������ܣ�% �Լ���ѹ��ϵ�����߽ڵ��Ž���������ʵʱ�����ǰͬ���ִ�ѹ���ֲ�����
            % �������ɹ��̡��ڵ���������ȫ����Ա��ɼ���󱣴����߰�����
% �����������˵��
% ���������shape_upper ���� ����������ԭʼ������ɢ����Ϣ
%           shape_lower ���� ����������ԭʼ������ɢ����Ϣ
%           com_data ����  ��������·��
% ���������shape_upper_arr_sor ���� ���������ļ������ݣ������棩
%           shape_lower_arr_sor ���� ���������ļ������ݣ������棩

function [shape_upper_num,shape_lower_num,shape_upper_arr_sor,shape_lower_arr_sor,n_com]=DA_Va_Airfoil_Ref_Adjustment_and_Visulization(com_data,fid_com_adj_path,fid_tar_new_path)

global shape_upper shape_lower shape_combined_arr_init;

% �������ݲ���
% ��������ԭʼ�������ݣ��������������Ͳ�ֵ
xx_upper_com=(0:0.001:1);
yy_upper_com=interp1(shape_upper(:,1),shape_upper(:,2),xx_upper_com,'linear');
xx_lower_com=(0:0.001:1);
yy_lower_com=interp1(shape_lower(:,1),shape_lower(:,2),xx_lower_com,'linear');

% �����������
com_x=com_data(:,1);
com_y=com_data(:,2);
com_Cp=com_data(:,3);
com_Cf=com_data(:,4);
com_y_plus=com_data(:,5);
% ��¼����ڵ����
n_com=length(com_data);
% �жϼ��������б߽����������λ�ã��ϻ��£�
% 0 ���������棻1 ����������
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
% �ֿ����������������ݵ㼰�������
% ������������������ݵ����
shape_lower_num=sum(flag_shape_pos);
shape_upper_num=n_com-shape_lower_num;
% �����洢��������������ݵ������
shape_upper_arr=zeros(shape_upper_num,7);
shape_lower_arr=zeros(shape_lower_num,7);
% arr�����е�һ�д洢ԭʼ�������ݵ�ı��
%          �ڶ��д洢�������ݵ��x����
%          �����д洢�������ݵ��y����
%          �����д洢�������ݵ��Ӧ��ʵ��ֵ����ֵ������
cnt_upper=0;
cnt_lower=0;
for i=1:n_com
    % ���õ�Ϊ���������ݵ�
    if (flag_shape_pos(i)==0)
        cnt_upper=cnt_upper+1;
        shape_upper_arr(cnt_upper,1)=i;          % ����ԭʼ�������ݵ���
        shape_upper_arr(cnt_upper,2)=com_x(i);   % ����������ݵ��x����
        shape_upper_arr(cnt_upper,3)=com_y(i);   % ����������ݵ��y����
        shape_upper_arr(cnt_upper,4)=com_Cp(i);   % ����������ݵ��y����
        shape_upper_arr(cnt_upper,5)=com_Cf(i);   % ����������ݵ��y����
        shape_upper_arr(cnt_upper,6)=com_y_plus(i);   % ����������ݵ��y����
    else
        cnt_lower=cnt_lower+1;
        shape_lower_arr(cnt_lower,1)=i;          % ����ԭʼ�������ݵ���
        shape_lower_arr(cnt_lower,2)=com_x(i);   % ����������ݵ��x����
        shape_lower_arr(cnt_lower,3)=com_y(i);   % ����������ݵ��y����
        shape_lower_arr(cnt_lower,4)=com_Cp(i);   % ����������ݵ��y����
        shape_lower_arr(cnt_lower,5)=com_Cf(i);   % ����������ݵ��y����
        shape_lower_arr(cnt_lower,6)=com_y_plus(i);   % ����������ݵ��y����
    end
end
% ���������������������x�����������ָ��[0,1]���䣩�������͵�������������
shape_upper_arr_sor=sortrows(shape_upper_arr,2,'descend');
shape_lower_arr_sor=sortrows(shape_lower_arr,2,'ascend');

% �ϲ��������������ʵ��ֵ����ֵ��������
shape_combined_arr_vis=[shape_lower_arr_sor;shape_upper_arr_sor];
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
