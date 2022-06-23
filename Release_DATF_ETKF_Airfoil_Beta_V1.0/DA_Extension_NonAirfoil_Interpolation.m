
%% �������������Σ��������ֵ����
function DA_Extension_NonAirfoil_Interpolation(NNDA,XDA,YDA)
    
    %% ����ȫ�ֱ���
    global savename new_output_folder fid_all fid_info PPNS initPath CFL_Setting Resi_target Resi_target_lev JX_NS STD_Comp JX_Turb L2NORM_Setting flag_Va_init;
    global S_Va_turb_Name S_Va_turb_Type S_Va_turb_Num S_Va_turb_FilePath S_Va_turb_ModiLine S_Va_turb_ModiFormat S_Va_turb_Std;
    global S_Va_target_Name S_Va_target_Num S_Va_target_TypeTurb S_Va_target_TypeDA S_Va_turb_target S_Va_Flow_length PT PS;
    global S_Va_ref_Name S_Va_ref_Num;
    global shape_upper shape_lower shape_outline shape_combined_arr_init Va_Beta Va_NNP_Beta Va_NNP_Beta_new;
    global flag_da_relative flag_sigma_obs_adapt flag_sigma_obs_coeff D_obs_e_ori D_obs_e;
    global Resi_flag Resi_flag_label;

    %% ��������ԭʼ��ֵ������
    % ��������ԭʼ��������
    %     shape_upper=load('.\Input_Data\Input_Data_Shape_S809_upper.dat');
    %     shape_lower=load('.\Input_Data\Input_Data_Shape_S809_lower.dat');

    %% �ֱ�����������Ͳ�ֵ����ΧΪ[0,1]��
    % ���ò�ֵ����
    grid_gap=0.001;

    % ���ò�ֵ����
    xx_min=min(shape_outline(:,1));
    xx_max=max(shape_outline(:,1));
    xx_com=(xx_min:grid_gap:xx_max)';
    yy_com=interp1(shape_outline(:,1),shape_outline(:,2),xx_com,'spline');
    shape_com=[xx_com,yy_com];

    % �ϲ�����������������
    shape_com_sort=sortrows(shape_com,1,'ascend');

    % �����ֵ������ѹ������
    n_com_all=length(shape_com_sort);

    %% ���Ͳ�ֵ�����ָ��ͬ��������ӻ�
    figure;
    hold on
    plot(xx_com,yy_com,'b-','LineWidth',2);
    for i=1:(NNDA-1)
        PDA_pre=[XDA(i),XDA(i+1)];
        PDA_aft=[YDA(i),YDA(i+1)];
        hDA(i)=plot(PDA_pre,PDA_aft,'r-','LineWidth',2);
        set(hDA(i),'handlevisibility','off');
    end
    plot([XDA(NNDA),XDA(1)],[YDA(NNDA),YDA(1)],'r-','LineWidth',2);
    %     plot(XDA,YDA,'r-','LineWidth',2);
    % �����ʺϲ鿴�Ĵ��ڴ�С
    eDA=0.2; % �������ųߴ�
    XDA_min=min(XDA)-eDA*abs(max(XDA));
    XDA_max=max(XDA)+eDA*abs(max(XDA));
    YDA_min=min(YDA)-eDA*abs(max(YDA));
    YDA_max=max(YDA)+eDA*abs(max(YDA));
    axis([XDA_min,XDA_max,YDA_min,YDA_max]);
    xlabel('X Coordinate');
    ylabel('Y Coordinate');
    legend('Upper Surface','Lower Surface','DA Region Boundary');
    
    set(gcf,'renderer','painters'); % ���������ʽΪʸ��ͼ
    % ����ɱ༭ͼƬ
    saveas(gcf,['.\Output_Figures\',savename,'\Init_State_Output.fig']);
    
    grid on
    box on
    hold off

    %% �������Ͳ�ֵ�����.dat�ļ�
    % �����ʺ� Pointwise ���������ݸ�ʽ��z����ֵ����Ϊ��
    shape_com_all_output=[shape_com_sort,zeros(n_com_all,1)];

    fid_save=fopen('Input_Data_Shape_Airfoil_All.dat','w+');
    fprintf(fid_save,'%d\n',n_com_all);
    for i=1:n_com_all
        fprintf(fid_save,'%-20.10f %-20.10f %-20.10f\n',shape_com_all_output(i,1),shape_com_all_output(i,2),shape_com_all_output(i,3));
    end
    fclose(fid_save);

end
