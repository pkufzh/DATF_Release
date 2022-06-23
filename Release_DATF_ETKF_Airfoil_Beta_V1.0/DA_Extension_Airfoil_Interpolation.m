
%% 翼型外形（输出）插值函数
function DA_Extension_Airfoil_Interpolation(NNDA,XDA,YDA)
    
    %% 定义全局变量
    global savename new_output_folder fid_all fid_info PPNS initPath CFL_Setting Resi_target Resi_target_lev JX_NS STD_Comp JX_Turb L2NORM_Setting flag_Va_init;
    global S_Va_turb_Name S_Va_turb_Type S_Va_turb_Num S_Va_turb_FilePath S_Va_turb_ModiLine S_Va_turb_ModiFormat S_Va_turb_Std;
    global S_Va_target_Name S_Va_target_Num S_Va_target_TypeTurb S_Va_target_TypeDA S_Va_turb_target S_Va_Flow_length PT PS;
    global S_Va_ref_Name S_Va_ref_Num;
    global shape_upper shape_lower shape_outline shape_combined_arr_init Va_Beta Va_NNP_Beta Va_NNP_Beta_new;
    global flag_da_relative flag_sigma_obs_adapt flag_sigma_obs_coeff D_obs_e_ori D_obs_e;
    global Resi_flag Resi_flag_label;

    %% 导入翼型原始插值点数据
    % 导入翼型原始外形数据
    %     shape_upper=load('.\Input_Data\Input_Data_Shape_S809_upper.dat');
    %     shape_lower=load('.\Input_Data\Input_Data_Shape_S809_lower.dat');

    %% 分别进行上下翼型插值（范围为[0,1]）
    % 设置插值步长
    grid_gap=0.001;
    grid_gap_front=0.0001;

    % 设置前缘分割点
    xx_front=0.05;
    xx_upper_min=min(shape_upper(:,1));
    xx_upper_max=max(shape_upper(:,1));
    xx_upper_com=[(xx_upper_min:grid_gap_front:xx_front),((xx_front+grid_gap):grid_gap:xx_upper_max)]';
    yy_upper_com=interp1(shape_upper(:,1),shape_upper(:,2),xx_upper_com,'spline');
    shape_upper_com=[xx_upper_com,yy_upper_com];
    xx_lower_min=min(shape_lower(:,1));
    xx_lower_max=max(shape_lower(:,1));
    xx_lower_com=[(xx_lower_min:grid_gap_front:xx_front),((xx_front+grid_gap):grid_gap:xx_lower_max)]';
    yy_lower_com=interp1(shape_lower(:,1),shape_lower(:,2),xx_lower_com,'spline');
    shape_lower_com=[xx_lower_com,yy_lower_com];

    % 合并上下翼型外形数据
    shape_upper_com_sort=sortrows(shape_upper_com,1,'ascend');
    shape_lower_com_sort=sortrows(shape_lower_com,1,'descend');
    shape_com_all=[shape_upper_com_sort;shape_lower_com_sort];

    % 计算插值后翼型点总数
    n_com_all=length(shape_com_all);

    %% 翼型插值结果与指定同化区域可视化
    figure;
    hold on
    plot(xx_upper_com,yy_upper_com,'b-','LineWidth',2);
    plot(xx_lower_com,yy_lower_com,'b-','LineWidth',2);
    for i=1:(NNDA-1)
        PDA_pre=[XDA(i),XDA(i+1)];
        PDA_aft=[YDA(i),YDA(i+1)];
        hDA(i)=plot(PDA_pre,PDA_aft,'r-','LineWidth',2);
        set(hDA(i),'handlevisibility','off');
    end
    plot([XDA(NNDA),XDA(1)],[YDA(NNDA),YDA(1)],'r-','LineWidth',2);
    %     plot(XDA,YDA,'r-','LineWidth',2);
    % 设置适合查看的窗口大小
    eDA=0.2; % 坐标缩放尺寸
    XDA_min=min(XDA)-eDA*abs(max(XDA));
    XDA_max=max(XDA)+eDA*abs(max(XDA));
    YDA_min=min(YDA)-eDA*abs(max(YDA));
    YDA_max=max(YDA)+eDA*abs(max(YDA));
    axis([XDA_min,XDA_max,YDA_min,YDA_max]);
    xlabel('X Coordinate');
    ylabel('Y Coordinate');
    legend('Upper Surface','Lower Surface','DA Region Boundary');
    
    set(gcf,'renderer','painters'); % 设置输出格式为矢量图
    % 保存可编辑图片
    saveas(gcf,['.\Output_Figures\',savename,'\Init_State_Output.fig']);
    
    grid on
    box on
    hold off

    %% 保存翼型插值结果至.dat文件
    % 调整适合 Pointwise 的输入数据格式
    shape_com_all_output=[shape_com_all,zeros(n_com_all,1)];

    fid_save=fopen('Input_Data_Shape_Airfoil_All.dat','w+');
    fprintf(fid_save,'%d\n',n_com_all);
    for i=1:n_com_all
        fprintf(fid_save,'%-20.10f %-20.10f %-20.10f\n',shape_com_all_output(i,1),shape_com_all_output(i,2),shape_com_all_output(i,3));
    end
    fclose(fid_save);

end
