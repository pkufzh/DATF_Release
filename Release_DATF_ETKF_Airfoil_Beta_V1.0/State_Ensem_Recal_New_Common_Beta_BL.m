
% 更新日期：2021/05/20
% 函数功能：实现数据同化流场重计算
% 输入输出变量说明
% 输入变量：p —— 当前内迭代轮次
%           S_Va_ana_ave —— 当前预测得到的状态变量矩阵
%           MDA —— 同化工作变量维度
%           MStepNS_recal ——  设置重计算计算步数
%           targetPath —— 内迭代结果保存目标路径
% 输出变量：无

function State_Ensem_Recal_New_Common_Beta_BL(p,S_Va_ana_ave,MStepNS_recal,targetPath)

% 定义全局变量
global savename new_output_folder fid_all fid_info PPNS initPath CFL_Setting Resi_target Resi_target_lev JX_NS STD_Comp JX_Turb L2NORM_Setting flag_Va_init;
global S_Va_turb_Name S_Va_turb_Type S_Va_turb_Num S_Va_turb_FilePath S_Va_turb_ModiLine S_Va_turb_ModiFormat S_Va_turb_Std;
global S_Va_target_Name S_Va_target_Num S_Va_target_TypeTurb S_Va_target_TypeDA S_Va_turb_target S_Va_Flow_length PT PS;
global S_Va_ref_Name S_Va_ref_Num;
global shape_upper shape_lower shape_outline shape_combined_arr_init Va_Beta Va_NNP_Beta Va_NNP_Beta_new;
global flag_da_relative flag_sigma_obs_adapt flag_sigma_obs_coeff D_obs_e_ori D_obs_e flag_flow_type;
global BL_setting N_BL L_BL BL_num BL_lab BL_dot BL_obs_num kb;
global Resi_flag Resi_flag_label;

% 修改NS方程最大迭代步数：MStepNS -- NS方程最大迭代步数
% 指定替换行，对应输入文件（InputNS.dat）第9行
rep_line_MStepNS=9;
rep_format_MStepNS='%-16d';

% 修改续算开关： JX_NS -- 是否进行流场变量续算（0, 1）
%                JX_Turb -- 是否进行流场变量续算（0, 1, 3）
% 指定替换行，分别对应输入文件（InputNS.dat）第10行，（InputTurb.dat）第7行
rep_line_JX_NS=10;
rep_format_JX_NS='%-12d';
rep_line_JX_Turb=7;
rep_format_JX_Turb='%-12d';

% 修改计算CFL数：CFL_Setting
% 指定替换行，对应输入文件（InputNS.dat）第8行
rep_line_CFL=8;
rep_format_CFL='%-12d';
CFL_Recal=CFL_Setting(3);

% 修改标准算例计算结果比较开关：STD_Comp
% 指定替换行，对应输入文件（InputNS.dat）第22行
rep_line_STD_Comp=22;
rep_format_STD_Comp='%-12d';
STD_Comp_Recal=STD_Comp(3);

% 修改期望残差： Resi_target
rep_line_Resi=15;
rep_format_Resi='%-10.2e';
Resi_Recal=Resi_target(3);

% 修改期望 L2NORM 水平： L2NORM_Setting
rep_line_L2NORM=8;
rep_format_L2NORM='%-12d';
L2NORM_Recal=L2NORM_Setting(3);

target_runPath='.\Program_CFD_Parallel_Recal';

% 1. 更新（关闭）续算操纵开关
% JX_NS=0; % !!!!!!!!!!! 考虑可用前一轮重计算结果续算
% JX_Turb=1;
fid2=fopen([target_runPath,'\Input\InputNS.dat'],'r+');
for i=1:(rep_line_JX_NS-1)
   fgetl(fid2);
end
fseek(fid2,0,'cof');
fprintf(fid2,rep_format_JX_NS,JX_NS(3));
fclose(fid2);
fid2=fopen([target_runPath,'\Input\InputTurb.dat'],'r+');
for i=1:(rep_line_JX_Turb-1)
   fgetl(fid2);
end
fseek(fid2,0,'cof');
fprintf(fid2,rep_format_JX_Turb,JX_Turb(3));
fclose(fid2);

% 2. 设置NS方程最大运行迭代步数
fid3=fopen([target_runPath,'\Input\InputNS.dat'],'r+');
for i=1:(rep_line_MStepNS-1)
   fgetl(fid3);
end
fseek(fid3,0,'cof');
fprintf(fid3,rep_format_MStepNS,MStepNS_recal);
fclose(fid3);

% 3. 设置更新同化变量
% 更新目标同化变量(S_Va_ana)；若目标同化变量非初始扰动变量，则需要在更新步将值改为原始 STD 值
S_Va_cnt_One=0;
for ns=1:S_Va_target_Num
    % 提取目标同化变量类型（TypeTurb）：k - 是第k个初始扰动变量（k~=0）；0 - 非初始扰动变量
    % 若目标同化变量为初始扰动变量，且为单同化变量（e.g. "Ma" "AoA"）
    if ((S_Va_target_TypeTurb(ns)>0) && (S_Va_target_TypeDA(ns)==1))

        k=S_Va_target_TypeTurb(ns);
        S_Va_cnt_One=S_Va_cnt_One+1;

        fid3=fopen([target_runPath,char(S_Va_turb_FilePath(k))],'r+');
        for i=1:(S_Va_turb_ModiLine(k)-1)
           fgetl(fid3);
        end
        fseek(fid3,0,'cof');
        fprintf(fid3,char(S_Va_turb_ModiFormat(k)),S_Va_ana_ave(S_Va_cnt_One));                
        fclose(fid3);

    end
end

% 4. 设置计算CFL数
fid4=fopen([target_runPath,'\Input\InputNS.dat'],'r+');
for i=1:(rep_line_CFL-1)
   fgetl(fid4);
end
fseek(fid4,0,'cof');
fprintf(fid4,rep_format_CFL,CFL_Recal);
fclose(fid4);

% 5. 设置与标准算例比较开关
fid5=fopen([target_runPath,'\Input\InputNS.dat'],'r+');
for i=1:(rep_line_STD_Comp-1)
   fgetl(fid5);
end
fseek(fid5,0,'cof');
fprintf(fid5,rep_format_STD_Comp,STD_Comp_Recal);
fclose(fid5);
            
% 6. 设置计算收敛精度
fid6=fopen([target_runPath,'\Input\InputNS.dat'],'r+');
for i=1:(rep_line_Resi-1)
   fgetl(fid6);
end
fseek(fid6,0,'cof');
fprintf(fid6,rep_format_Resi,Resi_Recal);
fclose(fid6);

% 7. 设置L2NORM收敛水平
fid7=fopen([target_runPath,'\Input\InputTurb.dat'],'r+');
for i=1:(rep_line_L2NORM-1)
   fgetl(fid7);
end
fseek(fid7,0,'cof');
fprintf(fid7,rep_format_L2NORM,L2NORM_Recal);
fclose(fid7);
            
% 8. 恢复指定扰动变量为 STD 值
if (flag_Va_init==0)
    for ns=1:S_Va_turb_Num
        % 若初始扰动变量非目标同化变量，且为单变量（e.g. "CKS"）
        % 则将该变量恢复原始 STD 值
        if (S_Va_turb_target(ns)==0)

            fid8=fopen([target_runPath,char(S_Va_turb_FilePath(ns))],'r+');
            for i=1:(S_Va_turb_ModiLine(ns)-1)
               fgetl(fid8);
            end
            fseek(fid8,0,'cof');
            fprintf(fid8,char(S_Va_turb_ModiFormat(ns)),S_Va_turb_Std(ns));
            fclose(fid8);

        end
    end
end

% 9. 更新流场工作变量文件

% % 更新流程工作变量 TMuData 文件，并覆盖 Output 文件夹下的续算文件
% % 更新当前流场数据（湍流工作变量）
% TMu_new=S_Va_ana_ave((PT+1):MDA);
% % 更新工作变量 TMuData 文件
% save('.\Program_CFD_Parallel_Recal\Output\TMuData.PLT','TMu_new','-ascii');

S_Va_cnt_Two=0;
SS=PT+1; % 数组索引起点
ST=PT; % 数组索引终点（初值，待后续迭代更新）
for i=1:S_Va_target_Num
    if (S_Va_target_TypeDA(i)==2)

        S_Va_cnt_Two=S_Va_cnt_Two+1;

        % 计算数组索引终点，迭代更新
        SL=S_Va_Flow_length(S_Va_cnt_Two);
        ST=ST+SL;

        if (strcmp(char(S_Va_target_Name(i)),'Vc')==1)
            Vc_new=S_Va_ana_ave(SS:ST);
            save('.\Program_CFD_Parallel_Recal\Output\VcData.PLT','Vc_new','-ascii');
            SS=ST+1;
        elseif (strcmp(char(S_Va_target_Name(i)),'TMu')==1)
            TMu_new=S_Va_ana_ave(SS:ST);
            save('.\Program_CFD_Parallel_Recal\Output\TMuData.PLT','TMu_new','-ascii');
            SS=ST+1;
        elseif (strcmp(char(S_Va_target_Name(i)),'TMuP')==1)
            TMuP_new=S_Va_ana_ave(SS:ST);
            save('.\Program_CFD_Parallel_Recal\Output\TMuPData.PLT','TMuP_new','-ascii');
            SS=ST+1;
        elseif (strcmp(char(S_Va_target_Name(i)),'Beta')==1)
            Beta_new=S_Va_ana_ave(SS:ST);
            save('.\Program_CFD_Parallel_Recal\Output\Beta_plus.PLT','Beta_new','-ascii'); % 将同化得到的单元修正项值信息写入指定文件
            save('.\Program_CFD_Parallel_Recal\Output\Beta_plus_NNP.PLT','Va_NNP_Beta_new','-ascii'); % 将同化并插值得到的节点修正项值信息写入指定文件
            SS=ST+1;
        end

    end
end

% 删除已有程序完成标记文件
if (exist([target_runPath,'\ProgramExitFlag.txt'],'file'))~=0
    delete([target_runPath,'\ProgramExitFlag.txt']);
end

% 10. 运行 CFD 程序，生成初始集合成员样本并保存
% 运行CFD求解主程序
fprintf(fid_all,'**************************************************\n');
fprintf(fid_all,'|                 PROGRAM STATUS                  \n');
fprintf(fid_all,'**************************************************\n');
fprintf(fid_all,'|             RE-CALCULATION PROCEDURE            \n');
fprintf(fid_all,'**************************************************\n');
fprintf(fid_all,'|                                                 \n');
fprintf(fid_all,'|  The calculation of DA_COM sample is RUNNING!   \n');
fprintf(fid_all,'|                                                 \n');
fprintf(fid_all,'**************************************************\n');
fprintf(fid_all,'\n');

fprintf(['集合成员重计算 ',num2str(p,'%03d'),' 轮：状态']);
fprintf(fid_info,['集合成员重计算 ',num2str(p,'%03d'),' 轮：状态']);
% 输出时仅显示 目标同化变量中的单变量
S_Va_cnt_One=0;
for ns=1:S_Va_target_Num
    % 若目标同化变量为初始扰动变量，且为单同化变量（e.g. "Ma" "AoA"）
    if ((S_Va_target_TypeTurb(ns)>0) && (S_Va_target_TypeDA(ns)==1))
        S_Va_cnt_One=S_Va_cnt_One+1;
        fprintf([' ',char(S_Va_target_Name(ns)),' = %.5f'],S_Va_ana_ave(S_Va_cnt_One));
        fprintf(fid_info,[' ',char(S_Va_target_Name(ns)),' = %.5f'],S_Va_ana_ave(S_Va_cnt_One));
    end
end
fprintf(' 已启动计算\n');
fprintf(fid_info,' 已启动计算\n');

% cmd=('.\Program_CFD_Parallel_Recal\CFD_2D.exe');
% open(cmd);
[s,e]=dos('cd Program_CFD_Parallel_Recal && run_CFD_2D.bat&');
pause(5); % 暂停数秒

% 设置重迭代过程退出标志
exit_recal_flag=0;

while (exit_recal_flag~=1)

    if (exist([target_runPath,'\ProgramExitFlag.txt'],'file'))~=0
        
        % -------------------- 进行结果保存 START ------------------- %

        % 提取当前状态计算最终残差值 Resi(dual)_Final
        if (JX_Turb(3)==0)
            Resi_All=load([target_runPath,'\Output\Residual_1.PLT']);
        elseif ((JX_Turb(3)==1) || (JX_Turb(3)==2))
            Resi_All=load([target_runPath,'\Output\Residual_1_Turb_Const.PLT']);
        elseif (JX_Turb(3)==3)
            Resi_All=load([target_runPath,'\Output\Residual_1_Turb_Beta.PLT']);
        end
        % Resi_All=load([target_runPath,'\Output\Residual_1_Turb_Const.PLT']);
        Resi_Final=Resi_All(size(Resi_All,1),2);
        % 导入计算数据
        com_data=load('.\Program_CFD_Parallel_Recal\Output\CP_CF_Comp.PLT');
        % 指定调整后计算数据保存路径
        fid_com_adj_path='.\Program_CFD_Parallel_Recal\Output\CP_CF_Comp_New.PLT';
        fid_tar_new_path='.\Program_CFD_Parallel_Recal\Output\CP_New.PLT';

        % 对计算压力系数曲线节点编号进行梳理，并实时输出当前同化轮次压力分布曲线
        % 集合生成过程、内迭代过程中全部成员完成计算后保存曲线包络线
        if (flag_flow_type==0)
            [shape_upper_num,shape_lower_num,shape_upper_arr_sor,shape_lower_arr_sor,n_com]=DA_Va_Airfoil_Ref_Adjustment_and_Visulization(com_data,fid_com_adj_path,fid_tar_new_path);
            % 合并上下翼面计算与实验值（插值）点数据
            shape_combined_arr_recal=[shape_lower_arr_sor;shape_upper_arr_sor]; % 先下翼面，后上翼面
        elseif (flag_flow_type==1)
            [shape_num,shape_arr_sor,n_com]=DA_Va_NonAirfoil_Ref_Adjustment_and_Visulization(com_data,fid_com_adj_path,fid_tar_new_path);
            shape_combined_arr_recal=shape_arr_sor;
        end
        
        % 初始标准状态计算结果可视化
        figure;
        hold on
        plot(shape_combined_arr_init(:,2),shape_combined_arr_init(:,7),'ro','MarkerSize',4.0);
        plot(shape_combined_arr_recal(:,2),shape_combined_arr_recal(:,4),'b-','LineWidth',2.0);
        xlabel('\itx'),ylabel('\itCp'),title(['重计算 ',num2str(p,'%03d'),' 轮计算结果']);
        legend('实验测量值','重计算翼型压力系数分布');
        grid on
        box on
        hold off
        % 保存可编辑图片
        saveas(gcf,['.\Output_Figures\',savename,'\Recal_Round_',num2str(p,'%03d'),'_Output.fig']);
        
        % 输出命令行提示
        fprintf(fid_all,'**************************************************\n');
        fprintf(fid_all,'|                 PROGRAM STATUS                  \n');
        fprintf(fid_all,'**************************************************\n');
        fprintf(fid_all,'|             RE-CALCULATION PROCEDURE            \n');
        fprintf(fid_all,'**************************************************\n');
        fprintf(fid_all,'|                                                 \n');
        fprintf(fid_all,'|  The calculation of DA_COM sample is FINISHED!  \n');
        fprintf(fid_all,'|                                                 \n');
        fprintf(fid_all,'**************************************************\n');
        fprintf(fid_all,'\n');
        fprintf(fid_all,'\n');
        
        fprintf(['集合成员重计算 ',num2str(p,'%03d'),' 轮：状态已保存完毕, Total Step = %d，Final Residual = %.5f\n'],(size(Resi_All,1)-1),Resi_Final);
        fprintf(fid_info,['集合成员重计算 ',num2str(p,'%03d'),' 轮：状态已保存完毕, Total Step = %d，Final Residual = %.5f\n'],(size(Resi_All,1)-1),Resi_Final);

        if (strcmp(char(S_Va_ref_Name(1)),'BL')==1)
            % 运行提取边界层信息批处理文件
            fprintf(['集合成员重计算 ',num2str(p,'%03d'),' 边界层数据提取开始\n']);
            fprintf(fid_info,['集合成员重计算 ',num2str(p,'%03d'),' 边界层数据提取开始\n']);
            % [s,e]=dos(['cd ',target_runPath,' && run_BL_Hump_Extract_Tecplot.bat']);
            [s,e]=dos('cd Program_CFD_Parallel_Recal && run_BL_Hump_Extract_Tecplot.bat');
            fprintf(['集合成员重计算 ',num2str(p,'%03d'),' 边界层数据提取完毕\n']);
            fprintf(fid_info,['集合成员重计算 ',num2str(p,'%03d'),' 边界层数据提取完毕\n']);
            pause(3); % 暂停数秒
        end
            
        % 输出文件另存，覆盖预测步之前的计算结果
        save_recal='.\Program_CFD_Parallel_Recal\Output';
        fileList=dir(save_recal);
        mkdir([targetPath,'\Recalculation_Round_',num2str(p,'%03d')]); % 利用 mkdir() 函数创建新文件夹
        filename_s=cell(length(fileList),1);
        filename_t=cell(length(fileList),1);

        for i=3:length(fileList)
              filename_s{i}=[save_recal,'\',fileList(i).name];
              filename_t{i}=[targetPath,'\Recalculation_Round_',num2str(p,'%03d')];
              copyfile(filename_s{i},filename_t{i});
        end

        fprintf('\n');
        fprintf(fid_info,'\n');

        % -------------------- 进行结果保存 END -------------------- %
        
        exit_recal_flag=1;  % 退出标志置一
        
    end
    
    pause(10);
    
end

end
