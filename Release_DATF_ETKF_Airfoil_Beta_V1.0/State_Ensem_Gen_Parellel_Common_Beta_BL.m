
% 更新日期：2021/05/27
% 函数功能：实现数据同化集合成员的初始化（并窗口运行）
% 输入输出变量说明
% 输入变量：S_Va_init —— 初始扰动状态变量矩阵（马赫数，迎角，卡门常数）
%           N —— 集合成员数目
%           MStepNS_cov ——  设置定常NS求解收敛最大迭代步数
%           targetPath —— 保存目标路径
% 输出变量：S_Va_init_reg —— 经过集合成员初始化后的状态变量矩阵，代入后面的内迭代环节
%           NNP —— 同化区域网格节点维度
%           NG —— 同化区域网格维度
%           MDA —— 同化工作变量维度
%           M —— 集合成员维度（状态数量）

% 主要更新：限定仅保留计算收敛的集合成员进入同化过程，修改有效集合成员 N_max → N_new

function [S_Va_init_Combined_eff,target_width,Val_Sen,N_new,NNP,NG,MDA,M]=State_Ensem_Gen_Parellel_Common_Beta_BL(S_Va_init,N_max,MStepNS_cov,targetPath)

% 定义全局变量
global savename new_output_folder fid_all fid_info PPNS initPath CFL_Setting Resi_target Resi_target_lev JX_NS STD_Comp JX_Turb L2NORM_Setting flag_Va_init;
global S_Va_turb_Name S_Va_turb_Type S_Va_turb_Num S_Va_turb_FilePath S_Va_turb_ModiLine S_Va_turb_ModiFormat S_Va_turb_Std;
global S_Va_target_Name S_Va_target_Num S_Va_target_TypeTurb S_Va_target_TypeDA S_Va_turb_target S_Va_Flow_length PT PS;
global S_Va_ref_Name S_Va_ref_Num;
global shape_upper shape_lower shape_outline shape_combined_arr_init Va_Beta Va_NNP_Beta Va_NNP_Beta_new;
global flag_da_relative flag_sigma_obs_adapt flag_sigma_obs_coeff D_obs_e_ori D_obs_e flag_flow_type;
global BL_setting N_BL L_BL BL_num BL_lab BL_dot BL_obs_num kb;
global Resi_flag Resi_flag_label;

% --------------- 修改参数预置+集合成员生成 START --------------- %

% 初始化target_width
target_width=-1;

% 初始化敏感性指标
TMu_RMSE_sum=0;

% 定义定常集合成员计算收敛条件
% Resi_target_gen_level=Resi_target_log(1)+1; % 规定集合生成截断残差水平为期望残差值 + 1（以 10 为底），以扩大集合成员捕捉范围
Resi_flag=zeros(1,N_max);  % Resi_flag_resi = 0 或 1，表示 相应编号的扰动状态计算是否收敛（是 - 1；否 - 2）；
Resi_label_pre=0;
Resi_flag_label=zeros(1,N_max);  % Resi_flag_label = k，表示对应数组编号的扰动状态计算的原始成员编号 k；

% 修改NS方程最大迭代步数：MStepNS -- NS方程最大迭代步数
% 指定替换行，对应输入文件（InputNS.dat）第9行
rep_line_MStepNS=9;
rep_format_MStepNS='%-16d';

% 修改续算开关： JX_NS -- 是否进行流场变量续算（0，1）
%                JX_Turb -- 是否进行流场变量续算（0，1）
% 指定替换行，分别对应输入文件（InputNS.dat）第10行，（InputTurb.dat）第7行
rep_line_JX_NS=10;
rep_format_JX_NS='%-12d';
rep_line_JX_Turb=7;
rep_format_JX_Turb='%-12d';

% 修改计算CFL数：CFL_Setting
% 指定替换行，对应输入文件（InputNS.dat）第8行
rep_line_CFL=8;
rep_format_CFL='%-12d';
CFL_Gen=CFL_Setting(1);

% 修改标准算例计算结果比较开关：STD_Comp
% 指定替换行，对应输入文件（InputNS.dat）第22行
rep_line_STD_Comp=22;
rep_format_STD_Comp='%-12d';
STD_Comp_Gen=STD_Comp(1);

% 修改期望残差： Resi_target
rep_line_Resi=15;
rep_format_Resi='%-10.2e';
Resi_Gen=Resi_target(1);

% 修改期望 L2NORM 水平： L2NORM_Setting
rep_line_L2NORM=8;
rep_format_L2NORM='%-12d';
L2NORM_Gen=L2NORM_Setting(1);

addPath_pre=[targetPath,'Pre_Ensemble_Generation'];

% 循环运行流场求解程序，实现修改参数+集合成员生成
% 初始化计算第k个集合成员
k_flag=0;
% run_flag(i)   保存第i个求解器当前计算的状态编号
%           =0  代表当前求解器处于空闲状态
run_flag=zeros(1,PPNS);
cnt_flag=zeros(1,PPNS);
% stop_flag(i) = 0 代表第i个求解器空闲
%              = 1 代表第i个求解器正在工作
stop_flag=zeros(1,PPNS);
% exit_flag = 0 代表并行组程序未结束
%           = 1 代表并行组程序已结束
exit_flag=0;

for kp=1:PPNS
    % 步骤①：检测 target_run_Path 目录下是否存在 ProgramExitFlag.txt
    %         若存在，则首先删除
    % 清理程序运行结束标志文件模块
    % 设置读写路径
    deletePath=['Program_CFD_Parallel_',num2str(kp,'%03d')];
    if (exist([deletePath,'\ProgramExitFlag.txt'],'file'))~=0
        delete([deletePath,'\ProgramExitFlag.txt']);
    end
end

while (exit_flag~=1)
    
    % 设置标志，若找不到空闲的求解器，则MATLAB程序暂停扫描数秒，提升计算效率
    % 初始状态所有求解器空闲
    % for i=1:PPNS
    %    stop_flag(i)=0;
    % end
    
    for kp=1:PPNS
    
        % 设置当前状态读写路径
        target_runPath=['Program_CFD_Parallel_',num2str(kp,'%03d')];
        sourcePath=[target_runPath,'\Output'];

        % 步骤②：执行并行窗口模块
        if ((run_flag(kp)==0) && (k_flag<N_max+1))
            
            k_flag=k_flag+1;
            if (k_flag==N_max+1)
                break;
            end
            run_flag(kp)=k_flag;  % 当前计算状态实际集合内编号 run_flag(kp)，即 k_flag
            
            % 求解参数修改
            % 1. 更新（关闭）续算操纵开关
            % JX_NS=0;
            % JX_Turb=0;
            fid0=fopen([target_runPath,'\Input\InputNS.dat'],'r+');
            for i=1:(rep_line_JX_NS-1)
               fgetl(fid0);
            end
            fseek(fid0,0,'cof');
            fprintf(fid0,rep_format_JX_NS,JX_NS(1));
            fclose(fid0);
            fid1=fopen([target_runPath,'\Input\InputTurb.dat'],'r+');
            for i=1:(rep_line_JX_Turb-1)
               fgetl(fid1);
            end
            fseek(fid1,0,'cof');
            fprintf(fid1,rep_format_JX_Turb,JX_Turb(1));
            fclose(fid1);
            
            % 2. 设置NS方程最大运行迭代步数
            fid2=fopen([target_runPath,'\Input\InputNS.dat'],'r+');
            for i=1:(rep_line_MStepNS-1)
               fgetl(fid2);
            end
            fseek(fid2,0,'cof');
            fprintf(fid2,rep_format_MStepNS,MStepNS_cov);
            fclose(fid2);
            
            % 3. 设置初始扰动变量
            if (flag_Va_init==0)
                % 若不进行修正项 Beta 同化
                for ns=1:S_Va_turb_Num
                    fid3=fopen([target_runPath,char(S_Va_turb_FilePath(ns))],'r+');
                    for i=1:(S_Va_turb_ModiLine(ns)-1)
                       fgetl(fid3);
                    end
                    fseek(fid3,0,'cof');
                    fprintf(fid3,char(S_Va_turb_ModiFormat(ns)),S_Va_init(run_flag(kp),ns));                
                    fclose(fid3);
                end
            else
                % 若进行修正项 Beta 同化，则写入相应信息（全场网格单元与节点 Beta 值；注意：本框架对网格单元的 Beta 进行同化）
                Beta_init=Va_Beta(:,run_flag(kp));
                save([sourcePath,'\Beta_plus.PLT'],'Beta_init','-ascii');
                Beta_init_NNP=Va_NNP_Beta(:,run_flag(kp));
                save([sourcePath,'\Beta_plus_NNP.PLT'],'Beta_init_NNP','-ascii');
            end
            
            % 4. 设置计算CFL数
            fid4=fopen([target_runPath,'\Input\InputNS.dat'],'r+');
            for i=1:(rep_line_CFL-1)
               fgetl(fid4);
            end
            fseek(fid4,0,'cof');
            fprintf(fid4,rep_format_CFL,CFL_Gen);
            fclose(fid4);
            
            % 5. 设置与标准算例比较开关
            fid5=fopen([target_runPath,'\Input\InputNS.dat'],'r+');
            for i=1:(rep_line_STD_Comp-1)
               fgetl(fid5);
            end
            fseek(fid5,0,'cof');
            fprintf(fid5,rep_format_STD_Comp,STD_Comp_Gen);
            fclose(fid5);
            
            % 6. 设置计算收敛精度
            fid6=fopen([target_runPath,'\Input\InputNS.dat'],'r+');
            for i=1:(rep_line_Resi-1)
               fgetl(fid6);
            end
            fseek(fid6,0,'cof');
            fprintf(fid6,rep_format_Resi,Resi_Gen);
            fclose(fid6);
            
            % 7. 设置L2NORM收敛水平
            fid7=fopen([target_runPath,'\Input\InputTurb.dat'],'r+');
            for i=1:(rep_line_L2NORM-1)
               fgetl(fid7);
            end
            fseek(fid7,0,'cof');
            fprintf(fid7,rep_format_L2NORM,L2NORM_Gen);
            fclose(fid7);

            % 8. 运行 CFD 程序，生成初始集合成员样本并保存
            % 运行CFD求解主程序
            fprintf(fid_all,'**************************************************\n');
            fprintf(fid_all,'|                 PROGRAM STATUS                  \n');
            fprintf(fid_all,'**************************************************\n');
            fprintf(fid_all,'|         PRE-ENSEMBLE GENERATOR PROCEDURE        \n');
            fprintf(fid_all,'**************************************************\n');
            fprintf(fid_all,'|                                                 \n');
            fprintf(fid_all,'|    The calculation of %03d sample is RUNNING!   \n',run_flag(kp));
            fprintf(fid_all,'|                                                 \n');
            fprintf(fid_all,'**************************************************\n');
            fprintf(fid_all,'\n');
                
            fprintf(['集合成员生成过程：状态 ',num2str(run_flag(kp),'%03d')]);
            fprintf(fid_info,['集合成员生成过程：状态 ',num2str(run_flag(kp),'%03d')]);
            if (flag_Va_init==0)
                % 若不进行修正项 Beta 同化
                for ns=1:S_Va_turb_Num
                    fprintf([' ',char(S_Va_turb_Name(ns)),' = %.5f'],S_Va_init(run_flag(kp),ns));
                    fprintf(fid_info,[' ',char(S_Va_turb_Name(ns)),' = %.5f'],S_Va_init(run_flag(kp),ns));
                end
            else
                % 若进行修正项 Beta 同化
                fprintf(' 初始扰动 Beta 场');
                fprintf(fid_info,' 初始扰动 Beta 场');
            end
            fprintf(' 已启动计算\n');
            fprintf(fid_info,' 已启动计算\n');

            % 运行CFD程序，启动并行窗口（并且最小化），分组运行程序
            % cmd=[target_runPath,'\CFD_2D.exe'];
            % open(cmd);
            [s,e]=dos(['cd ',target_runPath,' && run_CFD_2D.bat&']);
            stop_flag(kp)=1;
            pause(3); % 暂停数秒
            
        end
        
        % 步骤③：判断是否程序已经运行完毕
        % 若找到 ProgramExitFlag.txt，则程序已经运行完毕，进行保存
        if (exist([target_runPath,'\ProgramExitFlag.txt'],'file'))~=0
            
            % 保存初始集合成员信息
            % 保存成员指定单元处的涡粘系数

            % -------------------- 进行结果保存 START ------------------- %
            
            % 提取当前状态计算最终残差值 Resi(dual)_Final
            if (JX_Turb(1)==0)
                Resi_All=load([target_runPath,'\Output\Residual_1.PLT']);
            elseif ((JX_Turb(1)==1) || (JX_Turb(1)==2))
                Resi_All=load([target_runPath,'\Output\Residual_1_Turb_Const.PLT']);
            elseif (JX_Turb(1)==3)
                Resi_All=load([target_runPath,'\Output\Residual_1_Turb_Beta.PLT']);
            end
            Resi_Final=Resi_All(size(Resi_All,1),2);
            % 导入计算数据
            com_data=load([sourcePath,'\CP_CF_Comp.PLT']);
            % 指定调整后计算数据保存路径
            fid_com_adj_path=[sourcePath,'\CP_CF_Comp_New.PLT'];
            fid_tar_new_path=[sourcePath,'\CP_New.PLT'];
        
            % 对计算压力系数曲线节点编号进行梳理，并实时输出当前同化轮次压力分布曲线
            % 集合生成过程、内迭代过程中全部成员完成计算后保存曲线包络线
            if (flag_flow_type==0)
                [shape_upper_num,shape_lower_num,shape_upper_arr_sor,shape_lower_arr_sor,n_com]=DA_Va_Airfoil_Ref_Adjustment_and_Visulization(com_data,fid_com_adj_path,fid_tar_new_path);
            elseif (flag_flow_type==1)
                [shape_num,shape_arr_sor,n_com]=DA_Va_NonAirfoil_Ref_Adjustment_and_Visulization(com_data,fid_com_adj_path,fid_tar_new_path);
            end
            
            % 输出命令行提示
            fprintf(fid_all,'**************************************************\n');
            fprintf(fid_all,'|                 PROGRAM STATUS                  \n');
            fprintf(fid_all,'**************************************************\n');
            fprintf(fid_all,'|         PRE-ENSEMBLE GENERATOR PROCEDURE        \n');
            fprintf(fid_all,'**************************************************\n');
            fprintf(fid_all,'|                                                 \n');
            fprintf(fid_all,'|   The calculation of %03d sample is FINISHED!   \n',run_flag(kp));
            fprintf(fid_all,'|                                                 \n');
            fprintf(fid_all,'**************************************************\n');
            fprintf(fid_all,'\n');
            fprintf(fid_all,'\n');
            
            fprintf(['集合成员生成过程：状态 ',num2str(run_flag(kp),'%03d'),' 已保存完毕，Total Step = %d，Final Residual = %.5f\n'],(size(Resi_All,1)-1),Resi_Final);
            fprintf(fid_info,['集合成员生成过程：状态 ',num2str(run_flag(kp),'%03d'),' 已保存完毕，Total Step = %d，Final Residual = %.5f\n'],(size(Resi_All,1)-1),Resi_Final);
            
            % 标记计算达到收敛的算例编号 Resi_flag = 1
            if (Resi_Final<=Resi_target_lev(1))
                Resi_flag(run_flag(kp))=1; % 收敛状态设为 1，非收敛状态设为 0
                Resi_label_pre=Resi_label_pre+1; % 收敛状态数 +1
                Resi_flag_label(Resi_label_pre)=run_flag(kp); % 保存当前收敛状态编号 Resi_flag_label，累计数 = Resi_label_pre，总数 = N_new
            end
            
            if (strcmp(char(S_Va_ref_Name(1)),'BL')==1)
                
                % 运行提取边界层信息批处理文件
                fprintf(['集合成员生成过程：状态 ',num2str(run_flag(kp),'%03d'),' 边界层数据提取开始\n']);
                fprintf(fid_info,['集合成员生成过程：状态 ',num2str(run_flag(kp),'%03d'),' 边界层数据提取开始\n']);
                [s,e]=dos(['cd ',target_runPath,' && run_BL_Hump_Extract_Tecplot.bat']);
                fprintf(['集合成员生成过程：状态 ',num2str(run_flag(kp),'%03d'),' 边界层数据提取完毕\n']);
                fprintf(fid_info,['集合成员生成过程：状态 ',num2str(run_flag(kp),'%03d'),' 边界层数据提取完毕\n']);
                pause(3); % 暂停数秒
                
                % 提取边界层速度型数据
                cnt_BL_cal=0;
                BL_Cal_Vel=zeros(kb,1);
                for i=1:N_BL
                    pre_BL=[sourcePath,'\BL_Output\BL_Hump_xc_',num2str(BL_lab(i),'%03d'),'.dat'];
                    BL_info=load(pre_BL);
                    BL_Vel=BL_info(:,3); % 提取完整边界层速度信息
                    for j=1:BL_num(i)
                        cnt_BL_cal=cnt_BL_cal+1;
                        BL_Cal_Vel(cnt_BL_cal)=(BL_Vel(BL_dot(i,j)))/0.1; % 提取用于同化的边界层速度信息
                    end
                end
                if (sum(cnt_flag)==0)
                    BLSC=zeros(length(BL_Cal_Vel),N_max);
                    BLSC(:,run_flag(kp))=BL_Cal_Vel;
                else
                    BLSC(:,run_flag(kp))=BL_Cal_Vel;
                end
                
            end
             
            % 输出文件另存（防止覆盖）
            pre_Path=[addPath_pre,'\',initPath,'_',num2str(run_flag(kp),'%03d')];
            fileList=dir(sourcePath);
            mkdir(pre_Path); % 利用 mkdir() 函数创建新文件夹
            filename_s=cell(length(fileList),1);
            filename_t=cell(length(fileList),1);

            for i=3:length(fileList)
                  filename_s{i}=[sourcePath,'\',fileList(i).name];
                  filename_t{i}=pre_Path;
                  copyfile(filename_s{i},filename_t{i});
            end
            
%             if (strcmp(char(S_Va_ref_Name(1)),'BL')==1)
%                 % 提取边界层速度型数据
%                 cnt_BL_cal=0;
%                 BL_Cal_Vel=zeros(kb,1);
%                 for i=1:N_BL
%                     pre_BL=[pre_Path,'\BL_Output\BL_Hump_xc_',num2str(BL_lab(i),'%03d'),'.dat'];
%                     BL_info=load(pre_BL);
%                     BL_Vel=BL_info(:,3); % 提取完整边界层速度信息
%                     for j=1:BL_num(i)
%                         cnt_BL_cal=cnt_BL_cal+1;
%                         BL_Cal_Vel(cnt_BL_cal)=(BL_Vel(BL_dot(i,j)))/0.1; % 提取用于同化的边界层速度信息
%                     end
%                 end
%                 if (sum(cnt_flag)==0)
%                     BLSC=zeros(length(BL_Cal_Vel),N_max);
%                     BLSC(:,run_flag(kp))=BL_Cal_Vel;
%                 else
%                     BLSC(:,run_flag(kp))=BL_Cal_Vel;
%                 end
%             end
            
            % 记录当前保存文件夹路径，添加所需文件名称
            pre_VcPath=[pre_Path,'\VcData.PLT'];
            pre_TMuPath=[pre_Path,'\TMuData.PLT'];
            pre_TMuPPath=[pre_Path,'\TMuPData.PLT'];
            pre_Cp_Cf_yplus_Path=[pre_Path,'\CP_CF_Comp_New.PLT']; % 注意：已经修改为调整顺序后的计算结果文件
            pre_Beta=[pre_Path,'\Beta_plus.PLT'];
            pre_ErrInfo=[pre_Path,'\FlowField_Err_Info.txt'];

            % 提取湍流模型数据变量 Vc 与 TMu
            VcS=load(pre_VcPath);
            TMuS=load(pre_TMuPath);
            TMuPS=load(pre_TMuPPath);
            BetaS=load(pre_Beta);
            Err_val_all=load(pre_ErrInfo);
            
            Cp_Cf_yplus_S=load(pre_Cp_Cf_yplus_Path);
            Xp=Cp_Cf_yplus_S(:,1);
            Yp=Cp_Cf_yplus_S(:,2);
            CpS=Cp_Cf_yplus_S(:,3);
            CfS=Cp_Cf_yplus_S(:,4);
            yplusS=Cp_Cf_yplus_S(:,5);
            
            TMu_MSE=Err_val_all(1,7);
            TMu_RMSE=Err_val_all(2,7);
            % 若计算收敛才计算敏感性指标
            if (Resi_flag(run_flag(kp))==1)
                TMu_RMSE_sum=TMu_RMSE_sum+TMu_RMSE;
            end

            % 保存顺序需要特别注意！！！待修改
            if (sum(cnt_flag)==0)
                VcSC=zeros(length(VcS),N_max);
                TMuSC=zeros(length(TMuS),N_max);
                TMuPSC=zeros(length(TMuPS),N_max);
                BetaSC=zeros(length(BetaS),N_max);
                CpSC=zeros(length(CpS),N_max);
                CfSC=zeros(length(CfS),N_max);
                yplusSC=zeros(length(yplusS),N_max);
                
                VcSC(:,run_flag(kp))=VcS;
                TMuSC(:,run_flag(kp))=TMuS;
                TMuPSC(:,run_flag(kp))=TMuPS;
                BetaSC(:,run_flag(kp))=BetaS;
                CpSC(:,run_flag(kp))=CpS;
                CfSC(:,run_flag(kp))=CfS;
                yplusSC(:,run_flag(kp))=yplusS;
            else
                VcSC(:,run_flag(kp))=VcS;
                TMuSC(:,run_flag(kp))=TMuS;
                TMuPSC(:,run_flag(kp))=TMuPS;
                BetaSC(:,run_flag(kp))=BetaS;
                
                CpSC(:,run_flag(kp))=CpS;
                CfSC(:,run_flag(kp))=CfS;
                yplusSC(:,run_flag(kp))=yplusS;
            end

            % -------------------- 进行结果保存 END -------------------- %
            
            run_flag(kp)=0;  % 将求解器重置为空闲状态
            delete([target_runPath,'\ProgramExitFlag.txt']);
            cnt_flag(kp)=cnt_flag(kp)+1; % 累计各组计算个数
            stop_flag(kp)=0;
            pause(3); % 保存文件后，暂停数秒
            
        end
        
        if (sum(stop_flag)==PPNS)
            % 找不到空闲的求解器，MATLAB程序暂停扫描数秒
            pause(10);
        end
        
        if ((sum(stop_flag)==0) && (sum(cnt_flag)==N_max))
            % 若求解器均空闲，且满足退出条件，则完成并行程序组运行
            exit_flag=1;
            pause(10);
            break;
        end
    
    end 
    
end

% 所有状态均计算结束并保存后，计算新的集合成员个数 N_new
N_new=sum(Resi_flag);
% 将收敛标记数组 Resi_flag_label 的尾端无关项删去
Resi_flag_label((N_new+1):end)=[];
% 计算该组计算结果敏感性指标
Val_Sen=TMu_RMSE_sum/N_new;

% !!!!!!!!!!!!!!!!!! 须修改模块，对每部分状态变量根据计算收敛情况进行修改
% % 保存初始更新后状态变量
% S_Va_init_reg=[(S_Va_init(:,(1:PT))),TMuSC',CpSC']';
% NNP=length(TMuPS);   % 计算全流场网格节点维度
% NG=length(TMuS);   % 计算全流场网格维度
% MDA=length(TMuS)+PT;   % 计算全流场工作变量维度
% M=length(S_Va_init_reg);  % 计算集合成员维度（状态数量）

% 同化状态变量由三部分组成：① 单同化变量；② 流场同化变量；③ 观测变量
% 以下分别命名为 S_Va_init_Part_One(Two,Three)
% 最终得到的合并后的同化状态变量命名为 S_Va_init_combined

% ① 提取单同化变量
S_Va_cnt_One=0;
S_Va_init_Part_One=zeros(N_max,PT);
for i=1:S_Va_target_Num
    if ((S_Va_target_TypeTurb(i)>0) && (S_Va_target_TypeDA(i)==1))
        S_Va_cnt_One=S_Va_cnt_One+1;
        S_Va_init_Part_One(:,S_Va_cnt_One)=S_Va_init(:,i);
    end
end
if (S_Va_cnt_One==0)
    S_Va_init_Part_One=[];
end

% ② 提取流场同化变量
S_Va_cnt_Two=0;
% 按顺序保存各流场同化变量维度（长度）
S_Va_Flow_length=zeros(1,sum((S_Va_target_TypeDA==2)));
for i=1:S_Va_target_Num
    if (S_Va_target_TypeDA(i)==2)
        S_Va_cnt_Two=S_Va_cnt_Two+1;
        if (S_Va_cnt_Two==1)
            if (strcmp(char(S_Va_target_Name(i)),'Vc')==1)
                S_Va_init_Part_Two=VcSC';
                S_Va_Flow_length(S_Va_cnt_Two)=length(VcS);
            elseif (strcmp(char(S_Va_target_Name(i)),'TMu')==1)
                S_Va_init_Part_Two=TMuSC';
                S_Va_Flow_length(S_Va_cnt_Two)=length(TMuS);
            elseif (strcmp(char(S_Va_target_Name(i)),'TMuP')==1)
                S_Va_init_Part_Two=TMuPSC';
                S_Va_Flow_length(S_Va_cnt_Two)=length(TMuPS);
           elseif (strcmp(char(S_Va_target_Name(i)),'Beta')==1)
                S_Va_init_Part_Two=BetaSC';
                S_Va_Flow_length(S_Va_cnt_Two)=length(BetaS);
            end
        else
            if (strcmp(char(S_Va_target_Name(i)),'Vc')==1)
                S_Va_init_Part_Two_F=[S_Va_init_Part_Two,VcSC'];
                S_Va_init_Part_Two=S_Va_init_Part_Two_F;
                S_Va_Flow_length(S_Va_cnt_Two)=length(VcS);
            elseif (strcmp(char(S_Va_target_Name(i)),'TMu')==1)
                S_Va_init_Part_Two_F=[S_Va_init_Part_Two,TMuSC'];
                S_Va_init_Part_Two=S_Va_init_Part_Two_F;
                S_Va_Flow_length(S_Va_cnt_Two)=length(TMuS);
            elseif (strcmp(char(S_Va_target_Name(i)),'TMuP')==1)
                S_Va_init_Part_Two_F=[S_Va_init_Part_Two,TMuPSC'];
                S_Va_init_Part_Two=S_Va_init_Part_Two_F;
                S_Va_Flow_length(S_Va_cnt_Two)=length(TMuPS);
            elseif (strcmp(char(S_Va_target_Name(i)),'Beta')==1)
                S_Va_init_Part_Two_F=[S_Va_init_Part_Two,BetaSC'];
                S_Va_init_Part_Two=S_Va_init_Part_Two_F;
                S_Va_Flow_length(S_Va_cnt_Two)=length(BetaS);
            end
        end
    end
end
if (S_Va_cnt_Two==0)
    S_Va_init_Part_Two=[];
end

% ③ 提取同化观测变量
S_Va_cnt_Three=0;
for i=1:S_Va_ref_Num
    S_Va_cnt_Three=S_Va_cnt_Three+1;
    if (S_Va_cnt_Three==1)
        if (strcmp(char(S_Va_ref_Name(i)),'Cp')==1)
            S_Va_init_Part_Three=CpSC';
        elseif (strcmp(char(S_Va_ref_Name(i)),'Cf')==1)
            S_Va_init_Part_Three=CfSC';
        elseif (strcmp(char(S_Va_ref_Name(i)),'yplus')==1)
            S_Va_init_Part_Three=yplusSC';
        elseif (strcmp(char(S_Va_ref_Name(1)),'BL')==1)
            S_Va_init_Part_Three=BLSC';
        end
    else
        if (strcmp(char(S_Va_ref_Name(i)),'Cp')==1)
            S_Va_init_Part_Three_F=[S_Va_init_Part_Three,CpSC'];
            S_Va_init_Part_Three=S_Va_init_Part_Three_F;
        elseif (strcmp(char(S_Va_ref_Name(i)),'Cf')==1)
            S_Va_init_Part_Three_F=[S_Va_init_Part_Three,CfSC'];
            S_Va_init_Part_Three=S_Va_init_Part_Three_F;
        elseif (strcmp(char(S_Va_ref_Name(i)),'yplus')==1)
            S_Va_init_Part_Three_F=[S_Va_init_Part_Three,yplusSC'];
            S_Va_init_Part_Three=S_Va_init_Part_Three_F;
        elseif (strcmp(char(S_Va_ref_Name(1)),'BL')==1)
            S_Va_init_Part_Three_F=[S_Va_init_Part_Three,BLSC'];
            S_Va_init_Part_Three=S_Va_init_Part_Three_F;
        end
    end
end
if (S_Va_cnt_Three==0)
    S_Va_init_Part_Three=[];
end

% 保存计算数据信息
save([targetPath,'Com_Ensem_Gen.dat'],'S_Va_init_Part_Three','-ascii');

if ((strcmp(char(S_Va_ref_Name(1)),'Cp')==1) || (strcmp(char(S_Va_ref_Name(1)),'Cf')==1) || (strcmp(char(S_Va_ref_Name(1)),'yplus')==1))

    % !!!!! 修改部分，提取有效部分的观测信息 !!!!!
    S_Va_init_Part_Three_eff=zeros(N_new,size(S_Va_init_Part_Three,2));
    n_combined_eff=0;
    for i=1:N_max
        if (Resi_flag(i)==1)
            n_combined_eff=n_combined_eff+1; % 收敛判断计数器加一
            S_Va_init_Part_Three_eff(n_combined_eff,:)=S_Va_init_Part_Three(i,:);
        end
    end

    % 对所有集合状态计算得到的目标同化变量（分布）的包络线进行绘制
    S_Va_init_Part_Three_min=min(S_Va_init_Part_Three_eff); % 下边界
    S_Va_init_Part_Three_max=max(S_Va_init_Part_Three_eff); % 上边界
    S_Va_init_Part_Three_mid=(S_Va_init_Part_Three_min+S_Va_init_Part_Three_max)/2; % 上下包络线中线
    target_width=S_Va_init_Part_Three_max-S_Va_init_Part_Three_min; % 上下包络边界宽度

    if (flag_flow_type==0)

        y1=S_Va_init_Part_Three_max(1:shape_lower_num);
        y2=S_Va_init_Part_Three_min(1:shape_lower_num);
        y3=S_Va_init_Part_Three_max(shape_lower_num+1:end);
        y4=S_Va_init_Part_Three_min(shape_lower_num+1:end);
        % 新建可视化图层
        % !!! 增加包络区域中线 !!! Middle_line
        figure;
        hold on
        [hha1,hhb1,hhc1]=shadedplot(Xp(1:shape_lower_num),y1,y2,'[0.8706 0.9216 0.9804]');
        [hha2,hhb2,hhc2]=shadedplot(Xp(shape_lower_num+1:end),y3,y4,'[0.8706 0.9216 0.9804]');
        hh1=plot(shape_combined_arr_init(:,2),shape_combined_arr_init(:,7),'ro','MarkerSize',4.0); % 绘制标准实验值
        hh2=plot(shape_combined_arr_init(:,2),S_Va_init_Part_Three_mid,'-','LineWidth',0.8,'Color',[0.7451,0.7451,0.7451]); % 绘制包络区间中线
        % plot(Xp,S_Va_init_Part_Three_max); % 绘制上区间包络线
        % plot(Xp,S_Va_init_Part_Three_min); % 绘制下区间包络线
        legend([hhb1,hhc1,hha1(2),hh2,hh1],{'计算上区间包络线','计算下区间包络线','上下边界所围区域','计算包络区间中线','实验值测量值'}); % 添加部分图例
        xlabel('\itx'),ylabel('\itCp'),title('集合生成计算结果');
        grid on
        box on
        hold off
        % 保存可编辑图片
        saveas(gcf,['.\Output_Figures\',savename,'\Init_Ensem_Gen_Output.fig']);

    elseif (flag_flow_type==1)

        y1=S_Va_init_Part_Three_max;
        y2=S_Va_init_Part_Three_min;
        % 新建可视化图层
        % !!! 增加包络区域中线 !!! Middle_line
        figure;
        hold on
        [hha1,hhb1,hhc1]=shadedplot(Xp,y1,y2,'[0.8706 0.9216 0.9804]');
        hh1=plot(shape_combined_arr_init(:,2),shape_combined_arr_init(:,7),'ro','MarkerSize',4.0); % 绘制标准实验值
        hh2=plot(shape_combined_arr_init(:,2),S_Va_init_Part_Three_mid,'-','LineWidth',0.8,'Color',[0.7451,0.7451,0.7451]); % 绘制包络区间中线
        legend([hhb1,hhc1,hha1(2),hh2,hh1],{'计算上区间包络线','计算下区间包络线','上下边界所围区域','计算包络区间中线','实验值测量值'}); % 添加部分图例
        xlabel('\itx'),ylabel('\itCp'),title('集合生成计算结果');
        grid on
        box on
        hold off
        % 保存可编辑图片
        saveas(gcf,['.\Output_Figures\',savename,'\Init_Ensem_Gen_Output.fig']);

    end
    
end

% 保存初始更新后全状态变量（三部分叠加）
if (flag_da_relative==1)
    S_Va_init_Combined=[S_Va_init_Part_One,S_Va_init_Part_Two,S_Va_init_Part_Three]';
elseif (flag_da_relative==0)
    % !!!! 2021/03/03构思模块测试 !!!!
    % 对观测变量的值进行转换，目标是得到与实验值（真值）的相对误差值
    % 可将实验值（真值）设置为零
    % 同化目标即转变为减小观测值与实验值（真值），即零的相对误差
    S_Va_init_Part_Three_err=(S_Va_init_Part_Three./D_obs_e_ori')-ones(size(D_obs_e_ori'));
    S_Va_init_Combined=[S_Va_init_Part_One,S_Va_init_Part_Two,S_Va_init_Part_Three_err]';
end

% -- 备份操作 -- %
save([new_output_folder,savename,'_S_Va_init_Combined.mat'],'S_Va_init_Combined'); % 保存生成原始成员结果至工作区
% -- 备份操作 -- %

% % !!!!! 修改部分，提取有效部分的状态变量 !!!!!
% S_Va_init_Combined_eff=zeros(N_new,size(S_Va_init_Combined,1));
% n_combined_eff=0;
% for i=1:N_max
%     if (Resi_flag(i)==1)
%         n_combined_eff=n_combined_eff+1;
%         S_Va_init_Combined_eff(n_combined_eff,:)=S_Va_init_Combined(i,:);
%     end
% end

% !!!!! 修改部分，提取有效部分的状态变量 !!!!!
S_Va_init_Combined_eff=zeros(size(S_Va_init_Combined,1),N_new);
n_combined_eff=0;
for i=1:N_max
    if (Resi_flag(i)==1)
        n_combined_eff=n_combined_eff+1;
        S_Va_init_Combined_eff(:,n_combined_eff)=S_Va_init_Combined(:,i);
    end
end

% -- 备份操作 -- %
save([new_output_folder,savename,'_S_Va_init_Combined_eff.mat'],'S_Va_init_Combined_eff'); % 保存生成原始成员结果至工作区
% -- 备份操作 -- %

% 保存初始更新后全状态变量 - 特征数值
NNP=length(TMuPS);   % 计算全流场网格节点维度
NG=length(TMuS);   % 计算全流场网格维度
% 同化相关维度计算，并传递出函数
MDA=PT+size(S_Va_init_Part_Two,2);   % 计算全部同化工作变量维度（单值变量 + 流场多值变量）
M=size(S_Va_init_Combined_eff,1);  % 计算集合成员维度（状态数量）

% --------------- 修改参数预置+集合成员生成 END --------------- %

fprintf('\n');
fprintf(fid_info,'\n');

end
