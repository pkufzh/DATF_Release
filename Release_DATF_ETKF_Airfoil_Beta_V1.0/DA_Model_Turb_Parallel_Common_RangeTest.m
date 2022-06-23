%% 数据同化湍流计算辅助程序（Data Assimilation & Turbulent Flows）
%% 集合变换卡尔曼滤波实验（Ensemble Transform Kalman Filter（ETKF））
%% 利用ETKF方法结合实验数据，实现对翼型绕流等算例的数据同化
%% 程序主要功能：能够针对期望变量进行指定范围、指定步长的扰动，并保存标准计算值（算例）
%%              有助于对指定算例的初始集合扰动区间设计

%% 初始化程序
clear all
clc
close all

% 导入初始成员数据
% load('Pre_Cal_Path.mat');

format long;

% 定义全局变量
global savename fid_all fid_info PPNS initPath CFL_Setting;
global S_Va_turb_Name S_Va_turb_Num S_Va_turb_FilePath S_Va_turb_ModiLine S_Va_turb_ModiFormat S_Va_turb_Std S_Va_turb_Det_base_mark;
global S_Va_target_Name S_Va_target_Num S_Va_target_TypeTurb S_Va_target_TypeDA S_Va_turb_target S_Va_Flow_length PT PS;
global S_Va_ref_Name S_Va_ref_Num;
global shape_upper shape_lower shape_combined_arr_init;

%% 设置主要数据同化过程参数

% ------------  设置主要文件交互与并行相关参数 START ------------ %

% 设置辅助程序数据交互路径
savename='Results_Model_Turb_Parallel_NACA4412_RangeTest_Cb1_Sigma_Cv1_CKS';
new_output_folder=['.\',savename,'\'];
mkdir(new_output_folder);
targetPath=new_output_folder;
initPath='Init_Ensem_Member';

% 打开记录程序信息文件
fid_info=fopen([new_output_folder,'DA_CompInfo_Brief.dat'],'w+');
fid_all=fopen([new_output_folder,'DA_CompInfo_Detailed.dat'],'w+');

% 设置并行线程/分组数（可根据工作机线程设置）
PPNS=7;

% 导入翼型原始外形数据，并进行上下翼型插值
shape_upper=load('.\Input_Data\Input_Data_Shape_NACA4412_upper.dat');
shape_lower=load('.\Input_Data\Input_Data_Shape_NACA4412_lower.dat');
% 导入上下翼面的实验数据
exp_upper_data=load('.\Input_Data\Input_Data_Exp_NACA4412_upper.dat');
exp_lower_data=load('.\Input_Data\Input_Data_Exp_NACA4412_lower.dat');

% -------------  设置主要文件交互与并行相关参数 END ------------- %

% -------------  设置生成初始集合成员相关参数 START ------------- %

% 定义集合成员数量
N=100;

% ----------  指定初始扰动变量、目标同化变量、同化观测变量模块 START ---------- %

% ① 指定初始扰动变量模块（目前默认为单变量扰动）
% 算例标准状态设置
% 定义初始扰动变量名称、文件路径、修改行数、修改格式（均用数组表示，一一对应）
% 初始扰动变量名称
S_Va_turb_Name=["Cb1","Sigma","Cv1","CKS"];
S_Va_turb_CharName=char(S_Va_turb_Name);
% 初始扰动变量个数
S_Va_turb_Num=length(S_Va_turb_Name);
% 初始扰动变量文件路径（省略运行程序根目录）
% S_Va_turb_FilePath_raw=char('\Input\InputNS.dat','\Input\InputNS.dat','\Input\InputTurb.dat');
% S_Va_turb_FilePath=strip(S_Va_turb_FilePath_raw,'right');
S_Va_turb_FilePath=["\Input\InputSAconst.dat","\Input\InputSAconst.dat","\Input\InputSAconst.dat","\Input\InputSAconst.dat"];
% 初始扰动变量修改行数（对应文件内）
S_Va_turb_ModiLine=[3,5,6,7];
% 初始扰动变量修改格式（对应文件内）
S_Va_turb_ModiFormat=["%-16.5f","%-16.5f","%-16.5f","%-16.5f"];
% 初始扰动变量扰动区间
% 设置扰动方式：0 - 绝对值上下扰动；1 - 比例上下扰动（例如：设置50代表扰动区间为标准值的[50%,150%]）
S_Va_turb_way=[1,1,1,1];
S_Va_turb_Std=[0.1355,(2/3),7.1,0.41];  % 初始扰动变量标准值
S_Va_turb_Det_base_lower=[10,40,80,10];  % 初始扰动变量下区间
S_Va_turb_Det_base_upper=[80,80,80,150];  % 初始扰动变量上区间

% S_Gap=8; % 定义最宽区间状态扰动步长（百分比，即数值 5 代表具有最大扰动区间的变量以 5% 步长变化）
% S_Va_turb_Det_base_range=S_Va_turb_Det_base_lower+S_Va_turb_Det_base_upper;
% N=round(max(S_Va_turb_Det_base_range)/S_Gap)+1;
% S_Va_turb_Det_base_ratio=S_Va_turb_Det_base_range./(N-1); % 计算每一个变量实际的等步长变化率（百分比）
% S_Va_turb_Det_base_mark=zeros(N,S_Va_turb_Num);
% for i=1:S_Va_turb_Num
%     S_Va_turb_Det_base_mark(:,i)=linspace(((-1)*S_Va_turb_Det_base_lower(i)),S_Va_turb_Det_base_upper(i),N);
% end

S_Va_turb_Det_real_lower=S_Va_turb_Det_base_lower;
S_Va_turb_Det_real_upper=S_Va_turb_Det_base_upper;
for i=1:length(S_Va_turb_Std)
    % 若采用比例上下扰动方式
    if (S_Va_turb_way(i)==1)
        S_Va_turb_Det_real_lower(i)=S_Va_turb_Std(i)*(S_Va_turb_Det_base_lower(i)/100);
        S_Va_turb_Det_real_upper(i)=S_Va_turb_Std(i)*(S_Va_turb_Det_base_upper(i)/100);
    end
end
S_Va_turb_min=S_Va_turb_Std-S_Va_turb_Det_real_lower;
S_Va_turb_max=S_Va_turb_Std+S_Va_turb_Det_real_upper;

% S_Va_init=zeros(N,length(S_Va_turb_min));
% % 等步长初始扰动变量采样
% for i=1:length(S_Va_turb_min)
%     Va_min=S_Va_turb_min(i);
%     Va_max=S_Va_turb_max(i);
%     Va_turb=linspace(Va_min,Va_max,N)';
%     S_Va_init(:,i)=Va_turb;
% end
% S_Va_init_ave=mean(S_Va_init);

% 利用拉丁超立方进行初始扰动变量采样
Va_min=S_Va_turb_min;
Va_max=S_Va_turb_max;
% 初始同化变量个数
PN=length(Va_min);
S_Va_coef=lhsdesign(N,PN);
% 生成N*P列扰动变量向量
S_Va_init=(ones(N,1)*Va_min)+S_Va_coef.*(ones(N,1)*(Va_max-Va_min));
S_Va_init_ave=mean(S_Va_init);

S_Va_turb_Det_base_mark=zeros(N,S_Va_turb_Num);
for i=1:N
    for j=1:S_Va_turb_Num
        S_Va_turb_Det_base_mark(i,j)=((S_Va_init(i,j)-S_Va_turb_Std(j))/S_Va_turb_Std(j))*100;
    end
end

% ② 指定目标同化变量模块
% 全部目标同化变量名称
S_Va_target_Name=["TMu"];
S_Va_target_CharName=char(S_Va_target_Name);
% 全部目标同化变量数量
S_Va_target_Num=length(S_Va_target_Name);
% 判断目标同化变量类型①（TypeTurb）：k - 是第k个初始扰动变量（k~=0）；0 - 非初始扰动变量
% 注：若目标同化变量非初始扰动变量，则需要在更新步将值改为原始 STD 值
S_Va_target_TypeTurb=[0];
% 判断目标同化变量类型②（TypeDA）：1 - 单（值）同化变量；2 - 流场（多值）同化变量
S_Va_target_TypeDA=[2];
% 目标单值同化（且为初始扰动）的变量个数为 PT；流场（多值）同化变量个数 PS，即
% 例如在本算例中，如欲同化马赫数 Ma、迎角 AoA、全（部分）流场涡粘 TMu，则 PT = 2；
%                   欲仅同化迎角 AoA、全（部分）流场涡粘 TMu，则 PT = 1；
% PT=sum((S_Va_target_TypeTurb>0) & (S_Va_target_TypeDA==1));
PT=sum((S_Va_target_TypeDA==1));
PS=sum((S_Va_target_TypeDA==2));
% 标记初始扰动变量中同时满足[为目标同化]和[单值]的变量（e.g. "Ma" "AoA"），并保存在S_Va_target_NonTurb数组
Lia=ismember(S_Va_turb_Name,S_Va_target_Name);
S_Va_turb_target=double(Lia);

% ③ 指定同化观测变量模块
% 全部同化观测变量名称
S_Va_ref_Name=["Cp"];
S_Va_ref_CharName=char(S_Va_ref_Name);
% 全部同化观测变量数量
S_Va_ref_Num=length(S_Va_ref_Name);

% -----------  指定初始扰动变量、目标同化变量、同化观测变量模块 END ----------- %

% --------------  设置生成初始集合成员相关参数 END -------------- %

% -----------------  设置主要求解迭代参数 START ----------------- %

MStepNS_cov=50000;   % 定常NS求解收敛最大迭代步数（集合生成过程）
CFL_Setting=5; % 分别设置集合成员生成、内迭代以及重计算过程的 CFL 数
    
% ------------------  设置主要求解迭代参数 END ------------------ %

% ---------------  输出算例主要设置信息 START --------------- %

% 命令行输出提示
fprintf('******************************************************************\n');
fprintf('                      集合扰动测试主要设置信息                    \n');
fprintf('******************************************************************\n');
fprintf('[程序保存路径 & 并行设置]\n');
fprintf(['结果保存路径：',savename,'\n']);
fprintf(['程序信息文件：DA_CompInfo_Brief.dat\n']);
fprintf('程序并行线程/分组数 PPNS = %d\n',PPNS);
fprintf('******************************************************************\n');
fprintf('[同化主要参数]\n');
fprintf('集合成员数量 N = %d\n',N);
fprintf('******************************************************************\n');
fprintf('[初始扰动变量信息]\n');
for i=1:S_Va_turb_Num
    fprintf(['初始扰动变量 ',num2str(i,'%03d'),'：',S_Va_turb_CharName(:,:,i),' ∈ [ %.5f , %.5f ], STD = %.5f\n'],S_Va_turb_min(i),S_Va_turb_max(i),S_Va_turb_Std(i));
end
fprintf('******************************************************************\n');
fprintf('[同化流程控制]\n');
fprintf('设置集合生成过程最大迭代步数：%d, CFL = %d\n',MStepNS_cov,CFL_Setting(1));
fprintf('******************************************************************\n');
fprintf('\n');

% 文件输出提示
fprintf(fid_info,'******************************************************************\n');
fprintf(fid_info,'                      集合扰动测试主要设置信息                    \n');
fprintf(fid_info,'******************************************************************\n');
fprintf(fid_info,'[程序保存路径 & 并行设置]\n');
fprintf(fid_info,['结果保存路径：',savename,'\n']);
fprintf(fid_info,['程序信息文件：DA_CompInfo_Brief.dat\n']);
fprintf(fid_info,'程序并行线程/分组数 PPNS = %d\n',PPNS);
fprintf(fid_info,'******************************************************************\n');
fprintf(fid_info,'[同化主要参数]\n');
fprintf(fid_info,'集合成员数量 N = %d\n',N);
fprintf('******************************************************************\n');
fprintf(fid_info,'[初始扰动变量信息]\n');
for i=1:S_Va_turb_Num
    fprintf(fid_info,['初始扰动变量 ',num2str(i,'%03d'),'：',S_Va_turb_CharName(:,:,i),' ∈ [ %.5f , %.5f ], STD = %.5f\n'],S_Va_turb_min(i),S_Va_turb_max(i),S_Va_turb_Std(i));
end
fprintf(fid_info,'******************************************************************\n');
fprintf(fid_info,'[同化流程控制]\n');
fprintf(fid_info,'设置集合生成过程最大迭代步数：%d, CFL = %d\n',MStepNS_cov,CFL_Setting(1));
fprintf(fid_info,'******************************************************************\n');
fprintf(fid_info,'\n');

% ----------------  输出算例主要设置信息 END ---------------- %

%% 步骤一：导入参考实验与首轮计算数据，并进行插值操作
% 计算数据操作
% 导入计算数据
com_data=load('.\Program_CFD_Base\Output\CP_CF_Comp.PLT');
% 指定调整后计算数据保存路径
fid_com_adj_path='.\Program_CFD_Base\Output\CP_CF_Comp_New.PLT';
fid_tar_new_path='.\Program_CFD_Base\Output\CP_New.PLT';

% 对计算压力系数曲线节点编号进行梳理，并实时输出当前同化轮次压力分布曲线
% 集合生成过程、内迭代过程中全部成员完成计算后保存曲线包络线
[shape_upper_num,shape_lower_num,shape_upper_arr_sor,shape_lower_arr_sor,n_com]=DA_Va_Ref_Adjustment_and_Visulization(com_data,fid_com_adj_path,fid_tar_new_path);

% 实验数据操作
% 将上下翼型实验数据数组根据x坐标进行升序排序（指定[0,1]区间）
exp_upper_data_sor=sortrows(exp_upper_data,1);
exp_lower_data_sor=sortrows(exp_lower_data,1);
% 写入输出文件夹下的 InputExpData_all.dat 与 InputExpData.PLT 文件
exp_upper_data_sor_out=sortrows(exp_upper_data,1,'descend');
exp_lower_data_sor_out=sortrows(exp_lower_data,1,'ascend');
exp_data_out=[exp_upper_data_sor_out;exp_lower_data_sor_out];
fid_exp_plt=fopen('.\Program_CFD_Base\Output\InputExpData.PLT','w+');
fprintf(fid_exp_plt,'VARIABLES="X","CP"\n');
fprintf(fid_exp_plt,'ZONE T="ZONE EXP "\n');
for i=1:length(exp_data_out)
    fprintf(fid_exp_plt,'%20.8f %20.8f\n',exp_data_out(i,1),exp_data_out(i,2));
end
fclose(fid_exp_plt);
fid_exp_dat=fopen('.\Program_CFD_Base\Output\InputExpData_all.dat','w+');
for i=1:length(exp_data_out)
    fprintf(fid_exp_dat,'%20.8f %20.8f\n',exp_data_out(i,1),exp_data_out(i,2));
end
fclose(fid_exp_dat);
% 记录实验测点个数
exp_upper_num=length(exp_upper_data_sor);
exp_lower_num=length(exp_lower_data_sor);
exp_num_tot=exp_upper_num+exp_lower_num;
% 进行上下翼面实验数据值的插值
xx_upper_exp=(0:0.001:1);
yy_upper_exp=interp1(exp_upper_data_sor(:,1),exp_upper_data_sor(:,2),xx_upper_exp,'linear','extrap');
xx_lower_exp=(0:0.001:1);
yy_lower_exp=interp1(exp_lower_data_sor(:,1),exp_lower_data_sor(:,2),xx_lower_exp,'linear','extrap');

% 利用插值方法求出计算数据点坐标对应的实验值
for i=1:shape_upper_num
    x_com_upper_pre=shape_upper_arr_sor(i,2);
    shape_upper_arr_sor(i,7)=interp1(exp_upper_data_sor(:,1),exp_upper_data_sor(:,2),x_com_upper_pre,'linear','extrap');
end
for i=1:shape_lower_num
    x_com_lower_pre=shape_lower_arr_sor(i,2);
    shape_lower_arr_sor(i,7)=interp1(exp_lower_data_sor(:,1),exp_lower_data_sor(:,2),x_com_lower_pre,'linear','extrap');
end

% 合并上下翼面计算与实验值（插值）点数据
shape_combined_arr_init=[shape_lower_arr_sor;shape_upper_arr_sor]; % 先下翼面，后上翼面

% % 按计算数据点原始编号进行升序排序
% shape_combined_arr_sor=sortrows(shape_combined_arr,1);

% % 建立由实验值至计算点的映射关系
% % 记录参与同化的实验数据点信息，可导出文件
% exp_for_DA=zeros(k,2);
% for i=1:k
%     
%     exp_for_DA(i,1)=shape_combined_arr_init(num_da_arr(i),2);
%     
%     if (strcmp(char(S_Va_ref_Name(1)),'Cp')==1)
%         exp_for_DA(i,2)=shape_combined_arr_init(num_da_arr(i),7);
%     elseif (strcmp(char(S_Va_ref_Name(1)),'Cf')==1)
%         exp_for_DA(i,2)=shape_combined_arr_init(num_da_arr(i),7);
%     elseif (strcmp(char(S_Va_ref_Name(1)),'yplus')==1)
%         exp_for_DA(i,2)=shape_combined_arr_init(num_da_arr(i),7);
%     end
%     
% end
% 
% % 输出进行同化的观测实验值分布
% fid_exp_da=fopen('.\Program_CFD_Base\Output\InputExpData_for_DA.PLT','w+');
% fprintf(fid_exp_da,'VARIABLES="X","CP"\n');
% fprintf(fid_exp_da,'ZONE T="ZONE EXP for DA"\n');
% for i=1:k
%     fprintf(fid_exp_da,'%20.8f %20.8f\n',exp_for_DA(i,1),exp_for_DA(i,2));
% end
% fclose(fid_exp_da);

% 复制标准湍流模型计算结果，方便模板化显示
copyfile('.\Program_CFD_Base\Output\CP_New.PLT','.\Program_CFD_Base\Output\CP_SA.PLT');

% 初始标准状态计算结果可视化
figure;
hold on
plot(shape_combined_arr_init(:,2),shape_combined_arr_init(:,7),'ro','MarkerSize',4.0);
plot(shape_combined_arr_init(:,2),shape_combined_arr_init(:,4),'b-','LineWidth',2.0);
xlabel('\itx'),ylabel('\itCp'),title('初始标准状态计算结果');
legend('实验测量值','初始翼型压力系数分布');
grid on
box on
hold off
% 保存可编辑图片
saveas(gcf,['.\Output_Figures\',savename,'_RangeTest_Init_State_Output.fig']);

%% 步骤一：并行程序组创建

% 提取翼型网格参数
% 运行网格转换主程序
fprintf('翼型绕流数据同化辅助程序（集合成员设计）启动！\n');
fprintf(fid_info,'翼型绕流数据同化辅助程序（集合成员设计）启动！\n');
fprintf('\n');
fprintf(fid_info,'\n');
cmd_grid=('.\Program_CFD_Base\GridTreat_.exe');
open(cmd_grid);
pause(3);
fprintf('基础程序网格转换完成！\n');
fprintf(fid_info,'基础程序网格转换完成！\n');

% 创建程序运行文件夹（数量为PPNS）
for p=1:PPNS
    ddir=['Program_CFD_Parallel_',num2str(p,'%03d')];
    copyfile('Program_CFD_Base',ddir);
    % 首次运行网格转换
    cmd_grid=([ddir,'\GridTreat_.exe']);
    open(cmd_grid);
    pause(3);
    fprintf('并行程序组 %03d 网格转换完成！\n',p);
    fprintf(fid_info,'并行程序组 %03d 网格转换完成！\n',p);
end

fprintf('\n');
fprintf(fid_info,'\n');

%% 步骤三：初始化集合成员

% -------------- 修改参数预置 + 集合成员生成 START -------------- %

S_Va_init_Combined=State_Ensem_Gen_Parellel_Common_Assist_RangeTest(S_Va_init,N,MStepNS_cov,targetPath);

% 构建全状态变量矩阵，初始集合成员矩阵
S_Va_ana=S_Va_init_Combined;

% --------------- 修改参数预置 + 集合成员生成 END --------------- %

% 保存初始集合信息
% 保存可供后续直接调用或续算的工作区数据
save([new_output_folder,savename,'_Pre.mat']);  % 保存集合成员生成过程所有结果
save([new_output_folder,savename,'_S_Va_init.mat'],'S_Va_ana'); % 保存更新成员结果至工作区
save([new_output_folder,'S_Va_init.dat'],'S_Va_ana','-ascii');
save([new_output_folder,'S_Va_init.mat'],'S_Va_ana');

fclose(fid_info);
fclose(fid_all);
