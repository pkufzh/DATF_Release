%% 数据同化湍流计算程序（Data Assimilation & Turbulent Flows）
%% 集合变换卡尔曼滤波实验（Ensemble Transform Kalman Filter（ETKF））
%% 利用ETKF方法结合实验数据，实现对翼型绕流等算例的数据同化
%% 【更新】2021/06/04：主程序——指定初始扰动变量、目标同化变量，进行数据同化（通用程序 Beta_0510）
%% 【更新】2021/06/04：框架 & 流程性检查
%% 主要功能：重新生成集合成员并进行数据同化（完整）
%% 修复BUG：① 对于在内迭代过程需要恢复标准值的变量——进行判断语句修改；
%%          ② 正态分布函数设置问题 X ~ N(mu,sigma)；
%% 更新功能：
%%          ① 涡粘系数取对数方法（log）测试（有BUG）；
%%          ② 增加“设置同化对象/数量”模块；
%%          ③ 增加区域数据同化模块（适用于Beta或者TMu同化 ）；
%%          ④ 增加指定初始扰动变量、目标同化变量模块；
%%          ⑤ 增加同化滤波过程采用相对误差，替代绝对误差；
%%          ⑥ 增加对计算压力系数曲线节点编号进行梳理，并实时输出当前同化轮次压力分布曲线（包络线）；
%%          ⑦ 增加限定初始集合成员【必须满足收敛条件】才能进入后续同化过程，真实同化集合成员数量 N_new ≤ N_max；
%%          ⑧ 增加对区域 S-A模型修正项 Beta 设计扰动策略，初步测试收敛性并添加相应模块；
%%          ⑨ 更新可使用相对误差进行同化，相当于将同化目标转变为使相对误差趋向零（更新重要BUG）；
%%          ⑩ 增加Beta场可视化（同化单元数据→插值→同化节点数据）；
%%          ⑪ 增加便捷全同化程序续算开关；
%%          ⑫ 增加自适应观测点标准差调整模型（待测试）；
%%          ⑬ 增加对非翼型算例的前后处理支持；
%%          ⑭ 增加敏感性分析指标计算（2021/05/03）；
%%          ⑮ 增加边界层速度型同化功能（2021/05/18 初始 ）；
%% 重要更新：可选择使用集合卡尔曼滤波（EnKF）或集合变换卡尔曼滤波（ETKF）
%% 算例一：RAE2822翼型（耦合 CFD 程序） 
%% 标准状态：雷诺数 Re = 650W；马赫数 Ma = 0.729；迎角 AoA = 0.310deg
%% 算例二：S809翼型（耦合 CFD 程序）
%% 标准状态：雷诺数 Re = 200W；马赫数 Ma = 0.200；迎角 AoA = 10.2,11.2,12.2,14.0deg
%% 算例三：NACA4412翼型（耦合 CFD 程序）
%% 标准状态：雷诺数 Re = 152W；马赫数 Ma = 0.09（不可压状态）；迎角 AoA = 13.87deg

%% 程序框架实现
%% 设置同化程序全局续算开关
flag_da_xusuan=0; % 设置是否进行全同化程序续算
                  % 0 - 不开启续算；
                  % 1 - 开启中断续算，可导入任意中断后保存.mat文件后自由调整同化参数，开始进行续算；
                  % 2 - 直接定位至指定同化轮数 p 的保存文件，开始续算；
%% 初始化操作
% 若不开启续算，则清理所有工作区变量
if (flag_da_xusuan==0)
    clear all
    clc
    clf
    close all
    flag_da_xusuan=0;
    p=0; % 初始化同化轮数 p
elseif (flag_da_xusuan==2)
    clear all
    flag_da_xusuan=2;
    p=1; % 定位至指定同化轮数 p，继续同化
end

% 导入初始成员数据
% load('Pre_Cal_Path.mat');

format long;

% 定义全局变量
global savename new_output_folder fid_all fid_info PPNS initPath CFL_Setting Resi_target Resi_target_lev JX_NS STD_Comp JX_Turb L2NORM_Setting flag_Va_init;
global S_Va_turb_Name S_Va_turb_Type S_Va_turb_Num S_Va_turb_FilePath S_Va_turb_ModiLine S_Va_turb_ModiFormat S_Va_turb_Std;
global S_Va_target_Name S_Va_target_Num S_Va_target_TypeTurb S_Va_target_TypeDA S_Va_turb_target S_Va_Flow_length PT PS;
global S_Va_ref_Name S_Va_ref_Num;
global shape_upper shape_lower shape_outline shape_combined_arr_init Va_Beta Va_NNP_Beta Va_NNP_Beta_new;
global flag_da_relative flag_sigma_obs_adapt flag_sigma_obs_coeff D_obs_e_ori D_obs_e flag_flow_type;
global BL_setting N_BL L_BL BL_num BL_lab BL_dot BL_obs_num kb;
global Resi_flag Resi_flag_label;

%% 设置主要数据同化过程参数

% ------------  设置主要文件交互与并行相关参数 START ------------ %

% 设置数据交互路径
savename='Results_20210709_Case_RAE2822_Nmax_128_P_420_DA_81';
ddir_recal=('Program_CFD_Parallel_Recal');  % 设置重计算路径
ddir_post=('Program_CFD_Parallel_Post');    % 设置结果后处理路径
new_output_folder=['.\',savename,'\'];

% 若不开启续算开关，则直接建立新文件夹，保存计算结果
if (flag_da_xusuan==0)
    mkdir(new_output_folder);
    mkdir(['.\Output_Figures\',savename]);
end
% 若开启指定续算开关：直接定位至指定同化轮数 p 的保存文件，开始续算（注意：指定同化轮数必须已经计算完毕）；
if (flag_da_xusuan==2)
    load([new_output_folder,savename,'_S_Va_pre_restarat_',num2str(p,'%03d'),'.mat']);
end
% 重启/续算模块设置完毕

targetPath=new_output_folder;
initPath='Init_Ensemble_Member';

% 打开记录程序信息文件
fid_info=fopen([new_output_folder,'DA_CompInfo_Brief.dat'],'w+');
fid_all=fopen([new_output_folder,'DA_CompInfo_Detailed.dat'],'w+');

% 设置并行线程/分组数（可根据工作机线程设置）
PPNS=8;

% -------------  设置主要文件交互与并行相关参数 END ------------- %

% -------------  设置生成初始集合成员相关参数 START ------------- %

% 定义集合成员最大数量（注意：实际的 N_new ≤ N）
N_max=128;

% ----------  指定初始扰动变量、目标同化变量、同化观测变量模块 START ---------- %

% ① 指定初始扰动变量模块（目前默认为单变量扰动）
% 算例标准状态设置
% 定义初始扰动变量名称、文件路径、修改行数、修改格式（均用数组表示，一一对应）
% 初始扰动变量名称
S_Va_turb_Name=["Ma","AoA","CKS"];
% S_Va_turb_Name=["Beta"];
S_Va_turb_CharName=char(S_Va_turb_Name);
% 初始扰动变量类型：0 - 常数变量及其组合；1 - 空间分布变量
S_Va_turb_Type=0;
% 初始扰动变量个数
S_Va_turb_Num=length(S_Va_turb_Name);
% 初始扰动变量文件路径（省略运行程序根目录）
% S_Va_turb_FilePath_raw=char('\Input\InputNS.dat','\Input\InputNS.dat','\Input\InputTurb.dat');
% S_Va_turb_FilePath=strip(S_Va_turb_FilePath_raw,'right');
S_Va_turb_FilePath=["\Input\InputNS.dat","\Input\InputNS.dat","\Input\InputSAconst.dat"];
% 初始扰动变量修改行数（对应文件内）
S_Va_turb_ModiLine=[1,2,7];
% 初始扰动变量修改格式（对应文件内）
S_Va_turb_ModiFormat=["%-16.5f","%-16.5f","%-16.5f"];
% 初始扰动变量扰动区间
% 设置扰动方式：0 - 绝对值上下扰动；1 - 比例上下扰动（例如：设置50代表扰动区间为标准值的[50%,150%]）
S_Va_turb_way=[0,0,1];
S_Va_turb_Std=[0.729,2.310,0.410];  % 初始扰动变量标准值
S_Va_turb_Det_base_lower=[0.029,0.310,20];  % 初始扰动变量下区间
S_Va_turb_Det_base_upper=[0.029,0.310,20];  % 初始扰动变量上区间

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
% 利用拉丁超立方进行初始扰动变量采样
Va_min=S_Va_turb_min;
Va_max=S_Va_turb_max;
% 初始同化变量个数
PN=length(Va_min);
S_Va_coef=lhsdesign(N_max,PN);
% 生成N*P列扰动变量向量
S_Va_init=(ones(N_max,1)*Va_min)+S_Va_coef.*(ones(N_max,1)*(Va_max-Va_min));
S_Va_init_ave=mean(S_Va_init);

% ② 指定目标同化变量模块
% 全部目标同化变量名称
S_Va_target_Name=["Ma","AoA","TMu"];
% S_Va_target_Name=["Beta"];
S_Va_target_CharName=char(S_Va_target_Name);
% 全部目标同化变量数量
S_Va_target_Num=length(S_Va_target_Name);
% 判断目标同化变量类型①（TypeTurb）：k - 是第k个初始扰动变量或分布（k~=0）；0 - 非初始扰动变量
% 注：若目标同化变量非初始扰动变量，则需要在更新步将值改为原始 STD 值
S_Va_target_TypeTurb=[1,2,0];
% 判断目标同化变量类型②（TypeDA）：1 - 单（值）同化变量；2 - 流场（多值）同化变量
S_Va_target_TypeDA=[1,1,2];
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
% S_Va_ref_Name=["BL"];
S_Va_ref_CharName=char(S_Va_ref_Name);
% 全部同化观测变量数量
S_Va_ref_Num=length(S_Va_ref_Name);

% 若目标观测变量为Cp,Cf,yplus    
    % 设置参与同化过程的翼型数据点编号
    % num_da_min=1;    % 取样最小值
    % num_da_max=97;  % 取样最大值
    % num_da_gap=2;    % 取样间隔
    % 建立同化数据点序号数组
    % num_da_arr=(num_da_min:num_da_gap:num_da_max);
    % 可根据重点区域进行标志
    % num_da_arr=[(10:2:42),(48:64),(65:73),(74:2:92)];
    % num_da_arr=[(10:2:42),(48:92)];
    % num_da_arr=[(10:2:42),(48:73),(74:2:92)]; % RAE2822
    % num_da_arr=[(18:3:63),(90:2:168)]; % 56 Points Selected for S809 （P = 199）
    % num_da_arr=[(20:2:60),(94:2:166)]; % 58 Points Selected for S809 （P = 199）
    num_da_arr=[(12:7:75),(76:5:186),(208:5:288),(290:3:320),(321:5:376),(377:4:405)]; % 81 Points Selected for RAE2822 (P = 420)
    % 记录同化数据点个数
    k=length(num_da_arr);

% 若目标观测变量为边界层速度型，则需要指定边界层站位数量、导入各站位观测信息（y/c,u/u_ref）、各站位处边界层同化位置
    % N_BL=4;
    % BL_Type=[1,1,1,1];  % 设置各站位边界层取点模式：0 - 指定起始、终止编号以及间隔；1 - 指定同化点编号
    % BL_S=[]; % 设置各站位计算值起始编号
    % BL_T=[]; % 设置各站位计算值终止编号
    % BL_G=[]; % 设置各站位计算值编号间隔
    % global BL_setting N_BL L_BL BL_num BL_lab BL_dot
    BL_setting=textread('BL_DA_Setting.dat'); % 导入边界层同化设置文件
    N_BL=size(BL_setting,1);
    L_BL=size(BL_setting,2);
    % BL_dot=(-1)*ones(N_BL,100); % 保存边界层同化位置信息
    BL_num=zeros(1,N_BL); % 保存各站位边界层同化点数
    BL_lab=zeros(1,N_BL); % 保存参与同化的边界层编号
    for i=1:N_BL
        BL_lab(i)=BL_setting(i,1); % 保存边界层站位编号
        BL_type=BL_setting(i,2); % 设置各站位边界层取点模式：0 - 指定起始、终止编号以及间隔；1 - 指定同化点编号
        if (BL_type==0)
            BL_S=BL_setting(i,3);
            BL_T=BL_setting(i,4);
            BL_G=BL_setting(i,5);
            cnt_dot=0;
            for j=BL_S:BL_G:BL_T
                cnt_dot=cnt_dot+1;
                BL_dot(i,cnt_dot)=j;
            end
            BL_num(i)=cnt_dot;
        elseif (BL_type==1)
            cnt_dot=0;
            for kkk=3:L_BL
                if (BL_setting(i,kkk)~=0)
                    cnt_dot=cnt_dot+1;
                    BL_dot(i,cnt_dot)=BL_setting(i,kkk);
                end
            end
            BL_num(i)=cnt_dot;
        end
    end
    
    % Pos_Obs=(-1)*ones(N_BL,100); % 保存边界层同化位置信息
    % 导入指定边界层的速度型信息（实验数据 + 位置信息）
    kb=sum(BL_num); % 统计参与同化的边界层速度型点数量
    Exp_Obs_Arr=zeros(1,kb);
    BL_obs_num=0;
    for i=1:N_BL
        Vel_Exp_Data=load(['.\Input_Data\Input_Exp_Hump_BL_Vel_',num2str(BL_lab(i),'%03d'),'.dat']);
        Vel_Shape_Data=load(['.\Input_Data\Input_Hump_BL_Profile_',num2str(BL_lab(i),'%03d'),'.dat']);
        
        % 整理指定编号同化点的位置信息（y/c）
        for j=1:BL_num(i)
            Pos_Obs(i,j)=Vel_Shape_Data(BL_dot(i,j));
            % 插值得到边界层指定位置的速度实验数据
            Exp_Obs(i,j)=interp1(Vel_Exp_Data(:,1),Vel_Exp_Data(:,2),Pos_Obs(i,j),'linear','extrap');
            BL_obs_num=BL_obs_num+1; % 同化观测点 + 1
            Exp_Obs_Arr(BL_obs_num)=Exp_Obs(i,j); % 保存所有插值出的实验观测数据，数量 = k
        end
        
    end
    

% -----------  指定初始扰动变量、目标同化变量、同化观测变量模块 END ----------- %

% --------------  设置生成初始集合成员相关参数 END -------------- %

% -----------------  设置主要求解迭代参数 START ----------------- %

MStepNS_cov=12000;   % 定常NS求解收敛最大迭代步数（集合生成过程）
Max_inner_iter=10;  % 设置内迭代操作最大次数
MStepNS_pre=10000;   % 设置内迭代预测步数
MStepNS_recal_ensemble=18000;  % 设置每代成员重计算步数（集合重计算过程）
CFL_Setting=[5,5,5]; % 分别设置集合成员生成、内迭代以及重计算过程的 CFL 数
Resi_target_log=[-9.0,-9.0,-9.0]; % 分别定义集合成员生成、内迭代、重计算过程的收敛精度（以 10 为底）
Resi_target=10.^(Resi_target_log);
Resi_target_lev=[-4.0,-4.0,-4.0]; % 规定各阶段（集合生成、内迭代、重计算）的截断残差水平，以确定集合成员捕捉范围
% Resi_target_gen_level=Resi_target_log(1)+1.5; % 规定集合生成截断残差水平为期望残差值 + f ，以扩大集合成员捕捉范围
JX_NS=[0,1,0]; % 更新（关闭）JX_NS 续算操纵开关，即设置流场是否续算：0 - 否（从初始流场开始计算）；1 - 是 （从前一轮计算流场开始计算）

JX_Turb=[0,1,1]; % 更新（关闭）JX_Turb 续算操纵开关，即设置涡粘场是否续算：0 - 否（不固定涡粘，直接采用湍流模型计算）；1 - 是 （固定涡粘进行计算）；3 - 引入 Beta 项进行计算
STD_Comp=[1,1,1]; % 是否进行当前算例与标准算例计算结果的比较（流场各变量 RU、U、V、P、TMu 的绝对误差场，以及总MSE、RMSE）
                  % 需要在Output文件夹中导入 FlowField_STD.PLT 文件，输出 FlowField_Err.PLT 和 FlowField_Err_Info.txt 文件
L2NORM_Setting=[0,0,0]; % 设置每轮 L2NORM 收敛要求水平，达到即退出当前计算（注：5 - 代表收敛截断水平为 5%；当且仅当湍流续算开关 JX_Turb = 1 时生效）


% 设置算例类型（0 - 翼型绕流：需导入上下翼型外形数据以及实验数据，共4个输入文件，并进行插值排序）
%               1 - 单壁内流，如后台阶（Backstep）、驼峰（Hump）等：需导入壁面外形参数以及实验数据，共2个输入文件，并进行插值排序
flag_flow_type=0;

flag_std=0;  % 设置是否进行初始标准算例：0 - 否；1 - 是
flag_da_part=0;  % 设置是否进行部分区域数据同化：0 - 否；1 - 是
flag_da_method=2; % 设置数据同化核心算法：1 - 集合卡尔曼滤波（EnKF）；2 - 集合变换卡尔曼滤波（ETKF）
flag_da_relative=1; % 设置数据同化观测值类型：0 - 相对误差，同化目标变为使相对误差变为零；1 - 绝对误差（传统方式）
flag_beta_turb_type=0; % 设置网格取样点规则（对于 Beta修正项 同化情形）：0 - 每次扰动点固定；1 - 每次扰动点随机；

% 如果初始扰动与同化目标中出现 "Beta" 项，则进行扰动设置
if ((strcmp(char(S_Va_target_Name(1)),'Beta')==1) || (strcmp(char(S_Va_target_Name(1)),'Beta')==1))
    flag_da_part=1; % 开启区域数据同化
    flag_da_beta=1; % 开启扰动项Beta(x)同化
    flag_Va_init=1; % 开启集合成员导入流场多维工作变量标志
    flag_sel_dot=2; % 设置扰动Beta(x)的方式
                    % 拟定规则 1 - 直接选取相关编号间隔，即在选取区域内每隔一定编号进行单点采样，选取散点并进行平面插值；
                    % 拟定规则 2 - 精细化采样均匀采样；
                    % 拟定规则 3 - 采样二维光滑函数进行采样（待开发......）；
    % 设置 Beta 扰动规则：方差 sigma_2，满足分布 Beta ~ N(std_beta,sigma_2)
    % 规则①取点间隔信息
    DA_NG_gap=200;
    % 规则①~③取点距离信息
    std_beta=1.0;
    sigma_beta=3.0;
    dis_ratio=5.0;
else
    % 若不进行 Beta扰动项 同化
    flag_da_beta=0; % 关闭扰动项Beta(x)同化
    flag_Va_init=0; % 关闭集合成员导入流场多维工作变量，使用普通单变量组合扰动
end

% 区域数据同化设置
NNDA=4; % NNDA -- 同化区域节点个数，可为任意多边形（须满足>=3），下面输入 NNDA 个坐标（XDA,YDA）
% XDA=[1.05,1.2,1.2,1.05]; % XDA(NNDA) -- 同化区域节点横坐标
% YDA=[0.1,0.1,0.0,0.0]; % YDA(NNDA) -- 同化区域节点纵坐标
% XDA=[0.8,1.5,1.5,0.8]; % XDA(NNDA) -- 同化区域节点横坐标
% YDA=[0.1,0.1,0.02,0.02]; % YDA(NNDA) -- 同化区域节点纵坐标
XDA=[0.5499846528026755,1.3977168074273609,1.3977168074273609,0.5499846528026755]; % XDA(NNDA) -- 同化区域节点横坐标
YDA=[0.11283900866252453,0.16714063248086963,0.005,0.005]; % YDA(NNDA) -- 同化区域节点纵坐标

if (flag_flow_type==0)
    % 导入上下翼面的外形数据（自动插值）
    shape_upper=load('.\Input_Data\Input_Data_Shape_RAE2822_upper.dat');
    shape_lower=load('.\Input_Data\Input_Data_Shape_RAE2822_lower.dat');
    % 导入上下翼面的实验数据
    exp_upper_data=load('.\Input_Data\Input_Data_Exp_RAE2822_upper.dat');
    exp_lower_data=load('.\Input_Data\Input_Data_Exp_RAE2822_lower.dat');
elseif (flag_flow_type==1)
    % 导入外形数据（自动插值）
    shape_outline=load('.\Input_Data\Input_Data_Shape_Hump.dat');
    % 导入实验数据
    exp_outline=load('.\Input_Data\Input_Data_Exp_Hump.dat');
end
    
% ------------------  设置主要求解迭代参数 END ------------------ %

% -----------  设置数据同化实验值与误差估计参数 START ----------- %

% 设置观测值协方差（方差），并假定这些测量误差不相关
% !!!!!!! 方差具体设置值需要再考虑（大一些？自适应方差调整模型？） !!!!!!!
% 2021/04/27 更新模块
% 设置是否开启 自适应观测值标准差调节模型 - Adaptive observation standard deviation adjustment model 
%               abbr. AOSDAM（否 - 0；固定放缩模式 - 1；区间平均模式 - 2；区间适应模式 - 3）
% 算法说明：若关闭自适应模型，则设定所有观测值标准差为一固定值 sigma_obs；
%           若开启自适应模型（1 - 固定放缩模式）：则每轮同化前所有观测值标准差设置为上一轮标准差乘以放缩因子；
%           若开启自适应模型（2 - 区间平均模式），则每轮同化前根据上一轮所有观测值对应的计算值包络区间平均宽度，将所有观测值的标准差调整为平均期望标准差水平（乘以放缩因子）；
%           若开启自适应模型（3 - 区间适应模式），则每轮同化前根据上一轮各观测值对应的计算值包络区间宽度，将各观测值的标准差调整为相应期望标准差水平（乘以放缩因子）；
% 变量说明：flag_sigma_obs_adapt - 模式开关
%           flag_sigma_obs_coeff - 放缩因子
%           sigma_obs_init       - 设定初始观测值标准差水平
%           sigma_obs_min        - 设定观测值标准差下限，即调整前后的标准差均不可低于下限
%           sigma_obs_max        - 设定观测值标准差上限，即调整前后的标准差均不可高于上限
flag_sigma_obs_adapt=0;
flag_sigma_obs_coeff=0.80;
sigma_obs_ini=0.03; % 设置初始观测标准差水平
sigma_obs_min=0.02; % 指定观测标准差下限为0.05，对应方差为0.0025
sigma_obs_max=0.80; % 指定观测标准差上限为0.80，对应方差为0.6400

% 计算数据同化初始集合分析值误差协方差
% C_pp_e_a=cov(S_Va_ana');

% ----------  设置数据同化实验值与误差估计参数 END ---------- %

% ---------------  设置部分数据同化区域 START --------------- %

% 区域数据模块同化设置
if (flag_da_part==1)
    
    % 将设置参数写入 InputDAArea.dat 文件
    % 指定替换行，分别对应输入文件第1、2、3行
    rep_line_NNDA=1;
    rep_line_XDA=2;
    rep_line_YDA=3;
    rep_format_NNDA='%-6d\n';
    rep_format_XDA_init='%-6.3f';
    rep_format_YDA_init='%-6.3f';
    for i=1:(NNDA-1)
        rep_format_XDA=[rep_format_XDA_init,' %-6.3f'];
        rep_format_YDA=[rep_format_YDA_init,' %-6.3f'];
        rep_format_XDA_init=rep_format_XDA;
        rep_format_YDA_init=rep_format_YDA;
    end
    rep_format_XDA=[rep_format_XDA_init,'\n'];
    rep_format_YDA=[rep_format_YDA_init,'\n'];

    fid_DA_Path='.\Program_CFD_Base\Input\InputDAArea.dat';
    fid_DA=fopen(fid_DA_Path,'w+');
    fprintf(fid_DA,rep_format_NNDA,NNDA);
    fprintf(fid_DA,rep_format_XDA,XDA);
    fprintf(fid_DA,rep_format_YDA,YDA);
    fclose(fid_DA);
    
    % 若开启区域数据同化，则设 Input 文件中的对应行中 JX_DA = 1
    % 指定替换行，对应输入文件（InputNS.dat）第21行
    rep_line_JX_DA=21;
    rep_format_JX_DA='%-12d';
    % 更新（关闭）操纵区域同化开关
    JX_DA=2;
    fid_DA=fopen('.\Program_CFD_Base\Input\InputNS.dat','r+');
    for i=1:(rep_line_JX_DA-1)
       fgetl(fid_DA);
    end
    fseek(fid_DA,0,'cof');
    fprintf(fid_DA,rep_format_JX_DA,JX_DA);
    fclose(fid_DA);
    
    % 若进行部分区域数据同化，则绘制（非）翼型外形与指定同化区域示意图
    if (flag_flow_type==0)
        DA_Extension_Airfoil_Interpolation(NNDA,XDA,YDA);
    elseif (flag_flow_type==1)
        DA_Extension_NonAirfoil_Interpolation(NNDA,XDA,YDA);
    end
    
end

% ----------------  设置部分数据同化区域 END ---------------- %

% ---------------  输出算例主要设置信息 START --------------- %

% 命令行输出提示
fprintf('******************************************************************\n');
fprintf('                      数据同化算例主要设置信息                    \n');
fprintf('******************************************************************\n');
fprintf('[算例类型设置]\n');
fprintf('设置算例所属类型（0 - 翼型绕流；1 - 固壁内流） - %d\n',flag_flow_type);
fprintf('[程序保存路径 & 并行设置]\n');
fprintf(['结果保存路径：',savename,'\n']);
fprintf(['程序信息文件：DA_CompInfo_Brief.dat\n']);
fprintf('程序并行线程/分组数 PPNS = %d\n',PPNS);
fprintf('******************************************************************\n');
fprintf('[同化主要参数]\n');
fprintf('最大集合成员数量 N_max = %d\n',N_max);
fprintf('同化观测点数 k = %d, 具体分布请查看 num_da_arr 变量\n',k);
fprintf('******************************************************************\n');
fprintf('[初始扰动变量信息]\n');
if (flag_Va_init==0)
    for i=1:S_Va_turb_Num
        fprintf(['初始扰动变量 ',num2str(i,'%03d'),'：',S_Va_turb_CharName(:,:,i),' ∈ [ %.5f , %.5f ], STD = %.5f\n'],S_Va_turb_min(i),S_Va_turb_max(i),S_Va_turb_Std(i));
    end
else
    for i=1:S_Va_turb_Num
        fprintf(['初始扰动分布 ',num2str(i,'%03d'),'：',S_Va_turb_CharName(:,:,i),' \n']);
    end
end
fprintf('******************************************************************\n');
fprintf('[目标同化变量信息]\n');
for i=1:S_Va_target_Num
    fprintf(['目标同化变量 ',num2str(i,'%03d'),'：',S_Va_target_CharName(:,:,i),' 初始扰动变量 - %d, 单/多（值）同化变量 - %d\n'],S_Va_target_TypeTurb(i),S_Va_target_TypeDA(i));
end
fprintf('******************************************************************\n');
fprintf('[指定同化观测变量信息]\n');
for i=1:S_Va_ref_Num
    fprintf(['指定同化观测变量 ',num2str(i,'%03d'),'：',S_Va_ref_CharName(:,:,i),'\n']);
end
fprintf('******************************************************************\n');
fprintf('[同化流程控制]\n');
fprintf('设置集合生成过程最大迭代步数：%d, CFL = %d\n',MStepNS_cov,CFL_Setting(1));
fprintf('设置集合内迭代预测计算步数：%d, CFL = %d\n',MStepNS_pre,CFL_Setting(2));
fprintf('设置集合成员重计算步数：%d, CFL = %d\n',MStepNS_recal_ensemble,CFL_Setting(3));
fprintf('设置集合内迭代操作最大迭代轮次：%d\n',Max_inner_iter);
fprintf('是否进行初始标准算例 - %d\n',flag_std);
fprintf('是否进行区域数据同化 - %d\n',flag_da_part);
fprintf('设置数据同化核心算法（1 - EnKF；2 - ETKF） - %d\n',flag_da_method);
fprintf('设置数据同化观测值类型（0 - 相对误差；1 - 绝对误差） - %d\n',flag_da_relative);
fprintf('是否进行湍流源项修正项 Beta 同化 - %d\n',flag_da_beta);
fprintf('是否在集合生成过程导入流场多维工作变量 - %d\n',flag_Va_init);
fprintf('是否开启 自适应观测值标准差调节模型（AOSDAM） 选择模式 - %d\n',flag_sigma_obs_adapt);
fprintf('******************************************************************\n');
fprintf('\n');

% 文件输出提示
fprintf(fid_info,'******************************************************************\n');
fprintf(fid_info,'                      数据同化算例主要设置信息                    \n');
fprintf(fid_info,'******************************************************************\n');
fprintf(fid_info,'[算例类型设置]\n');
fprintf(fid_info,'设置算例所属类型（0 - 翼型绕流；1 - 固壁内流） - %d\n',flag_flow_type);
fprintf(fid_info,'[程序保存路径 & 并行设置]\n');
fprintf(fid_info,['结果保存路径：',savename,'\n']);
fprintf(fid_info,['程序信息文件：DA_CompInfo_Brief.dat\n']);
fprintf(fid_info,'程序并行线程/分组数 PPNS = %d\n',PPNS);
fprintf(fid_info,'******************************************************************\n');
fprintf(fid_info,'[同化主要参数]\n');
fprintf(fid_info,'最大集合成员数量 N_max = %d\n',N_max);
fprintf(fid_info,'同化观测点数 k = %d, 具体分布请查看 num_da_arr 变量\n',k);
fprintf(fid_info,'******************************************************************\n');
fprintf(fid_info,'[初始扰动变量信息]\n');
if (flag_Va_init==0)
    for i=1:S_Va_turb_Num
        fprintf(fid_info,['初始扰动变量 ',num2str(i,'%03d'),'：',S_Va_turb_CharName(:,:,i),' ∈ [ %.5f , %.5f ], STD = %.5f\n'],S_Va_turb_min(i),S_Va_turb_max(i),S_Va_turb_Std(i));
    end
else
    for i=1:S_Va_turb_Num
        fprintf(fid_info,['初始扰动分布 ',num2str(i,'%03d'),'：',S_Va_turb_CharName(:,:,i),' \n']);
    end
end
fprintf(fid_info,'******************************************************************\n');
fprintf(fid_info,'[目标同化变量信息]\n');
for i=1:S_Va_target_Num
    fprintf(fid_info,['目标同化变量 ',num2str(i,'%03d'),'：',S_Va_target_CharName(:,:,i),' 初始扰动变量 - %d, 单/多（值）同化变量 - %d\n'],S_Va_target_TypeTurb(i),S_Va_target_TypeDA(i));
end
fprintf(fid_info,'******************************************************************\n');
fprintf(fid_info,'[指定同化观测变量信息]\n');
for i=1:S_Va_ref_Num
    fprintf(fid_info,['指定同化观测变量 ',num2str(i,'%03d'),'：',S_Va_ref_CharName(:,:,i),'\n']);
end
fprintf(fid_info,'******************************************************************\n');
fprintf(fid_info,'[同化流程控制]\n');
fprintf(fid_info,'设置集合生成过程最大迭代步数：%d, CFL = %d\n',MStepNS_cov,CFL_Setting(1));
fprintf(fid_info,'设置集合内迭代预测计算步数：%d, CFL = %d\n',MStepNS_pre,CFL_Setting(2));
fprintf(fid_info,'设置集合成员重计算步数：%d, CFL = %d\n',MStepNS_recal_ensemble,CFL_Setting(3));
fprintf(fid_info,'设置集合内迭代操作最大迭代轮次：%d\n',Max_inner_iter);
fprintf(fid_info,'是否进行初始标准算例 - %d\n',flag_std);
fprintf(fid_info,'是否进行区域数据同化 - %d\n',flag_da_part);
fprintf(fid_info,'设置数据同化核心算法（1 - EnKF；2 - ETKF） - %d\n',flag_da_method);
fprintf(fid_info,'设置数据同化观测值类型（0 - 相对误差；1 - 绝对误差） - %d\n',flag_da_relative);
fprintf(fid_info,'是否进行湍流源项修正项 Beta 同化 - %d\n',flag_da_beta);
fprintf(fid_info,'是否在集合生成过程导入流场多维工作变量 - %d\n',flag_Va_init);
fprintf(fid_info,'是否开启 自适应观测值标准差调节模型（AOSDAM） 选择模式 - %d\n',flag_sigma_obs_adapt);
fprintf(fid_info,'******************************************************************\n');
fprintf(fid_info,'\n');

% ----------------  输出算例主要设置信息 END ---------------- %

% -----------------  首轮执行标准算例 START ----------------- %

% 运行一次标准算例（Program_CFD_Base）
if (flag_std==1)
    fprintf('标准算例程序组 开始运行\n');
    fprintf(fid_info,'标准算例程序组 开始运行\n');
    delete_std_Path='Program_CFD_Base';
    if (exist([delete_std_Path,'\ProgramExitFlag.txt'],'file'))~=0
        delete([delete_std_Path,'\ProgramExitFlag.txt']);
    end
    cmd_std=('.\Program_CFD_Base\CFD_2D.exe');
    open(cmd_std);
    pause(3);
    exit_std_flag=0;
    while (exit_std_flag~=1)
        if (exist([delete_std_Path,'\ProgramExitFlag.txt'],'file'))~=0
            exit_std_flag=1;
            delete([delete_std_Path,'\ProgramExitFlag.txt']);
            break;
        end
        pause(10);
    end
    fprintf('标准算例程序组 运行结束\n');
    fprintf(fid_info,'标准算例程序组 运行结束\n');
    fprintf('\n');
    fprintf(fid_info,'\n');
else
    fprintf('已存在标准算例程序组\n');
    fprintf(fid_info,'已存在标准算例程序组\n');
    fprintf('\n');
    fprintf(fid_info,'\n');
end

% ------------------  首轮执行标准算例 END ------------------ %

%% 步骤一：导入参考实验与首轮计算数据，并进行插值操作
% 计算数据操作
% 导入翼型原始外形数据，并进行上下翼型插值
% 导入计算数据
com_data=load('.\Program_CFD_Base\Output\CP_CF_Comp.PLT');
% 指定调整后计算数据保存路径
fid_com_adj_path='.\Program_CFD_Base\Output\CP_CF_Comp_New.PLT';
fid_tar_new_path='.\Program_CFD_Base\Output\CP_New.PLT';

% 翼型算例实验值处理
if (flag_flow_type==0)
    
    % 对计算压力系数曲线节点编号进行梳理，并实时输出当前同化轮次压力分布曲线
    % 集合生成过程、内迭代过程中全部成员完成计算后保存曲线包络线
    [shape_upper_num,shape_lower_num,shape_upper_arr_sor,shape_lower_arr_sor,n_com]=DA_Va_Airfoil_Ref_Adjustment_and_Visulization(com_data,fid_com_adj_path,fid_tar_new_path);

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
    
elseif (flag_flow_type==1)
    
    % 对计算压力系数曲线节点编号进行梳理，并实时输出当前同化轮次压力分布曲线
    % 集合生成过程、内迭代过程中全部成员完成计算后保存曲线包络线
    [shape_num,shape_arr_sor,n_com]=DA_Va_NonAirfoil_Ref_Adjustment_and_Visulization(com_data,fid_com_adj_path,fid_tar_new_path);
    % 实验数据操作
    % 将实验数据数组根据x坐标进行升序排序（指定[0,1]区间）
    exp_data_sor=sortrows(exp_outline,1);
    % 写入输出文件夹下的 InputExpData_all.dat 与 InputExpData.PLT 文件
    exp_data_out=sortrows(exp_outline,1,'ascend');
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
    exp_num_tot=length(exp_data_sor);
    % 进行上下翼面实验数据值的插值
    xx_lower_exp=(0:0.001:1);
    yy_lower_exp=interp1(exp_data_sor(:,1),exp_data_sor(:,2),xx_lower_exp,'linear','extrap');
    % 利用插值方法求出计算数据点坐标对应的实验值
    for i=1:shape_num
        x_com_pre=shape_arr_sor(i,2);
        shape_arr_sor(i,7)=interp1(exp_data_sor(:,1),exp_data_sor(:,2),x_com_pre,'linear','extrap');
    end
    
    shape_combined_arr_init=shape_arr_sor;

end

% 建立由实验值至计算点的映射关系
% 记录参与同化的实验数据点信息，可导出文件
exp_for_DA=zeros(k,2);
for i=1:k

    exp_for_DA(i,1)=shape_combined_arr_init(num_da_arr(i),2);

    if (strcmp(char(S_Va_ref_Name(1)),'Cp')==1)
        exp_for_DA(i,2)=shape_combined_arr_init(num_da_arr(i),7);
    elseif (strcmp(char(S_Va_ref_Name(1)),'Cf')==1)
        exp_for_DA(i,2)=shape_combined_arr_init(num_da_arr(i),7);
    elseif (strcmp(char(S_Va_ref_Name(1)),'yplus')==1)
        exp_for_DA(i,2)=shape_combined_arr_init(num_da_arr(i),7);
    end

end

% 输出进行同化的观测实验值分布
fid_exp_da=fopen('.\Program_CFD_Base\Output\InputExpData_for_DA.PLT','w+');
fprintf(fid_exp_da,'VARIABLES="X","CP"\n');
fprintf(fid_exp_da,'ZONE T="ZONE EXP for DA"\n');
for i=1:k
    fprintf(fid_exp_da,'%20.8f %20.8f\n',exp_for_DA(i,1),exp_for_DA(i,2));
end
fclose(fid_exp_da);

% 复制标准湍流模型计算结果，方便模板化显示
copyfile('.\Program_CFD_Base\Output\CP_New.PLT','.\Program_CFD_Base\Output\CP_SA.PLT');
% 复制标准湍流模型计算结果，方便结果比较（进行误差及敏感性分析）
copyfile('.\Program_CFD_Base\Output\FlowField.PLT','.\Program_CFD_Base\Output\FlowField_Ori.PLT');
copyfile('.\Program_CFD_Base\Output\FlowField_Info.PLT','.\Program_CFD_Base\Output\FlowField_STD.PLT');

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
saveas(gcf,['.\Output_Figures\',savename,'_Init_State_Output.fig']);

%% 步骤二：并行程序组创建

% 提取翼型网格参数
% 运行网格转换主程序
fprintf('翼型绕流数据同化程序启动！\n');
fprintf(fid_info,'翼型绕流数据同化程序启动！\n');
fprintf('\n');
fprintf(fid_info,'\n');
cmd_grid=('.\Program_CFD_Base\GridTreat_.exe');
open(cmd_grid);
pause(3);
fprintf('基础程序网格转换完成！\n');
fprintf(fid_info,'基础程序网格转换完成！\n');

% 创建程序运行文件夹（数量为PPNS）
for po=1:PPNS
    ddir=['Program_CFD_Parallel_',num2str(po,'%03d')];
    copyfile('Program_CFD_Base',ddir);
    % 首次运行网格转换
    cmd_grid=([ddir,'\GridTreat_.exe']);
    open(cmd_grid);
    pause(3);
    fprintf('并行程序组 %03d 网格转换完成！\n',po);
    fprintf(fid_info,'并行程序组 %03d 网格转换完成！\n',po);
end

% 创建重计算过程程序组
copyfile('Program_CFD_Base',ddir_recal);
cmd_grid=([ddir_recal,'\GridTreat_.exe']);
open(cmd_grid);
pause(3);
fprintf('重计算程序组 网格转换完成！\n');
fprintf(fid_info,'重计算程序组 网格转换完成！\n');

% 创建结果后处理程序组
copyfile('Program_CFD_Base',ddir_post);
cmd_grid=([ddir_post,'\GridTreat_.exe']);
open(cmd_grid);
pause(3);
fprintf('结果后处理程序组 网格转换完成！\n');
fprintf(fid_info,'结果后处理程序组 网格转换完成！\n');

fprintf('\n');
fprintf(fid_info,'\n');

%% 步骤三：初始化集合成员

% ---------------- 敏感区域数据修正系数生成 START ---------------- %

% 若进行湍流模型源项修正项 Beta 同化
if (flag_da_beta==1)
    
    % 提取区域数据同化信息
    DA_area_info=load('.\Program_CFD_Base\Input\InputDAgrid.dat');
    % 第一行表示同化区域内网格节点
    NDAV=DA_area_info(1);
    % 第二行表示同化区域内单元数量
    NDAG=DA_area_info(2);
    % 提取同化区域内网格节点与单元编号信息
    DA_VN=DA_area_info(3:(NDAV+2));
    DA_GN=DA_area_info((NDAV+3):(NDAV+NDAG+2));
    
    % 提取当前算例所有网格节点、单元格心坐标信息
    VN_grid=load('.\Program_CFD_Base\Input\Input_NNP_Coordinate.dat');
    GN_grid=load('.\Program_CFD_Base\Input\Input_NG_Coordinate.dat');
    % 提取所有网格节点坐标
    VN_grid_x=VN_grid(:,1);
    VN_grid_y=VN_grid(:,2);
    NNP_tot=length(VN_grid);
    % 提取所有网格单元格心坐标
    GN_grid_x=GN_grid(:,1);
    GN_grid_y=GN_grid(:,2);
    NG_tot=length(GN_grid);
    
    % 提取选取区域内网格节点、单元格心坐标信息
    DA_area_coor_info=load('.\Program_CFD_Base\Input\InputDAgrid_Coordinate.dat');
    % 提取区域内网格节点个数
    NDAV_coor=DA_area_coor_info(1,1);
    % 提取区域内网格单元个数
    NDAG_coor=DA_area_coor_info(1,2);
    % 提取选取区域内网格节点与单元编号信息
    DA_VN_x=DA_area_coor_info((2:(NDAV_coor+1)),1); % 选取区域内网格节点 x 坐标
    DA_VN_y=DA_area_coor_info((2:(NDAV_coor+1)),2); % 选取区域内网格节点 y 坐标
    DA_GN_x=DA_area_coor_info(((NDAV_coor+2):(NDAV_coor+NDAG_coor+1)),1); % 选取区域内网格单元 x 坐标
    DA_GN_y=DA_area_coor_info(((NDAV_coor+2):(NDAV_coor+NDAG_coor+1)),2); % 选取区域内网格单元 x 坐标
    
    % 标记出现在同化区域内的网格单元编号，是 - 1；否 - 0
    DA_GN_flag_tot=zeros(NG_tot,1);
    for i=1:NDAG
        DA_GN_flag_tot(DA_GN(i))=1;
    end
    
    % 初始化 修正项Beta 矩阵
    Va_Beta=zeros(NG_tot,N_max);
    Va_NNP_Beta=zeros(NNP_tot,N_max);
    
    % 设置网格取样点规则（每次扰动点随机 or 固定？）
    %     if (flag_sel_dot==1)
    %         % 拟定规则①：直接选取相关编号间隔，即在选取区域内每隔一定编号进行单点采样，选取散点
    %         DA_NG_gap=200;
    %     elseif (flag_sel_dot==2)
    %         % 拟定规则②：圆形区域？
    %         XY_dots_num_max=round(0.1*NDAG_coor);
    %         XY_dots_dis_lim=max(((max(XDA)-min(XDA))/XY_dots_num_max),((max(YDA)-min(YDA))/XY_dots_num_max));
    %         [XY_cnt,Xp_sca,Yp_sca,XY_sel_flag]=Region_Beta_Dots_Selected(NDAG_coor,DA_GN_x,DA_GN_y,XY_dots_num_max,XY_dots_dis_lim);
    %     end
    %     figure;
    %     scatter(Xp_sca,Yp_sca);
    %     xlabel('x'),ylabel('y');
    %     XY_sel_flag_real=zeros(XY_cnt,1);
    %     for i=1:XY_cnt
    %         XY_sel_flag_real(i)=DA_GN(XY_sel_flag(i));
    %     end
    
    % 提取边界层第一层网格单元数，并释放临时内存
    % 目的：将选定区域之外，以及边界层第一层网格单元的 Beta扰动值
    % 设置为标准值 1.0，且（如果有）将边界层第一层网格单元排除于选定区域之外
    Grid_info=textread('.\Program_CFD_Base\Input\US2D.dat');
    BD_num=Grid_info(1,1); % 提取边界层第一层的网格单元数量
    clear Grid_info; % 清除网格信息
    
    for pp=1:N_max

        % 根据拟定规则选取区域散点，并进行扰动Beta值的赋值
        XY_NG_Beta=ones(NG_tot,1); % 初始化 Beta 矩阵
        XY_NNP_Beta=ones(NNP_tot,1); % 初始化 Beta 矩阵
        
        if (flag_sel_dot==1)
            % 拟定规则①：直接选取相关编号间隔，即在选取区域内每隔一定编号进行单点采样，选取散点
            % DA_NG_gap=200;
        
            XY_dots_num=0;
            Xp_sca=zeros(length(1:DA_NG_gap:NDAG_coor),1);
            Yp_sca=zeros(length(1:DA_NG_gap:NDAG_coor),1);
            Vp_sca=zeros(length(1:DA_NG_gap:NDAG_coor),1);
            DA_NG_note=zeros(NDAG_coor,1);
            for i=1:DA_NG_gap:NDAG_coor
                XY_dots_num=XY_dots_num+1;
                Xp_sca(XY_dots_num)=DA_GN_x(i);
                Yp_sca(XY_dots_num)=DA_GN_y(i);
                Vp_sca(XY_dots_num)=normrnd(std_beta,sigma_beta);
                XY_NG_Beta(DA_GN(i))=Vp_sca(XY_dots_num);
                DA_NG_note(i)=1;
            end
            % 利用MATLAB自带的 scatteredInterpolant() 函数对选取的散点进行曲面插值，构建 V = F(x,y) 的拟合曲面
            F=scatteredInterpolant(Xp_sca,Yp_sca,Vp_sca);
            F.Method='natural';
            F.ExtrapolationMethod='nearest';
            Vp_NG_fit=F(DA_GN_x,DA_GN_y);
            Vp_NNP_fit=F(DA_VN_x,DA_VN_y);
            % figure;
            % hold on
            for i=1:NDAG_coor
                if (DA_NG_note(i)==0)
                    Xp_pre=DA_GN_x(i);
                    Yp_pre=DA_GN_y(i);
                    Vp_pre=F(Xp_pre,Yp_pre);
                    XY_NG_Beta(DA_GN(i))=Vp_pre;
                    % plot3(Xp_pre,Yp_pre,Vp_pre,'r.');
                end
            end
            % hold off
            % Va_Beta(:,pp)=XY_Beta;

        elseif (flag_sel_dot==2)
            
            % 拟定规则②：精细化采样 + 扰动
            % 设置扰动限制条件
            XY_dots_num_max=round(0.1*NDAG_coor);
            XY_dots_dis_lim=dis_ratio*max(((max(XDA)-min(XDA))/XY_dots_num_max),((max(YDA)-min(YDA))/XY_dots_num_max));
            
            % 若 flag_beta_turb_type==1，则打开随机扰动网格单元开关，则每次均会随机产生 Beta初始扰动点，均会更新 XY_cnt 与 XY_sel_flag 值
            % 若 flag_beta_turb_type==0，则关闭随机扰动网格单元开关，则仅产生一轮 Beta初始扰动点，之后固定
            if ((flag_beta_turb_type==1) || ((flag_beta_turb_type==0) && (pp==1)))
                [XY_cnt_pre,Xp_sca,Yp_sca,XY_sel_flag_pre]=Region_Beta_Dots_Selected(NDAG_coor,DA_GN_x,DA_GN_y,XY_dots_num_max,XY_dots_dis_lim);
            end
            
            XY_sel_flag_real=zeros(XY_cnt_pre,1);
            for i=1:XY_cnt_pre
                XY_sel_flag_real(i)=DA_GN(XY_sel_flag_pre(i));
            end
            
            % 区域 Beta 扰动值初始化
            % 将全部 Beta 值设置为标准值 1.0（固定边界值）
            Vp_NG_std_all=ones(NG_tot,1);
            % 将指定网格单元进行 Beta 扰动
            for i=1:XY_cnt_pre
                Vp_NG_std_all(XY_sel_flag_real(i))=normrnd(std_beta,sigma_beta);
            end
            % 将边界层第一层单元Beta值设置为标准值 1.0 并覆盖
            for i=1:BD_num
                Vp_NG_std_all(i)=1.0;
            end
            % 提取非平面插值部分信息（除抽取的扰动点）
            cnt_non_fit=0;
            for i=1:NG_tot
                if ((i<=BD_num) || (DA_GN_flag_tot(i)==0) || (sum(find(XY_sel_flag_real==i))~=0))
                    cnt_non_fit=cnt_non_fit+1;
                    Xp_non_all(cnt_non_fit)=GN_grid_x(i);
                    Yp_non_all(cnt_non_fit)=GN_grid_y(i);
                    Vp_non_all(cnt_non_fit)=Vp_NG_std_all(i);
                end
            end
            
            % 利用MATLAB自带的 scatteredInterpolant() 函数对选取的散点进行曲面插值，构建 V = F(x,y) 的拟合曲面
            F=scatteredInterpolant(Xp_non_all',Yp_non_all',Vp_non_all');
            F.Method='natural';
            F.ExtrapolationMethod='nearest';
            
            % 插值得到所有网格单元 Beta 值
            Vp_NG_fit=F(DA_GN_x,DA_GN_y); % 指定区域单元插值
            Xp_NG_all=GN_grid_x;
            Yp_NG_all=GN_grid_y;
            XY_NG_Beta=F(Xp_NG_all,Yp_NG_all); % 指定区域单元插值
            
            % 利用已有网格单元 Beta 值对指定区域内的单元节点进行插值，并组合得到全场网格节点的 Beta 值，然后实现可视化
            % 将全部 Beta 值（全场网格节点）设置为标准值 1.0
            % F=scatteredInterpolant(Xp_non_all',Yp_non_all',XY_NG_Beta');
            Vp_NNP_std_all=ones(NNP_tot,1);
            Vp_NNP_fit=F(DA_VN_x,DA_VN_y); % 指定区域网格节点插值
            Xp_NNP_all=VN_grid_x;
            Yp_NNP_all=VN_grid_y;
            XY_NNP_Beta=F(Xp_NNP_all,Yp_NNP_all); % 插值得到所有网格节点 Beta 值，利用网格单元值插值？
            
        end
        
        % 完成初始 Beta 分布变量采样生成
        Va_Beta(:,pp)=XY_NG_Beta;
        Va_NNP_Beta(:,pp)=XY_NNP_Beta;
        
        % 清理内存
        clear Xp_non_all;
        clear Yp_non_all;
        clear Vp_non_all;
    
    end
    
    % 框选扰动区域可视化（基于网格单元 Beta 数据）
    figure;    
    scatter(Xp_sca,Yp_sca);
    xlabel('x'),ylabel('y');
    % 保存可编辑图片
    saveas(gcf,['.\Output_Figures\',savename,'\Init_Final_Beta_SampleScatter.fig']);
    
    figure;
    T=delaunay(DA_GN_x,DA_GN_y);
    TO=triangulation(T,DA_GN_x(:),DA_GN_y(:),Vp_NG_fit(:));
    % T=delaunay(DA_VN_x,DA_VN_y);
    % TO=triangulation(T,DA_VN_x(:),DA_VN_y(:),Vp_NNP_fit(:));
    % trimesh(TO)
    trisurf(TO)
    shading interp
    xlabel('x'),ylabel('y'),zlabel('\beta');
    % 保存可编辑图片
    saveas(gcf,['.\Output_Figures\',savename,'\Init_Final_Beta_SampleContour.fig']);
    
end

% ---------------- 敏感区域数据修正系数生成 END ---------------- %

% -------------- 修改参数预置 + 集合成员生成 START -------------- %

% [S_Va_init_Combined,NNP,NG,MDA,M]=State_Ensem_Gen_Parellel_Common(S_Va_init,N_max,MStepNS_cov,targetPath);
[S_Va_init_Combined_eff,target_width,Val_Sen,N_new,NNP,NG,MDA,M]=State_Ensem_Gen_Parellel_Common_Beta_BL(S_Va_init,N_max,MStepNS_cov,targetPath); % 更新支持 扰动项Beta(x) 的同化

fprintf('集合成员生成过程 结束：最大集合成员数 = %d，有效集合成员数 = %d，敏感性指标 TMu_re = %.6f\n',N_max,N_new,Val_Sen);
fprintf('\n');
fprintf(fid_info,'集合成员生成过程 结束：最大集合成员数 = %d，有效集合成员数 = %d，敏感性指标 TMu_re = %.6f\n',N_max,N_new,Val_Sen);
fprintf(fid_info,'\n');

% 构建全状态变量矩阵，初始集合成员矩阵
S_Va_ana=S_Va_init_Combined_eff;

% --------------- 修改参数预置 + 集合成员生成 END --------------- %

% ----------------- 区域数据同化信息处理 START ----------------- %

% 若进行区域数据同化
if (flag_da_part==1)

    % 规则定义：仅规定同化区域内网格单元对应的涡粘系数进行同化；
    %           同化区域外的剩余网格单元涡粘系数保持原值；

    % 提取区域数据同化信息
    DA_area_info=load('.\Program_CFD_Base\Input\InputDAgrid.dat');
    % 第一行表示同化区域内网格节点
    NDAV=DA_area_info(1);
    % 第二行表示同化区域内单元数量
    NDAG=DA_area_info(2);
    % 提取同化区域内网格节点与单元编号信息
    DA_VN=DA_area_info(3:(NDAV+2));
    DA_GN=DA_area_info((NDAV+3):(NDAV+NDAG+2));
    
    % 标记出现在同化区域内的网格单元编号
    DA_GN_flag=zeros(NG,1);
    for i=1:NDAG
        DA_GN_flag(DA_GN(i))=1;
    end

    % 计算指定同化区域的新 MDA_new 与 M_new 值，其小于或等于原始值
    MDA_new=NDAG+PT;
    M_new=MDA_new+n_com;
else
    MDA_new=MDA;
    M_new=M;
end

% ------------------ 区域数据同化信息处理 END ------------------ %

% 保存初始集合信息
% 保存可供后续直接调用或续算的工作区数据
save([new_output_folder,savename,'_Pre.mat']);  % 保存集合成员生成过程所有结果
save([new_output_folder,savename,'_S_Va_init.mat'],'S_Va_ana'); % 保存更新成员结果至工作区
save([new_output_folder,'S_Va_init.dat'],'S_Va_ana','-ascii');
save([new_output_folder,'S_Va_init.mat'],'S_Va_ana');

%% 步骤四：建立实验与计算数据坐标映射关系
if ((strcmp(char(S_Va_ref_Name(1)),'Cp')==1) || (strcmp(char(S_Va_ref_Name(1)),'Cf')==1) || (strcmp(char(S_Va_ref_Name(1)),'yplus')==1))
    
    % 建立观测算子矩阵H
    H=zeros(k,M_new);
    for i=1:k
        H(i,(M_new-n_com)+num_da_arr(i))=1;  % !!! 通用程序中，观测算子 H 需修改
    end
    H_Cp=H(:,((M_new-n_com+1):M_new));

    % % !!! 重要修改，作用位置更改 !!!
    % 生成原始观测数据及观测数据误差集合
    D_obs_e_ori=zeros(k,N_max);
    Y_err_ori=normrnd(0,sigma_obs_ini,[k,N_max]);
    for i=1:k
        % D_obs_e(i,:)=(ones(1,N).*exp_data(i,2))+Y_err(i,:);
        % D_obs_e(i,:)=(ones(1,N).*shape_combined_arr_sor(num_da_arr(i),4))+Y_err(i,:);
        % D_obs_e_ori(i,:)=(ones(1,N_max).*shape_combined_arr_init(num_da_arr(i),7))+Y_err_ori(i,:);   % 真实观测值（可用插值得到）
        D_obs_e_ori(i,:)=(ones(1,N_max).*shape_combined_arr_init(num_da_arr(i),7));   % 真实观测值（可用插值得到）
    end

    % 设置模式误差协方差（本实验中可假设无过程噪声）
    std_q=0.0;
    sigma_2_q=std_q^2;
    C_qq_e=ones(M_new,M_new)*sigma_2_q;
    
elseif (strcmp(char(S_Va_ref_Name(1)),'BL')==1)
    
    % !!!!! %
    n_com=kb;
    % num_da_arr=ones(1,kb);
    num_da_arr=(1:kb); % !!! 注意原始问题BUG !!!
    
    % 建立观测算子矩阵H
    H=zeros(kb,M_new);
    for i=1:kb
        H(i,(M_new-n_com)+num_da_arr(i))=1;  % !!! 通用程序中，观测算子 H 需修改
    end
    H_Cp=H(:,((M_new-n_com+1):M_new));

    % % !!! 重要修改，作用位置更改 !!!
    % 生成原始观测数据及观测数据误差集合
    D_obs_e_ori=zeros(kb,N_max);
    Y_err_ori=normrnd(0,sigma_obs_ini,[kb,N_max]);
    for i=1:kb
        % D_obs_e(i,:)=(ones(1,N).*exp_data(i,2))+Y_err(i,:);
        % D_obs_e(i,:)=(ones(1,N).*shape_combined_arr_sor(num_da_arr(i),4))+Y_err(i,:);
        % D_obs_e_ori(i,:)=(ones(1,N_max).*shape_combined_arr_init(num_da_arr(i),7))+Y_err_ori(i,:);   % 真实观测值（可用插值得到）
        D_obs_e_ori(i,:)=ones(1,N_max).*Exp_Obs_Arr(i);   % 真实观测值（可用插值得到）
    end

    % 设置模式误差协方差（本实验中可假设无过程噪声）
    std_q=0.0;
    sigma_2_q=std_q^2;
    C_qq_e=ones(M_new,M_new)*sigma_2_q;
    
end

%% 步骤五：开始内迭代过程，数据同化算法实现
% 关键点：实现对CFD程序续算功能的控制
% Case 1：若关闭续算模式 JX_NS = 0，且 JX_Turb = 1，则仅代入 TMu(I) 网格单元数据
%                                      JX_Turb = 0，则不代入任何已经保存的流场网格单元数据
% Case 2：若开启续算模式 JX_NS = 1，且 JX_Turb = 1，则将代入 RU(I),U(I),V(I),P(I),TMu(I)（湍流工作变量） 网格单元数据
%                                      JX_Turb = 0，则将代入 RU(I),U(I),V(I),P(I) 网格单元数据
% 内迭代实现过程：
% Step 1：预测步 —— 代入初始的Ma,Alfa和TMu数据，计算更新 RU(I),U(I),V(I),P(I) 等变量，Ma,Alfa,TMu(I) 同化变量不变
% Step 2：滤波步 —— 代入滤波同化算法，根据翼型表面压力分布实验值与计算值的差异，更新同化变量 Ma,Alfa,TMu(I)
% Step 3：后处理 —— 对湍流工作变量中的病态数据进行修正
% 注：实际运行时，首先运行一定数量的预测步，后每隔一段时间进行滤波步和后处理，直至同化后的流场变量L2Norm收敛

% 数据交互模块编写
% 写出并保存状态变量，作为下一轮迭代输入

% ----------------------- 内迭代过程 START ------------------------ %

% p=0; % 无须赋值，程序初始化阶段已定义；但可手动调试，确定运行轮数后进行续算
while (p<=Max_inner_iter)

    % ------------- 修改参数 + 集合成员计算（预测步） START ------------- %
    
    if (p>0)
        N_pre=N_new;
        [S_Va_pre_eff,target_width,Val_Sen,N_new]=State_Ensem_Renew_Parellel_Common_Beta_BL(p,S_Va_ana,N_pre,MDA,MStepNS_pre,targetPath);
        S_Va_pre=S_Va_pre_eff;
        % 更新最新集合成员数
        fprintf(['集合成员内迭代第 ',num2str(p,'%03d'),' 轮结束：最大集合成员数 = %d，有效集合成员数 = %d，敏感性指标 TMu_re = %.6f\n'],N_max,N_new,Val_Sen);
        fprintf('\n');
        fprintf(fid_info,['集合成员内迭代第 ',num2str(p,'%03d'),' 轮结束：最大集合成员数 = %d，有效集合成员数 = %d，敏感性指标 TMu_re = %.6f\n'],N_max,N_new,Val_Sen);
        fprintf(fid_info,'\n');
    else
        S_Va_pre=S_Va_ana;
    end
    
    % ------------ 生成观测数据及其适应模型标准差矩阵 START ------------- %
    
    % !!! 重要修改，作用位置更改 !!!
    % 生成观测数据矩阵
    if ((strcmp(char(S_Va_ref_Name(1)),'Cp')==1) || (strcmp(char(S_Va_ref_Name(1)),'Cf')==1) || (strcmp(char(S_Va_ref_Name(1)),'yplus')==1))
        D_obs_e=zeros(k,N_new);
        D_obs_e_z=zeros(k,N_new); % 新的观测数据设置为零矩阵，目标同化为减小观测值与实验值（真值）之间的相对误差，理想为零（重要BUG）
        for i=1:k
            D_obs_e(i,:)=(ones(1,N_new).*shape_combined_arr_init(num_da_arr(i),7));   % 真实观测值（可用插值得到）
        end
    elseif (strcmp(char(S_Va_ref_Name(1)),'BL')==1)
        D_obs_e=zeros(kb,N_new);
        D_obs_e_z=zeros(kb,N_new); % 新的观测数据设置为零矩阵，目标同化为减小观测值与实验值（真值）之间的相对误差，理想为零（重要BUG）
        for i=1:kb
            D_obs_e(i,:)=ones(1,N_new).*Exp_Obs_Arr(i);   % 真实观测值（可用插值得到）
        end       
    end
    
    if ((strcmp(char(S_Va_ref_Name(1)),'Cp')==1) || (strcmp(char(S_Va_ref_Name(1)),'Cf')==1) || (strcmp(char(S_Va_ref_Name(1)),'yplus')==1))
        % 2021/04/27 植入自适应观测值标准差模型 - AOSDM（待测试）
        % 初始化模型参数
        sigma_obs_sum=0;
        % Case 0：关闭模型开关，固定观测值标准差为原始值
        if (flag_sigma_obs_adapt==0)
            sigma_obs_new=ones(k,1).*sigma_obs_ini;

        % Case 1：开启自适应模型（1 - 固定放缩模式）：则每轮同化前所有观测值标准差设置为上一轮标准差乘以放缩因子
        elseif (flag_sigma_obs_adapt==1)
            if (p==0)
                sigma_obs_new=ones(k,1).*sigma_obs_ini;
                sigma_obs_val=sigma_obs_ini;
            else
                sigma_obs_val=flag_sigma_obs_coeff*sigma_obs_val;
                sigma_obs_new=ones(k,1).*sigma_obs_val;
            end
            % 设置误差水平边界
            for i=1:k
                if (sigma_obs_new(i)<=sigma_obs_min)
                    sigma_obs_new(i)=sigma_obs_min;
                end
                if (sigma_obs_new(i)>=sigma_obs_max)
                    sigma_obs_new(i)=sigma_obs_max;
                end
            end

        % Case 2：开启自适应模型（2 - 区间平均模式），则每轮同化前根据上一轮所有观测值对应的计算值包络区间平均宽度，将所有观测值的标准差调整为平均期望标准差水平（乘以放缩因子）
        elseif (flag_sigma_obs_adapt==2)
            for i=1:k
                sigma_obs_sum=sigma_obs_sum+target_width(num_da_arr(i));
            end
            sigma_obs_ave=sigma_obs_sum/k;
            % 由 95% 置信区间对应 μ ± 1.96 σ，宽度设定为 3.92 σ
            % 计算得到新的观测值标准差取值
            sigma_obs_new=ones(k,1).*(flag_sigma_obs_coeff*(sigma_obs_ave/3.92));
            % 设置误差水平边界
            for i=1:k
                if (sigma_obs_new(i)<=sigma_obs_min)
                    sigma_obs_new(i)=sigma_obs_min;
                end
                if (sigma_obs_new(i)>=sigma_obs_max)
                    sigma_obs_new(i)=sigma_obs_max;
                end
            end

        % Case 3：若开启自适应模型（3 - 区间适应模式），则每轮同化前根据上一轮各观测值对应的计算值包络区间宽度，将各观测值的标准差调整为相应期望标准差水平（乘以放缩因子）
        elseif (flag_sigma_obs_adapt==3)
            sigma_obs_new=ones(k,1);
            for i=1:k
                sigma_obs_new(i)=flag_sigma_obs_coeff*(target_width(num_da_arr(i))/3.92);
                if (sigma_obs_new(i)<=sigma_obs_min)
                    sigma_obs_new(i)=sigma_obs_min;
                end
                if (sigma_obs_new(i)>=sigma_obs_max)
                    sigma_obs_new(i)=sigma_obs_max;
                end
            end

        end

        % 根据自适应误差计算调整结果，生成观测数据误差矩阵
        Y_err=zeros(k,N_new);
        for i=1:k
            for j=1:N_new
                Y_err(i,j)=normrnd(0,sigma_obs_new(i));
            end
        end
        % Y_err=normrnd(0,sigma_obs_init,[k,N_new]);
        C_ss_e=cov(Y_err'); % 观测误差协方差矩阵，维度为观测点数量k*k
        
    elseif (strcmp(char(S_Va_ref_Name(1)),'BL')==1)
        
        sigma_obs_new=ones(kb,1).*sigma_obs_ini;
        Y_err=zeros(kb,N_new);
        for i=1:kb
            for j=1:N_new
                Y_err(i,j)=normrnd(0,sigma_obs_new(i));
            end
        end
        % Y_err=normrnd(0,sigma_obs_init,[k,N_new]);
        C_ss_e=cov(Y_err'); % 观测误差协方差矩阵，维度为观测点数量k*k
        
    end
    
    % ------------- 生成观测数据及其适应模型标准差矩阵 END -------------- %
        
    % -- 备份操作 -- %
    save([new_output_folder,savename,'_S_Va_pre_restart_',num2str(p,'%03d'),'.mat']); % 保存截止目前（第 p+1 轮同化开始前）全部计算结果至工作区，可直接导入数据文件后开始续算
    % -- 备份操作 -- %
    
    % 将同化预测过程（Prior）获得的涡粘系数取对数处理（暂未启用）
    % for i=(PT+1):MDA
    %     for j=1:N
    %         S_Va_pre(i,j)=log10(S_Va_pre(i,j));
    %     end
    % end
    
    % ------------ 进行区域数据同化 [前处理]状态变量 START -------------- %
    
    % 若进行区域数据同化，则处理状态变量
    if (flag_da_part==1)
    
        % 处理全流场状态变量 S_Va_pre
        % !!!!! 修改 适应通用程序 !!!!!
        S_Va_pre_NDA=S_Va_pre; % 保存先验估计状态矩阵中非同化变量部分
        S_Va_pre_FDA=S_Va_pre; % 保存先验估计状态矩阵中需同化变量部分
        % S_Va_pre_NDA(1,:)=zeros(1,N_new);
        % S_Va_pre_NDA(2,:)=zeros(1,N_new);
        % S_Va_pre_FDA(1,:)=S_Va_pre(1,:);
        % S_Va_pre_FDA(2,:)=S_Va_pre(2,:);
        S_Va_pre_NDA((1:PT),:)=zeros(PT,N_new);  % 处理单同化变量
        S_Va_pre_FDA((1:PT),:)=S_Va_pre((1:PT),:);
        for i=1:NG
            % 若该网格节点位于指定同化区域，即 DA_GN_flag(i) = 1
            % 则将状态变量（NDA）中的对应涡粘系数行置零，剩余状态变量（FDA）从原始状态变量对应行复制
            if (DA_GN_flag(i)==1)
                S_Va_pre_NDA((i+PT),:)=zeros(1,N_new);
                S_Va_pre_FDA((i+PT),:)=S_Va_pre((i+PT),:);
            else
                S_Va_pre_NDA((i+PT),:)=S_Va_pre((i+PT),:);
                S_Va_pre_FDA((i+PT),:)=zeros(1,N_new);
            end
        end

        % 提取进行同化的状态变量数组
        S_Va_pre_FDA_sum_flag=sum(S_Va_pre_FDA,2); % 标记 S_Va_pre_FDA 位置数组：若整行加和为零，则为非同化变量行；否则为同化变量行
        S_Va_pre_DA=zeros(M_new,N_new);
        DA_pre_cnt=0; % 初始化计数器
        % 遍历原始状态变量矩阵
        for i=1:M
            if (S_Va_pre_FDA_sum_flag(i)~=0)
                DA_pre_cnt=DA_pre_cnt+1;
                S_Va_pre_DA(DA_pre_cnt,:)=S_Va_pre_FDA(i,:); % 提取进行同化的状态变量数组，预期行数 = M_new
            end
        end
        
    else
        
        % 若不进行区域数据同化，则提取得到的进行同化的状态变量数组仍为原始数组
        S_Va_pre_DA=S_Va_pre;
        
    end
    
    % ------------- 进行区域数据同化 [前处理]状态变量 END --------------- %
    
    % -------------- 修改参数 + 集合成员计算（预测步） END -------------- %
    
    if (flag_da_method==1)
        
        % -------------- 核心算法①：EnKF 滤波方法更新集合成员 START ---------------- %
        
        % EnKF算法实现
        % 注意：EnKF运行过程中，需要保存的矩阵的最大维度为[M * M]，故当状态变量维度过大时，EnKF无法使用！
        % 步骤一：预测操作
        % 计算该时刻预测的集合误差协方差 C_pp_e_f -- dim[M * M]
        C_pp_e_f=cov(S_Va_pre_DA')+C_qq_e;
        
        % 步骤二：更新操作，引入观测信息
        % ① 计算卡尔曼增益Ke
        Ke=C_pp_e_f*H'*pinv(H*C_pp_e_f*H'+C_ss_e);
        
        % ② 引入观测信息，更新集合预测值为分析值
        % 对集合成员进行更新
        S_Va_ana_DA=S_Va_pre_DA+Ke*(D_obs_e-H*S_Va_pre_DA);
        
        % ③ 更新误差协方差
        C_pp_e_a=(ones(M_new,M_new)-Ke*H)*C_pp_e_f;
        % C_pp_e_a=cov(S_Va_ana_DA');
        
        % --------------- 核心算法①：EnKF 滤波方法更新集合成员 END ----------------- %
            
    elseif (flag_da_method==2)
        
        % -------------- 核心算法②：ETKF 滤波方法更新集合成员 START ---------------- %

        % 同化操作：提供相应位置的观测信息，真实模型向右传播
        % ETKF算法实现
        % 计算该时刻预测的集合误差协方差
        % C_pp_e_f=cov(S_Va_pre')+C_qq_e;

        % 步骤一：奇异值分解
        S_Va_pre_DA_m=(S_Va_pre_DA*ones(N_new,N_new))/N_new;
        S_Va_pre_DA_det=S_Va_pre_DA-S_Va_pre_DA_m;
        R=eye(N_new)+(S_Va_pre_DA_det')*(H')*pinv((N_new-1).*C_ss_e)*H*S_Va_pre_DA_det; % 调用 Moore-Penrose 伪逆，保证在秩亏情况下同化的稳定性
        [ZL,G,ZR]=svd(R); % 奇异值分解

        % 步骤二：更新操作，引入实验观测信息
        % ① 计算卡尔曼增益Ke，维度为M*k
        Ke=S_Va_pre_DA_det*ZL*inv(G)*ZR'*(S_Va_pre_DA_det')*(H')*pinv((N_new-1).*C_ss_e);
        % ② 计算变换矩阵T
        T=ZL*sqrtm(inv(G))*ZR';
        % ③ 引入观测信息，更新集合预测值为分析值
        % 对状态集合成员进行更新
        if (flag_da_relative==1)
            S_Va_ana_DA_m=S_Va_pre_DA_m+Ke*(D_obs_e-H*S_Va_pre_DA_m);
        elseif (flag_da_relative==0)
            S_Va_ana_DA_m=S_Va_pre_DA_m+Ke*(D_obs_e_z-H*S_Va_pre_DA_m);
        end
        S_Va_ana_DA_det=S_Va_pre_DA_det*T;
        S_Va_ana_DA=S_Va_ana_DA_m+S_Va_ana_DA_det;
        % 更新误差协方差
        % C_pp_e_a=(eye(M)-Ke*H)*C_pp_e_f;

        % --------------- 核心算法②：ETKF 滤波方法更新集合成员 END ----------------- %
        
    end

    % ------------- 进行区域数据同化 [还原]状态变量 START --------------- %
    
    % 若进行区域数据同化，则须还原状态变量
    if (flag_da_part==1)
        
        % 还原同化后的全流场状态变量数组，长维度为原始状态变量矩阵对应维度
        S_Va_ana_FDA=zeros(M,N_new);
        DA_ana_cnt=0; % 初始化计数器
        for i=1:M
            if (S_Va_pre_FDA_sum_flag(i)~=0)
                DA_ana_cnt=DA_ana_cnt+1;
                S_Va_ana_FDA(i,:)=S_Va_ana_DA(DA_ana_cnt,:); % 提取进行同化的状态变量矩阵
            end
        end
        
        % !!! 同化状态变量矩阵 + 原始非同化状态变量矩阵 = 叠加得到最终全流场状态变量分析矩阵（维度为 M * N_new）!!!
        S_Va_ana=S_Va_ana_FDA+S_Va_pre_NDA;
        
    else
        
        % 若不进行区域数据同化，则同化后得到的状态变量数组即为全流场状态变量
        S_Va_ana=S_Va_ana_DA;
        
    end
    
    % -------------- 进行区域数据同化 [还原]状态变量 END ---------------- %
    
    % -------------- 后处理步骤 + 修正病态涡粘系数 START ---------------- %
    
    for kk=1:S_Va_target_Num
        
        % 注意：仅当同化涡粘时，即状态变量中存在涡粘矩阵，对病态数据进行修正，其他情况无需修正
        if (strcmp(char(S_Va_target_Name(kk)),'TMu')==1)
            
            % 对同化滤波后（analytical）的湍流涡粘系数中的病态数据进行修正，并恢复原始系数值
            for i=(PT+1):MDA
                for j=1:N_new

                    % 修正病态数据
                    if (S_Va_ana(i,j)<1E-10)
                        S_Va_ana(i,j)=1E-10;
                    end

                    % !!!! 取对数策略暂未启用
                    % if (S_Va_ana(i,j)<-10)
                    %     S_Va_ana(i,j)=-10;
                    % end
                    % S_Va_ana(i,j)=10^(S_Va_ana(i,j));  % 转换后有部分精度损失

                end
            end
            
        % elseif (strcmp(char(S_Va_target_Name(kk)),'Beta')==1)
             % 修正项 Beta场 可视化（Visualization）
             % 内迭代过程提取 S_Va_ana 每轮结果中的有效数据并进行可视化
        
        end
        
    end
    
    % --------------- 后处理步骤 + 修正病态涡粘系数 END ----------------- %
    
    % 步骤三：记录当前最优分析估计 S_Va_ana_ave 与 S_Va_ana_e
    S_Va_ana_ave=mean(S_Va_ana'); % 计算得到最优估计值
    S_Va_ana_ave_da=S_Va_ana_ave((PT+1):MDA)'; % 提取同化后除观测部分的状态变量（如TMu场、Beta场等）!!! 0418_修复BUG 设置MDA !!!
    
    % 每轮同化结果 修正项Beta场 可视化
    for kk=1:S_Va_target_Num
        if (strcmp(char(S_Va_target_Name(kk)),'Beta')==1)
            % 利用MATLAB自带的 scatteredInterpolant() 函数对同化得到的网格单元值进行曲面插值至网格节点，构建 V = F(x,y) 的拟合曲面
            V=scatteredInterpolant(GN_grid_x,GN_grid_y,S_Va_ana_ave_da);
            V.Method='natural';
            V.ExtrapolationMethod='nearest';
            Va_NNP_Beta_new=V(VN_grid_x,VN_grid_y); % 插值得到所有网格节点 Beta 值，利用网格单元值插值
        end
    end
    
    % 保存当前计算信息
    target_renew_dat=['S_Va_round_',num2str(p+1,'%03d'),'.dat'];
    target_renew_mat=['S_Va_round_',num2str(p+1,'%03d'),'.mat'];
    target_renew_da_dat=['S_Va_round_da_',num2str(p+1,'%03d'),'.dat'];
    target_renew_da_mat=['S_Va_round_da_',num2str(p+1,'%03d'),'.mat'];
    save([new_output_folder,target_renew_dat],'S_Va_ana','-ascii'); % 保存更新成员结果至工作区
    save([new_output_folder,target_renew_mat],'S_Va_ana'); % 保存更新成员结果至工作区
    save([new_output_folder,target_renew_da_dat],'S_Va_ana_DA','-ascii'); % 保存更新同化成员结果至工作区
    save([new_output_folder,target_renew_da_mat],'S_Va_ana_DA'); % 保存更新同化成员结果至工作区
    
    % 设置同化退出标志
    % 满足 L2Norm条件 或 内循环操作达到最大次数，即退出内循环过程
    % 待编写ING。。。
    
    % 内迭代数目计数
    p=p+1;
    
    % ----------------------- 重计算过程 START ------------------------ %

    S_Va_ana_final=S_Va_ana_ave';
    DA_Results_targetPath=[targetPath,'DA_Final_Results'];
    State_Ensem_Recal_New_Common_Beta_BL(p,S_Va_ana_final,MStepNS_recal_ensemble,DA_Results_targetPath);
    save([savename,'_DA_RUN_All.mat']);  % 保存当前同化过程所有结果

    % ----------------------- 重计算过程 END ------------------------ %
    
end

% ----------------------- 内迭代过程 END ------------------------ %

%% 步骤六：数据同化结果后处理
save([savename,'_All.mat']);  % 保存同化过程所有结果
fprintf('\n');
fprintf(fid_info,'\n');
fprintf('数据同化结果后处理开始\n');
fprintf(fid_info,'数据同化结果后处理开始\n');
fprintf('\n');
fprintf(fid_info,'\n');

% 绘制内迭代过程气动力收敛曲线
new_output_folder_cov=[new_output_folder,'Results_Recal_Covergence_Analysis'];
mkdir(new_output_folder_cov);
targetPath=new_output_folder_cov;

fid_cal_force_all=fopen([new_output_folder_cov,'\Recal_Covergence_Analysis_Ensemble_Force_ALL.PLT'],'w+');

fprintf('内迭代过程气动力收敛曲线结果正在导出\n');
fprintf(fid_info,'内迭代过程气动力收敛曲线结果正在导出\n');

for i=1:N_new
    Cal_force_init=load([savename,'\Pre_Ensemble_Generation\Init_Ensemble_Member_',num2str(i,'%03d'),'\Force.PLT']);
    Cal_force_inner_iter=Cal_force_init;
    for p=1:Max_inner_iter
        Cal_force_new=load([savename,'\Inner_Iteration_Round_',num2str(p,'%03d'),'\Renew_Ensemble_Member_',num2str(i,'%03d'),'\Force.PLT']);
        Cal_force_temp=[Cal_force_inner_iter;Cal_force_new];
        Cal_force_inner_iter=Cal_force_temp;
    end
    % 保存内迭代过程气动力收敛数据
    for j=1:length(Cal_force_inner_iter)
        Cal_force_inner_iter(j,1)=100*j;
    end
    
    fid_cal_force=fopen([new_output_folder_cov,'\Recal_Covergence_Analysis_Ensemble_Force_',num2str(i,'%03d'),'.PLT'],'w+');
    fprintf(fid_cal_force,'VARIABLES="Residual","Cl","Cd","Cm"\n');
    fprintf(fid_cal_force,'ZONE T="ZONE %03d"\n',i);
    for j=1:length(Cal_force_inner_iter)
        fprintf(fid_cal_force,'%12d %16.5f %16.5f %16.5f\n',Cal_force_inner_iter(j,1),Cal_force_inner_iter(j,2),Cal_force_inner_iter(j,3),Cal_force_inner_iter(j,4));
    end
    fclose(fid_cal_force);
    
    fprintf(fid_cal_force_all,'VARIABLES="Residual","Cl","Cd","Cm"\n');
    fprintf(fid_cal_force_all,'ZONE T="ZONE %03d"\n',i);
    for j=1:length(Cal_force_inner_iter)
        fprintf(fid_cal_force_all,'%12d %16.5f %16.5f %16.5f\n',Cal_force_inner_iter(j,1),Cal_force_inner_iter(j,2),Cal_force_inner_iter(j,3),Cal_force_inner_iter(j,4));
    end
    fprintf(fid_cal_force_all,'\n');
    
    fprintf('集合成员 %03d 气动力收敛曲线结果导出完成\n',i);
    fprintf(fid_info,'集合成员 %03d 气动力收敛曲线结果导出完成\n',i);
    pause(2);

end
fprintf('内迭代过程气动力收敛曲线结果已全部导出\n');
fprintf(fid_info,'内迭代过程气动力收敛曲线结果已全部导出\n');
fprintf('\n');
fprintf(fid_info,'\n');

fclose(fid_cal_force_all);

% 涡粘系数方差系数对比
% 定义方差系数（对同化区域内各网格节点）
% 读取最终同化后各成员涡粘系数数据
fprintf('涡粘系数方差系数结果正在导出\n');
fprintf(fid_info,'涡粘系数方差系数结果正在导出\n');
% 处理同化前涡粘数据
TMuSP=zeros(NNP,N_new);
TMu_coeff=zeros(NNP,1); % 保存涡粘系数方差系数
for i=1:N_new
    TMu_temp=load([savename,'\Pre_Ensemble_Generation\',initPath,'_',num2str(i,'%03d'),'\TMuPData.PLT']);
    for m=1:NNP
        if (TMu_temp(m)<1E-5)
            TMu_temp(m)=1E-5;
        end
    end
    TMuSP(:,i)=TMu_temp;
end
TMu_mean=(mean(TMuSP'))';
for m=1:NNP
    TMu_err_sum=0;
    for i=1:N_new
        TMu_err_sum=TMu_err_sum+((TMuSP(m,i)-TMu_mean(m))^2);
    end
    % 计算涡粘系数方差系数
    TMu_coeff(m)=(sqrt(TMu_err_sum/(N_new-1)))/TMu_mean(m);
end
% 替换原始网格节点涡粘系数数据文件
save([ddir_post,'\Input\TMuPData.PLT'],'TMu_coeff','-ascii');
% 调用 CFD_2D_Post_TMu 程序，并保存输出流场云图
cmd_post=[ddir_post,'\CFD_2D_Post_TMu.exe'];
open(cmd_post);
pause(5);
pre_post_TMu_Path=[ddir_post,'\Output\TMuPData_Post.PLT'];
copyfile(pre_post_TMu_Path,[new_output_folder,'Pre_Ensemble_Generation\TMuPData_Post.PLT']);
fprintf('同化前状态 涡粘系数方差系数结果导出完成\n');
fprintf(fid_info,'同化前状态 涡粘系数方差系数结果导出完成\n');
pause(2);
% 处理各内迭代步的涡粘数据
for p=1:Max_inner_iter
    TMuSP=zeros(NNP,N_new);
    TMu_coeff=zeros(NNP,1); % 保存涡粘系数方差系数
    for i=1:N_new
        TMu_temp=load([savename,'\Inner_Iteration_Round_',num2str(p,'%03d'),'\Renew_Ensemble_Member_',num2str(i,'%03d'),'\TMuPData.PLT']);
        for m=1:NNP
            if (TMu_temp(m)<1E-5)
                TMu_temp(m)=1E-5;
            end
        end
        TMuSP(:,i)=TMu_temp;
    end
    TMu_mean=(mean(TMuSP'))';
    for m=1:NNP
        TMu_err_sum=0;
        for i=1:N_new
            TMu_err_sum=TMu_err_sum+((TMuSP(m,i)-TMu_mean(m))^2);
        end
        % 计算涡粘系数方差系数
        TMu_coeff(m)=(sqrt(TMu_err_sum/(N_new-1)))/TMu_mean(m);
    end
    % 替换原始网格节点涡粘系数数据文件
    save([ddir_post,'\Input\TMuPData.PLT'],'TMu_coeff','-ascii');
    % 调用 CFD_2D_Post_TMu 程序，并保存输出流场云图
    cmd_post=[ddir_post,'\CFD_2D_Post_TMu.exe'];
    open(cmd_post);
    % system([ddir_post,'\CFD_2D_Post_TMu.exe']);
    pause(5);
    pre_post_TMu_Path=[ddir_post,'\Output\TMuPData_Post.PLT'];
    copyfile(pre_post_TMu_Path,[new_output_folder,'Inner_Iteration_Round_',num2str(p,'%03d'),'\TMuPData_Post.PLT']);
    fprintf('内迭代状态 %03d 涡粘系数方差系数结果导出完成\n',p);
    fprintf(fid_info,'内迭代状态 %03d 涡粘系数方差系数结果导出完成\n',p);
    pause(2);
end
fprintf('涡粘系数方差系数结果已全部导出\n');
fprintf(fid_info,'涡粘系数方差系数结果已全部导出\n');
fprintf('\n');
fprintf(fid_info,'\n');

fclose(fid_info);
fclose(fid_all);
