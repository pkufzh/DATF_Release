%% ����ͬ�������������Data Assimilation & Turbulent Flows��
%% ���ϱ任�������˲�ʵ�飨Ensemble Transform Kalman Filter��ETKF����
%% ����ETKF�������ʵ�����ݣ�ʵ�ֶ���������������������ͬ��
%% �����¡�2021/06/04�������򡪡�ָ����ʼ�Ŷ�������Ŀ��ͬ����������������ͬ����ͨ�ó��� Beta_0510��
%% �����¡�2021/06/04����� & �����Լ��
%% ��Ҫ���ܣ��������ɼ��ϳ�Ա����������ͬ����������
%% �޸�BUG���� �������ڵ���������Ҫ�ָ���׼ֵ�ı������������ж�����޸ģ�
%%          �� ��̬�ֲ������������� X ~ N(mu,sigma)��
%% ���¹��ܣ�
%%          �� ��ճϵ��ȡ����������log�����ԣ���BUG����
%%          �� ���ӡ�����ͬ������/������ģ�飻
%%          �� ������������ͬ��ģ�飨������Beta����TMuͬ�� ����
%%          �� ����ָ����ʼ�Ŷ�������Ŀ��ͬ������ģ�飻
%%          �� ����ͬ���˲����̲�����������������
%%          �� ���ӶԼ���ѹ��ϵ�����߽ڵ��Ž���������ʵʱ�����ǰͬ���ִ�ѹ���ֲ����ߣ������ߣ���
%%          �� �����޶���ʼ���ϳ�Ա�����������������������ܽ������ͬ�����̣���ʵͬ�����ϳ�Ա���� N_new �� N_max��
%%          �� ���Ӷ����� S-Aģ�������� Beta ����Ŷ����ԣ��������������Բ������Ӧģ�飻
%%          �� ���¿�ʹ�����������ͬ�����൱�ڽ�ͬ��Ŀ��ת��Ϊʹ�����������㣨������ҪBUG����
%%          �� ����Beta�����ӻ���ͬ����Ԫ���ݡ���ֵ��ͬ���ڵ����ݣ���
%%          �7�6 ���ӱ��ȫͬ���������㿪�أ�
%%          �7�7 ��������Ӧ�۲���׼�����ģ�ͣ������ԣ���
%%          �7�8 ���ӶԷ�����������ǰ����֧�֣�
%%          �7�9 ���������Է���ָ����㣨2021/05/03����
%%          �7�0 ���ӱ߽���ٶ���ͬ�����ܣ�2021/05/18 ��ʼ ����
%% ��Ҫ���£���ѡ��ʹ�ü��Ͽ������˲���EnKF���򼯺ϱ任�������˲���ETKF��
%% ����һ��RAE2822���ͣ���� CFD ���� 
%% ��׼״̬����ŵ�� Re = 650W������� Ma = 0.729��ӭ�� AoA = 0.310deg
%% ��������S809���ͣ���� CFD ����
%% ��׼״̬����ŵ�� Re = 200W������� Ma = 0.200��ӭ�� AoA = 10.2,11.2,12.2,14.0deg
%% ��������NACA4412���ͣ���� CFD ����
%% ��׼״̬����ŵ�� Re = 152W������� Ma = 0.09������ѹ״̬����ӭ�� AoA = 13.87deg

%% ������ʵ��
%% ����ͬ������ȫ�����㿪��
flag_da_xusuan=0; % �����Ƿ����ȫͬ����������
                  % 0 - ���������㣻
                  % 1 - �����ж����㣬�ɵ��������жϺ󱣴�.mat�ļ������ɵ���ͬ����������ʼ�������㣻
                  % 2 - ֱ�Ӷ�λ��ָ��ͬ������ p �ı����ļ�����ʼ���㣻
%% ��ʼ������
% �����������㣬���������й���������
if (flag_da_xusuan==0)
    clear all
    clc
    clf
    close all
    flag_da_xusuan=0;
    p=0; % ��ʼ��ͬ������ p
elseif (flag_da_xusuan==2)
    clear all
    flag_da_xusuan=2;
    p=1; % ��λ��ָ��ͬ������ p������ͬ��
end

% �����ʼ��Ա����
% load('Pre_Cal_Path.mat');

format long;

% ����ȫ�ֱ���
global savename new_output_folder fid_all fid_info PPNS initPath CFL_Setting Resi_target Resi_target_lev JX_NS STD_Comp JX_Turb L2NORM_Setting flag_Va_init;
global S_Va_turb_Name S_Va_turb_Type S_Va_turb_Num S_Va_turb_FilePath S_Va_turb_ModiLine S_Va_turb_ModiFormat S_Va_turb_Std;
global S_Va_target_Name S_Va_target_Num S_Va_target_TypeTurb S_Va_target_TypeDA S_Va_turb_target S_Va_Flow_length PT PS;
global S_Va_ref_Name S_Va_ref_Num;
global shape_upper shape_lower shape_outline shape_combined_arr_init Va_Beta Va_NNP_Beta Va_NNP_Beta_new;
global flag_da_relative flag_sigma_obs_adapt flag_sigma_obs_coeff D_obs_e_ori D_obs_e flag_flow_type;
global BL_setting N_BL L_BL BL_num BL_lab BL_dot BL_obs_num kb;
global Resi_flag Resi_flag_label;

%% ������Ҫ����ͬ�����̲���

% ------------  ������Ҫ�ļ������벢����ز��� START ------------ %

% �������ݽ���·��
savename='Results_20210709_Case_RAE2822_Nmax_128_P_420_DA_81';
ddir_recal=('Program_CFD_Parallel_Recal');  % �����ؼ���·��
ddir_post=('Program_CFD_Parallel_Post');    % ���ý������·��
new_output_folder=['.\',savename,'\'];

% �����������㿪�أ���ֱ�ӽ������ļ��У����������
if (flag_da_xusuan==0)
    mkdir(new_output_folder);
    mkdir(['.\Output_Figures\',savename]);
end
% ������ָ�����㿪�أ�ֱ�Ӷ�λ��ָ��ͬ������ p �ı����ļ�����ʼ���㣨ע�⣺ָ��ͬ�����������Ѿ�������ϣ���
if (flag_da_xusuan==2)
    load([new_output_folder,savename,'_S_Va_pre_restarat_',num2str(p,'%03d'),'.mat']);
end
% ����/����ģ���������

targetPath=new_output_folder;
initPath='Init_Ensemble_Member';

% �򿪼�¼������Ϣ�ļ�
fid_info=fopen([new_output_folder,'DA_CompInfo_Brief.dat'],'w+');
fid_all=fopen([new_output_folder,'DA_CompInfo_Detailed.dat'],'w+');

% ���ò����߳�/���������ɸ��ݹ������߳����ã�
PPNS=8;

% -------------  ������Ҫ�ļ������벢����ز��� END ------------- %

% -------------  �������ɳ�ʼ���ϳ�Ա��ز��� START ------------- %

% ���弯�ϳ�Ա���������ע�⣺ʵ�ʵ� N_new �� N��
N_max=128;

% ----------  ָ����ʼ�Ŷ�������Ŀ��ͬ��������ͬ���۲����ģ�� START ---------- %

% �� ָ����ʼ�Ŷ�����ģ�飨ĿǰĬ��Ϊ�������Ŷ���
% ������׼״̬����
% �����ʼ�Ŷ��������ơ��ļ�·�����޸��������޸ĸ�ʽ�����������ʾ��һһ��Ӧ��
% ��ʼ�Ŷ���������
S_Va_turb_Name=["Ma","AoA","CKS"];
% S_Va_turb_Name=["Beta"];
S_Va_turb_CharName=char(S_Va_turb_Name);
% ��ʼ�Ŷ��������ͣ�0 - ��������������ϣ�1 - �ռ�ֲ�����
S_Va_turb_Type=0;
% ��ʼ�Ŷ���������
S_Va_turb_Num=length(S_Va_turb_Name);
% ��ʼ�Ŷ������ļ�·����ʡ�����г����Ŀ¼��
% S_Va_turb_FilePath_raw=char('\Input\InputNS.dat','\Input\InputNS.dat','\Input\InputTurb.dat');
% S_Va_turb_FilePath=strip(S_Va_turb_FilePath_raw,'right');
S_Va_turb_FilePath=["\Input\InputNS.dat","\Input\InputNS.dat","\Input\InputSAconst.dat"];
% ��ʼ�Ŷ������޸���������Ӧ�ļ��ڣ�
S_Va_turb_ModiLine=[1,2,7];
% ��ʼ�Ŷ������޸ĸ�ʽ����Ӧ�ļ��ڣ�
S_Va_turb_ModiFormat=["%-16.5f","%-16.5f","%-16.5f"];
% ��ʼ�Ŷ������Ŷ�����
% �����Ŷ���ʽ��0 - ����ֵ�����Ŷ���1 - ���������Ŷ������磺����50�����Ŷ�����Ϊ��׼ֵ��[50%,150%]��
S_Va_turb_way=[0,0,1];
S_Va_turb_Std=[0.729,2.310,0.410];  % ��ʼ�Ŷ�������׼ֵ
S_Va_turb_Det_base_lower=[0.029,0.310,20];  % ��ʼ�Ŷ�����������
S_Va_turb_Det_base_upper=[0.029,0.310,20];  % ��ʼ�Ŷ�����������

S_Va_turb_Det_real_lower=S_Va_turb_Det_base_lower;
S_Va_turb_Det_real_upper=S_Va_turb_Det_base_upper;
for i=1:length(S_Va_turb_Std)
    % �����ñ��������Ŷ���ʽ
    if (S_Va_turb_way(i)==1)
        S_Va_turb_Det_real_lower(i)=S_Va_turb_Std(i)*(S_Va_turb_Det_base_lower(i)/100);
        S_Va_turb_Det_real_upper(i)=S_Va_turb_Std(i)*(S_Va_turb_Det_base_upper(i)/100);
    end
end
S_Va_turb_min=S_Va_turb_Std-S_Va_turb_Det_real_lower;
S_Va_turb_max=S_Va_turb_Std+S_Va_turb_Det_real_upper;
% �����������������г�ʼ�Ŷ���������
Va_min=S_Va_turb_min;
Va_max=S_Va_turb_max;
% ��ʼͬ����������
PN=length(Va_min);
S_Va_coef=lhsdesign(N_max,PN);
% ����N*P���Ŷ���������
S_Va_init=(ones(N_max,1)*Va_min)+S_Va_coef.*(ones(N_max,1)*(Va_max-Va_min));
S_Va_init_ave=mean(S_Va_init);

% �� ָ��Ŀ��ͬ������ģ��
% ȫ��Ŀ��ͬ����������
S_Va_target_Name=["Ma","AoA","TMu"];
% S_Va_target_Name=["Beta"];
S_Va_target_CharName=char(S_Va_target_Name);
% ȫ��Ŀ��ͬ����������
S_Va_target_Num=length(S_Va_target_Name);
% �ж�Ŀ��ͬ���������͢٣�TypeTurb����k - �ǵ�k����ʼ�Ŷ�������ֲ���k~=0����0 - �ǳ�ʼ�Ŷ�����
% ע����Ŀ��ͬ�������ǳ�ʼ�Ŷ�����������Ҫ�ڸ��²���ֵ��Ϊԭʼ STD ֵ
S_Va_target_TypeTurb=[1,2,0];
% �ж�Ŀ��ͬ���������͢ڣ�TypeDA����1 - ����ֵ��ͬ��������2 - ��������ֵ��ͬ������
S_Va_target_TypeDA=[1,1,2];
% Ŀ�굥ֵͬ������Ϊ��ʼ�Ŷ����ı�������Ϊ PT����������ֵ��ͬ���������� PS����
% �����ڱ������У�����ͬ������� Ma��ӭ�� AoA��ȫ�����֣�������ճ TMu���� PT = 2��
%                   ����ͬ��ӭ�� AoA��ȫ�����֣�������ճ TMu���� PT = 1��
% PT=sum((S_Va_target_TypeTurb>0) & (S_Va_target_TypeDA==1));
PT=sum((S_Va_target_TypeDA==1));
PS=sum((S_Va_target_TypeDA==2));
% ��ǳ�ʼ�Ŷ�������ͬʱ����[ΪĿ��ͬ��]��[��ֵ]�ı�����e.g. "Ma" "AoA"������������S_Va_target_NonTurb����
Lia=ismember(S_Va_turb_Name,S_Va_target_Name);
S_Va_turb_target=double(Lia);

% �� ָ��ͬ���۲����ģ��
% ȫ��ͬ���۲��������
S_Va_ref_Name=["Cp"];
% S_Va_ref_Name=["BL"];
S_Va_ref_CharName=char(S_Va_ref_Name);
% ȫ��ͬ���۲��������
S_Va_ref_Num=length(S_Va_ref_Name);

% ��Ŀ��۲����ΪCp,Cf,yplus    
    % ���ò���ͬ�����̵��������ݵ���
    % num_da_min=1;    % ȡ����Сֵ
    % num_da_max=97;  % ȡ�����ֵ
    % num_da_gap=2;    % ȡ�����
    % ����ͬ�����ݵ��������
    % num_da_arr=(num_da_min:num_da_gap:num_da_max);
    % �ɸ����ص�������б�־
    % num_da_arr=[(10:2:42),(48:64),(65:73),(74:2:92)];
    % num_da_arr=[(10:2:42),(48:92)];
    % num_da_arr=[(10:2:42),(48:73),(74:2:92)]; % RAE2822
    % num_da_arr=[(18:3:63),(90:2:168)]; % 56 Points Selected for S809 ��P = 199��
    % num_da_arr=[(20:2:60),(94:2:166)]; % 58 Points Selected for S809 ��P = 199��
    num_da_arr=[(12:7:75),(76:5:186),(208:5:288),(290:3:320),(321:5:376),(377:4:405)]; % 81 Points Selected for RAE2822 (P = 420)
    % ��¼ͬ�����ݵ����
    k=length(num_da_arr);

% ��Ŀ��۲����Ϊ�߽���ٶ��ͣ�����Ҫָ���߽��վλ�����������վλ�۲���Ϣ��y/c,u/u_ref������վλ���߽��ͬ��λ��
    % N_BL=4;
    % BL_Type=[1,1,1,1];  % ���ø�վλ�߽��ȡ��ģʽ��0 - ָ����ʼ����ֹ����Լ������1 - ָ��ͬ������
    % BL_S=[]; % ���ø�վλ����ֵ��ʼ���
    % BL_T=[]; % ���ø�վλ����ֵ��ֹ���
    % BL_G=[]; % ���ø�վλ����ֵ��ż��
    % global BL_setting N_BL L_BL BL_num BL_lab BL_dot
    BL_setting=textread('BL_DA_Setting.dat'); % ����߽��ͬ�������ļ�
    N_BL=size(BL_setting,1);
    L_BL=size(BL_setting,2);
    % BL_dot=(-1)*ones(N_BL,100); % ����߽��ͬ��λ����Ϣ
    BL_num=zeros(1,N_BL); % �����վλ�߽��ͬ������
    BL_lab=zeros(1,N_BL); % �������ͬ���ı߽����
    for i=1:N_BL
        BL_lab(i)=BL_setting(i,1); % ����߽��վλ���
        BL_type=BL_setting(i,2); % ���ø�վλ�߽��ȡ��ģʽ��0 - ָ����ʼ����ֹ����Լ������1 - ָ��ͬ������
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
    
    % Pos_Obs=(-1)*ones(N_BL,100); % ����߽��ͬ��λ����Ϣ
    % ����ָ���߽����ٶ�����Ϣ��ʵ������ + λ����Ϣ��
    kb=sum(BL_num); % ͳ�Ʋ���ͬ���ı߽���ٶ��͵�����
    Exp_Obs_Arr=zeros(1,kb);
    BL_obs_num=0;
    for i=1:N_BL
        Vel_Exp_Data=load(['.\Input_Data\Input_Exp_Hump_BL_Vel_',num2str(BL_lab(i),'%03d'),'.dat']);
        Vel_Shape_Data=load(['.\Input_Data\Input_Hump_BL_Profile_',num2str(BL_lab(i),'%03d'),'.dat']);
        
        % ����ָ�����ͬ�����λ����Ϣ��y/c��
        for j=1:BL_num(i)
            Pos_Obs(i,j)=Vel_Shape_Data(BL_dot(i,j));
            % ��ֵ�õ��߽��ָ��λ�õ��ٶ�ʵ������
            Exp_Obs(i,j)=interp1(Vel_Exp_Data(:,1),Vel_Exp_Data(:,2),Pos_Obs(i,j),'linear','extrap');
            BL_obs_num=BL_obs_num+1; % ͬ���۲�� + 1
            Exp_Obs_Arr(BL_obs_num)=Exp_Obs(i,j); % �������в�ֵ����ʵ��۲����ݣ����� = k
        end
        
    end
    

% -----------  ָ����ʼ�Ŷ�������Ŀ��ͬ��������ͬ���۲����ģ�� END ----------- %

% --------------  �������ɳ�ʼ���ϳ�Ա��ز��� END -------------- %

% -----------------  ������Ҫ���������� START ----------------- %

MStepNS_cov=12000;   % ����NS��������������������������ɹ��̣�
Max_inner_iter=10;  % �����ڵ�������������
MStepNS_pre=10000;   % �����ڵ���Ԥ�ⲽ��
MStepNS_recal_ensemble=18000;  % ����ÿ����Ա�ؼ��㲽���������ؼ�����̣�
CFL_Setting=[5,5,5]; % �ֱ����ü��ϳ�Ա���ɡ��ڵ����Լ��ؼ�����̵� CFL ��
Resi_target_log=[-9.0,-9.0,-9.0]; % �ֱ��弯�ϳ�Ա���ɡ��ڵ������ؼ�����̵��������ȣ��� 10 Ϊ�ף�
Resi_target=10.^(Resi_target_log);
Resi_target_lev=[-4.0,-4.0,-4.0]; % �涨���׶Σ��������ɡ��ڵ������ؼ��㣩�Ľضϲв�ˮƽ����ȷ�����ϳ�Ա��׽��Χ
% Resi_target_gen_level=Resi_target_log(1)+1.5; % �涨�������ɽضϲв�ˮƽΪ�����в�ֵ + f �������󼯺ϳ�Ա��׽��Χ
JX_NS=[0,1,0]; % ���£��رգ�JX_NS ������ݿ��أ������������Ƿ����㣺0 - �񣨴ӳ�ʼ������ʼ���㣩��1 - �� ����ǰһ�ּ���������ʼ���㣩

JX_Turb=[0,1,1]; % ���£��رգ�JX_Turb ������ݿ��أ���������ճ���Ƿ����㣺0 - �񣨲��̶���ճ��ֱ�Ӳ�������ģ�ͼ��㣩��1 - �� ���̶���ճ���м��㣩��3 - ���� Beta ����м���
STD_Comp=[1,1,1]; % �Ƿ���е�ǰ�������׼�����������ıȽϣ����������� RU��U��V��P��TMu �ľ��������Լ���MSE��RMSE��
                  % ��Ҫ��Output�ļ����е��� FlowField_STD.PLT �ļ������ FlowField_Err.PLT �� FlowField_Err_Info.txt �ļ�
L2NORM_Setting=[0,0,0]; % ����ÿ�� L2NORM ����Ҫ��ˮƽ���ﵽ���˳���ǰ���㣨ע��5 - ���������ض�ˮƽΪ 5%�����ҽ����������㿪�� JX_Turb = 1 ʱ��Ч��


% �����������ͣ�0 - �����������赼�������������������Լ�ʵ�����ݣ���4�������ļ��������в�ֵ����
%               1 - �������������̨�ף�Backstep�����շ壨Hump���ȣ��赼��������β����Լ�ʵ�����ݣ���2�������ļ��������в�ֵ����
flag_flow_type=0;

flag_std=0;  % �����Ƿ���г�ʼ��׼������0 - ��1 - ��
flag_da_part=0;  % �����Ƿ���в�����������ͬ����0 - ��1 - ��
flag_da_method=2; % ��������ͬ�������㷨��1 - ���Ͽ������˲���EnKF����2 - ���ϱ任�������˲���ETKF��
flag_da_relative=1; % ��������ͬ���۲�ֵ���ͣ�0 - �����ͬ��Ŀ���Ϊʹ�������Ϊ�㣻1 - ��������ͳ��ʽ��
flag_beta_turb_type=0; % ��������ȡ������򣨶��� Beta������ ͬ�����Σ���0 - ÿ���Ŷ���̶���1 - ÿ���Ŷ��������

% �����ʼ�Ŷ���ͬ��Ŀ���г��� "Beta" �������Ŷ�����
if ((strcmp(char(S_Va_target_Name(1)),'Beta')==1) || (strcmp(char(S_Va_target_Name(1)),'Beta')==1))
    flag_da_part=1; % ������������ͬ��
    flag_da_beta=1; % �����Ŷ���Beta(x)ͬ��
    flag_Va_init=1; % �������ϳ�Ա����������ά����������־
    flag_sel_dot=2; % �����Ŷ�Beta(x)�ķ�ʽ
                    % �ⶨ���� 1 - ֱ��ѡȡ��ر�ż��������ѡȡ������ÿ��һ����Ž��е��������ѡȡɢ�㲢����ƽ���ֵ��
                    % �ⶨ���� 2 - ��ϸ���������Ȳ�����
                    % �ⶨ���� 3 - ������ά�⻬�������в�����������......����
    % ���� Beta �Ŷ����򣺷��� sigma_2������ֲ� Beta ~ N(std_beta,sigma_2)
    % �����ȡ������Ϣ
    DA_NG_gap=200;
    % �����~��ȡ�������Ϣ
    std_beta=1.0;
    sigma_beta=3.0;
    dis_ratio=5.0;
else
    % �������� Beta�Ŷ��� ͬ��
    flag_da_beta=0; % �ر��Ŷ���Beta(x)ͬ��
    flag_Va_init=0; % �رռ��ϳ�Ա����������ά����������ʹ����ͨ����������Ŷ�
end

% ��������ͬ������
NNDA=4; % NNDA -- ͬ������ڵ��������Ϊ�������Σ�������>=3������������ NNDA �����꣨XDA,YDA��
% XDA=[1.05,1.2,1.2,1.05]; % XDA(NNDA) -- ͬ������ڵ������
% YDA=[0.1,0.1,0.0,0.0]; % YDA(NNDA) -- ͬ������ڵ�������
% XDA=[0.8,1.5,1.5,0.8]; % XDA(NNDA) -- ͬ������ڵ������
% YDA=[0.1,0.1,0.02,0.02]; % YDA(NNDA) -- ͬ������ڵ�������
XDA=[0.5499846528026755,1.3977168074273609,1.3977168074273609,0.5499846528026755]; % XDA(NNDA) -- ͬ������ڵ������
YDA=[0.11283900866252453,0.16714063248086963,0.005,0.005]; % YDA(NNDA) -- ͬ������ڵ�������

if (flag_flow_type==0)
    % ��������������������ݣ��Զ���ֵ��
    shape_upper=load('.\Input_Data\Input_Data_Shape_RAE2822_upper.dat');
    shape_lower=load('.\Input_Data\Input_Data_Shape_RAE2822_lower.dat');
    % �������������ʵ������
    exp_upper_data=load('.\Input_Data\Input_Data_Exp_RAE2822_upper.dat');
    exp_lower_data=load('.\Input_Data\Input_Data_Exp_RAE2822_lower.dat');
elseif (flag_flow_type==1)
    % �����������ݣ��Զ���ֵ��
    shape_outline=load('.\Input_Data\Input_Data_Shape_Hump.dat');
    % ����ʵ������
    exp_outline=load('.\Input_Data\Input_Data_Exp_Hump.dat');
end
    
% ------------------  ������Ҫ���������� END ------------------ %

% -----------  ��������ͬ��ʵ��ֵ�������Ʋ��� START ----------- %

% ���ù۲�ֵЭ�����������ٶ���Щ���������
% !!!!!!! �����������ֵ��Ҫ�ٿ��ǣ���һЩ������Ӧ�������ģ�ͣ��� !!!!!!!
% 2021/04/27 ����ģ��
% �����Ƿ��� ����Ӧ�۲�ֵ��׼�����ģ�� - Adaptive observation standard deviation adjustment model 
%               abbr. AOSDAM���� - 0���̶�����ģʽ - 1������ƽ��ģʽ - 2��������Ӧģʽ - 3��
% �㷨˵�������ر�����Ӧģ�ͣ����趨���й۲�ֵ��׼��Ϊһ�̶�ֵ sigma_obs��
%           ����������Ӧģ�ͣ�1 - �̶�����ģʽ������ÿ��ͬ��ǰ���й۲�ֵ��׼������Ϊ��һ�ֱ�׼����Է������ӣ�
%           ����������Ӧģ�ͣ�2 - ����ƽ��ģʽ������ÿ��ͬ��ǰ������һ�����й۲�ֵ��Ӧ�ļ���ֵ��������ƽ����ȣ������й۲�ֵ�ı�׼�����Ϊƽ��������׼��ˮƽ�����Է������ӣ���
%           ����������Ӧģ�ͣ�3 - ������Ӧģʽ������ÿ��ͬ��ǰ������һ�ָ��۲�ֵ��Ӧ�ļ���ֵ���������ȣ������۲�ֵ�ı�׼�����Ϊ��Ӧ������׼��ˮƽ�����Է������ӣ���
% ����˵����flag_sigma_obs_adapt - ģʽ����
%           flag_sigma_obs_coeff - ��������
%           sigma_obs_init       - �趨��ʼ�۲�ֵ��׼��ˮƽ
%           sigma_obs_min        - �趨�۲�ֵ��׼�����ޣ�������ǰ��ı�׼������ɵ�������
%           sigma_obs_max        - �趨�۲�ֵ��׼�����ޣ�������ǰ��ı�׼������ɸ�������
flag_sigma_obs_adapt=0;
flag_sigma_obs_coeff=0.80;
sigma_obs_ini=0.03; % ���ó�ʼ�۲��׼��ˮƽ
sigma_obs_min=0.02; % ָ���۲��׼������Ϊ0.05����Ӧ����Ϊ0.0025
sigma_obs_max=0.80; % ָ���۲��׼������Ϊ0.80����Ӧ����Ϊ0.6400

% ��������ͬ����ʼ���Ϸ���ֵ���Э����
% C_pp_e_a=cov(S_Va_ana');

% ----------  ��������ͬ��ʵ��ֵ�������Ʋ��� END ---------- %

% ---------------  ���ò�������ͬ������ START --------------- %

% ��������ģ��ͬ������
if (flag_da_part==1)
    
    % �����ò���д�� InputDAArea.dat �ļ�
    % ָ���滻�У��ֱ��Ӧ�����ļ���1��2��3��
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
    
    % ��������������ͬ�������� Input �ļ��еĶ�Ӧ���� JX_DA = 1
    % ָ���滻�У���Ӧ�����ļ���InputNS.dat����21��
    rep_line_JX_DA=21;
    rep_format_JX_DA='%-12d';
    % ���£��رգ���������ͬ������
    JX_DA=2;
    fid_DA=fopen('.\Program_CFD_Base\Input\InputNS.dat','r+');
    for i=1:(rep_line_JX_DA-1)
       fgetl(fid_DA);
    end
    fseek(fid_DA,0,'cof');
    fprintf(fid_DA,rep_format_JX_DA,JX_DA);
    fclose(fid_DA);
    
    % �����в�����������ͬ��������ƣ��ǣ�����������ָ��ͬ������ʾ��ͼ
    if (flag_flow_type==0)
        DA_Extension_Airfoil_Interpolation(NNDA,XDA,YDA);
    elseif (flag_flow_type==1)
        DA_Extension_NonAirfoil_Interpolation(NNDA,XDA,YDA);
    end
    
end

% ----------------  ���ò�������ͬ������ END ---------------- %

% ---------------  ���������Ҫ������Ϣ START --------------- %

% �����������ʾ
fprintf('******************************************************************\n');
fprintf('                      ����ͬ��������Ҫ������Ϣ                    \n');
fprintf('******************************************************************\n');
fprintf('[������������]\n');
fprintf('���������������ͣ�0 - ����������1 - �̱������� - %d\n',flag_flow_type);
fprintf('[���򱣴�·�� & ��������]\n');
fprintf(['�������·����',savename,'\n']);
fprintf(['������Ϣ�ļ���DA_CompInfo_Brief.dat\n']);
fprintf('�������߳�/������ PPNS = %d\n',PPNS);
fprintf('******************************************************************\n');
fprintf('[ͬ����Ҫ����]\n');
fprintf('��󼯺ϳ�Ա���� N_max = %d\n',N_max);
fprintf('ͬ���۲���� k = %d, ����ֲ���鿴 num_da_arr ����\n',k);
fprintf('******************************************************************\n');
fprintf('[��ʼ�Ŷ�������Ϣ]\n');
if (flag_Va_init==0)
    for i=1:S_Va_turb_Num
        fprintf(['��ʼ�Ŷ����� ',num2str(i,'%03d'),'��',S_Va_turb_CharName(:,:,i),' �� [ %.5f , %.5f ], STD = %.5f\n'],S_Va_turb_min(i),S_Va_turb_max(i),S_Va_turb_Std(i));
    end
else
    for i=1:S_Va_turb_Num
        fprintf(['��ʼ�Ŷ��ֲ� ',num2str(i,'%03d'),'��',S_Va_turb_CharName(:,:,i),' \n']);
    end
end
fprintf('******************************************************************\n');
fprintf('[Ŀ��ͬ��������Ϣ]\n');
for i=1:S_Va_target_Num
    fprintf(['Ŀ��ͬ������ ',num2str(i,'%03d'),'��',S_Va_target_CharName(:,:,i),' ��ʼ�Ŷ����� - %d, ��/�ֵࣨ��ͬ������ - %d\n'],S_Va_target_TypeTurb(i),S_Va_target_TypeDA(i));
end
fprintf('******************************************************************\n');
fprintf('[ָ��ͬ���۲������Ϣ]\n');
for i=1:S_Va_ref_Num
    fprintf(['ָ��ͬ���۲���� ',num2str(i,'%03d'),'��',S_Va_ref_CharName(:,:,i),'\n']);
end
fprintf('******************************************************************\n');
fprintf('[ͬ�����̿���]\n');
fprintf('���ü������ɹ���������������%d, CFL = %d\n',MStepNS_cov,CFL_Setting(1));
fprintf('���ü����ڵ���Ԥ����㲽����%d, CFL = %d\n',MStepNS_pre,CFL_Setting(2));
fprintf('���ü��ϳ�Ա�ؼ��㲽����%d, CFL = %d\n',MStepNS_recal_ensemble,CFL_Setting(3));
fprintf('���ü����ڵ��������������ִΣ�%d\n',Max_inner_iter);
fprintf('�Ƿ���г�ʼ��׼���� - %d\n',flag_std);
fprintf('�Ƿ������������ͬ�� - %d\n',flag_da_part);
fprintf('��������ͬ�������㷨��1 - EnKF��2 - ETKF�� - %d\n',flag_da_method);
fprintf('��������ͬ���۲�ֵ���ͣ�0 - �����1 - ������ - %d\n',flag_da_relative);
fprintf('�Ƿ��������Դ�������� Beta ͬ�� - %d\n',flag_da_beta);
fprintf('�Ƿ��ڼ������ɹ��̵���������ά�������� - %d\n',flag_Va_init);
fprintf('�Ƿ��� ����Ӧ�۲�ֵ��׼�����ģ�ͣ�AOSDAM�� ѡ��ģʽ - %d\n',flag_sigma_obs_adapt);
fprintf('******************************************************************\n');
fprintf('\n');

% �ļ������ʾ
fprintf(fid_info,'******************************************************************\n');
fprintf(fid_info,'                      ����ͬ��������Ҫ������Ϣ                    \n');
fprintf(fid_info,'******************************************************************\n');
fprintf(fid_info,'[������������]\n');
fprintf(fid_info,'���������������ͣ�0 - ����������1 - �̱������� - %d\n',flag_flow_type);
fprintf(fid_info,'[���򱣴�·�� & ��������]\n');
fprintf(fid_info,['�������·����',savename,'\n']);
fprintf(fid_info,['������Ϣ�ļ���DA_CompInfo_Brief.dat\n']);
fprintf(fid_info,'�������߳�/������ PPNS = %d\n',PPNS);
fprintf(fid_info,'******************************************************************\n');
fprintf(fid_info,'[ͬ����Ҫ����]\n');
fprintf(fid_info,'��󼯺ϳ�Ա���� N_max = %d\n',N_max);
fprintf(fid_info,'ͬ���۲���� k = %d, ����ֲ���鿴 num_da_arr ����\n',k);
fprintf(fid_info,'******************************************************************\n');
fprintf(fid_info,'[��ʼ�Ŷ�������Ϣ]\n');
if (flag_Va_init==0)
    for i=1:S_Va_turb_Num
        fprintf(fid_info,['��ʼ�Ŷ����� ',num2str(i,'%03d'),'��',S_Va_turb_CharName(:,:,i),' �� [ %.5f , %.5f ], STD = %.5f\n'],S_Va_turb_min(i),S_Va_turb_max(i),S_Va_turb_Std(i));
    end
else
    for i=1:S_Va_turb_Num
        fprintf(fid_info,['��ʼ�Ŷ��ֲ� ',num2str(i,'%03d'),'��',S_Va_turb_CharName(:,:,i),' \n']);
    end
end
fprintf(fid_info,'******************************************************************\n');
fprintf(fid_info,'[Ŀ��ͬ��������Ϣ]\n');
for i=1:S_Va_target_Num
    fprintf(fid_info,['Ŀ��ͬ������ ',num2str(i,'%03d'),'��',S_Va_target_CharName(:,:,i),' ��ʼ�Ŷ����� - %d, ��/�ֵࣨ��ͬ������ - %d\n'],S_Va_target_TypeTurb(i),S_Va_target_TypeDA(i));
end
fprintf(fid_info,'******************************************************************\n');
fprintf(fid_info,'[ָ��ͬ���۲������Ϣ]\n');
for i=1:S_Va_ref_Num
    fprintf(fid_info,['ָ��ͬ���۲���� ',num2str(i,'%03d'),'��',S_Va_ref_CharName(:,:,i),'\n']);
end
fprintf(fid_info,'******************************************************************\n');
fprintf(fid_info,'[ͬ�����̿���]\n');
fprintf(fid_info,'���ü������ɹ���������������%d, CFL = %d\n',MStepNS_cov,CFL_Setting(1));
fprintf(fid_info,'���ü����ڵ���Ԥ����㲽����%d, CFL = %d\n',MStepNS_pre,CFL_Setting(2));
fprintf(fid_info,'���ü��ϳ�Ա�ؼ��㲽����%d, CFL = %d\n',MStepNS_recal_ensemble,CFL_Setting(3));
fprintf(fid_info,'���ü����ڵ��������������ִΣ�%d\n',Max_inner_iter);
fprintf(fid_info,'�Ƿ���г�ʼ��׼���� - %d\n',flag_std);
fprintf(fid_info,'�Ƿ������������ͬ�� - %d\n',flag_da_part);
fprintf(fid_info,'��������ͬ�������㷨��1 - EnKF��2 - ETKF�� - %d\n',flag_da_method);
fprintf(fid_info,'��������ͬ���۲�ֵ���ͣ�0 - �����1 - ������ - %d\n',flag_da_relative);
fprintf(fid_info,'�Ƿ��������Դ�������� Beta ͬ�� - %d\n',flag_da_beta);
fprintf(fid_info,'�Ƿ��ڼ������ɹ��̵���������ά�������� - %d\n',flag_Va_init);
fprintf(fid_info,'�Ƿ��� ����Ӧ�۲�ֵ��׼�����ģ�ͣ�AOSDAM�� ѡ��ģʽ - %d\n',flag_sigma_obs_adapt);
fprintf(fid_info,'******************************************************************\n');
fprintf(fid_info,'\n');

% ----------------  ���������Ҫ������Ϣ END ---------------- %

% -----------------  ����ִ�б�׼���� START ----------------- %

% ����һ�α�׼������Program_CFD_Base��
if (flag_std==1)
    fprintf('��׼���������� ��ʼ����\n');
    fprintf(fid_info,'��׼���������� ��ʼ����\n');
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
    fprintf('��׼���������� ���н���\n');
    fprintf(fid_info,'��׼���������� ���н���\n');
    fprintf('\n');
    fprintf(fid_info,'\n');
else
    fprintf('�Ѵ��ڱ�׼����������\n');
    fprintf(fid_info,'�Ѵ��ڱ�׼����������\n');
    fprintf('\n');
    fprintf(fid_info,'\n');
end

% ------------------  ����ִ�б�׼���� END ------------------ %

%% ����һ������ο�ʵ�������ּ������ݣ������в�ֵ����
% �������ݲ���
% ��������ԭʼ�������ݣ��������������Ͳ�ֵ
% �����������
com_data=load('.\Program_CFD_Base\Output\CP_CF_Comp.PLT');
% ָ��������������ݱ���·��
fid_com_adj_path='.\Program_CFD_Base\Output\CP_CF_Comp_New.PLT';
fid_tar_new_path='.\Program_CFD_Base\Output\CP_New.PLT';

% ��������ʵ��ֵ����
if (flag_flow_type==0)
    
    % �Լ���ѹ��ϵ�����߽ڵ��Ž���������ʵʱ�����ǰͬ���ִ�ѹ���ֲ�����
    % �������ɹ��̡��ڵ���������ȫ����Ա��ɼ���󱣴����߰�����
    [shape_upper_num,shape_lower_num,shape_upper_arr_sor,shape_lower_arr_sor,n_com]=DA_Va_Airfoil_Ref_Adjustment_and_Visulization(com_data,fid_com_adj_path,fid_tar_new_path);

    % ʵ�����ݲ���
    % ����������ʵ�������������x���������������ָ��[0,1]���䣩
    exp_upper_data_sor=sortrows(exp_upper_data,1);
    exp_lower_data_sor=sortrows(exp_lower_data,1);
    % д������ļ����µ� InputExpData_all.dat �� InputExpData.PLT �ļ�
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
    % ��¼ʵ�������
    exp_upper_num=length(exp_upper_data_sor);
    exp_lower_num=length(exp_lower_data_sor);
    exp_num_tot=exp_upper_num+exp_lower_num;
    % ������������ʵ������ֵ�Ĳ�ֵ
    xx_upper_exp=(0:0.001:1);
    yy_upper_exp=interp1(exp_upper_data_sor(:,1),exp_upper_data_sor(:,2),xx_upper_exp,'linear','extrap');
    xx_lower_exp=(0:0.001:1);
    yy_lower_exp=interp1(exp_lower_data_sor(:,1),exp_lower_data_sor(:,2),xx_lower_exp,'linear','extrap');

    % ���ò�ֵ��������������ݵ������Ӧ��ʵ��ֵ
    for i=1:shape_upper_num
        x_com_upper_pre=shape_upper_arr_sor(i,2);
        shape_upper_arr_sor(i,7)=interp1(exp_upper_data_sor(:,1),exp_upper_data_sor(:,2),x_com_upper_pre,'linear','extrap');
    end
    for i=1:shape_lower_num
        x_com_lower_pre=shape_lower_arr_sor(i,2);
        shape_lower_arr_sor(i,7)=interp1(exp_lower_data_sor(:,1),exp_lower_data_sor(:,2),x_com_lower_pre,'linear','extrap');
    end

    % �ϲ��������������ʵ��ֵ����ֵ��������
    shape_combined_arr_init=[shape_lower_arr_sor;shape_upper_arr_sor]; % �������棬��������

    % % ���������ݵ�ԭʼ��Ž�����������
    % shape_combined_arr_sor=sortrows(shape_combined_arr,1);
    
elseif (flag_flow_type==1)
    
    % �Լ���ѹ��ϵ�����߽ڵ��Ž���������ʵʱ�����ǰͬ���ִ�ѹ���ֲ�����
    % �������ɹ��̡��ڵ���������ȫ����Ա��ɼ���󱣴����߰�����
    [shape_num,shape_arr_sor,n_com]=DA_Va_NonAirfoil_Ref_Adjustment_and_Visulization(com_data,fid_com_adj_path,fid_tar_new_path);
    % ʵ�����ݲ���
    % ��ʵ�������������x���������������ָ��[0,1]���䣩
    exp_data_sor=sortrows(exp_outline,1);
    % д������ļ����µ� InputExpData_all.dat �� InputExpData.PLT �ļ�
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
    % ��¼ʵ�������
    exp_num_tot=length(exp_data_sor);
    % ������������ʵ������ֵ�Ĳ�ֵ
    xx_lower_exp=(0:0.001:1);
    yy_lower_exp=interp1(exp_data_sor(:,1),exp_data_sor(:,2),xx_lower_exp,'linear','extrap');
    % ���ò�ֵ��������������ݵ������Ӧ��ʵ��ֵ
    for i=1:shape_num
        x_com_pre=shape_arr_sor(i,2);
        shape_arr_sor(i,7)=interp1(exp_data_sor(:,1),exp_data_sor(:,2),x_com_pre,'linear','extrap');
    end
    
    shape_combined_arr_init=shape_arr_sor;

end

% ������ʵ��ֵ��������ӳ���ϵ
% ��¼����ͬ����ʵ�����ݵ���Ϣ���ɵ����ļ�
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

% �������ͬ���Ĺ۲�ʵ��ֵ�ֲ�
fid_exp_da=fopen('.\Program_CFD_Base\Output\InputExpData_for_DA.PLT','w+');
fprintf(fid_exp_da,'VARIABLES="X","CP"\n');
fprintf(fid_exp_da,'ZONE T="ZONE EXP for DA"\n');
for i=1:k
    fprintf(fid_exp_da,'%20.8f %20.8f\n',exp_for_DA(i,1),exp_for_DA(i,2));
end
fclose(fid_exp_da);

% ���Ʊ�׼����ģ�ͼ�����������ģ�廯��ʾ
copyfile('.\Program_CFD_Base\Output\CP_New.PLT','.\Program_CFD_Base\Output\CP_SA.PLT');
% ���Ʊ�׼����ģ�ͼ��������������Ƚϣ������������Է�����
copyfile('.\Program_CFD_Base\Output\FlowField.PLT','.\Program_CFD_Base\Output\FlowField_Ori.PLT');
copyfile('.\Program_CFD_Base\Output\FlowField_Info.PLT','.\Program_CFD_Base\Output\FlowField_STD.PLT');

% ��ʼ��׼״̬���������ӻ�
figure;
hold on
plot(shape_combined_arr_init(:,2),shape_combined_arr_init(:,7),'ro','MarkerSize',4.0);
plot(shape_combined_arr_init(:,2),shape_combined_arr_init(:,4),'b-','LineWidth',2.0);
xlabel('\itx'),ylabel('\itCp'),title('��ʼ��׼״̬������');
legend('ʵ�����ֵ','��ʼ����ѹ��ϵ���ֲ�');
grid on
box on
hold off
% ����ɱ༭ͼƬ
saveas(gcf,['.\Output_Figures\',savename,'_Init_State_Output.fig']);

%% ����������г����鴴��

% ��ȡ�����������
% ��������ת��������
fprintf('������������ͬ������������\n');
fprintf(fid_info,'������������ͬ������������\n');
fprintf('\n');
fprintf(fid_info,'\n');
cmd_grid=('.\Program_CFD_Base\GridTreat_.exe');
open(cmd_grid);
pause(3);
fprintf('������������ת����ɣ�\n');
fprintf(fid_info,'������������ת����ɣ�\n');

% �������������ļ��У�����ΪPPNS��
for po=1:PPNS
    ddir=['Program_CFD_Parallel_',num2str(po,'%03d')];
    copyfile('Program_CFD_Base',ddir);
    % �״���������ת��
    cmd_grid=([ddir,'\GridTreat_.exe']);
    open(cmd_grid);
    pause(3);
    fprintf('���г����� %03d ����ת����ɣ�\n',po);
    fprintf(fid_info,'���г����� %03d ����ת����ɣ�\n',po);
end

% �����ؼ�����̳�����
copyfile('Program_CFD_Base',ddir_recal);
cmd_grid=([ddir_recal,'\GridTreat_.exe']);
open(cmd_grid);
pause(3);
fprintf('�ؼ�������� ����ת����ɣ�\n');
fprintf(fid_info,'�ؼ�������� ����ת����ɣ�\n');

% ����������������
copyfile('Program_CFD_Base',ddir_post);
cmd_grid=([ddir_post,'\GridTreat_.exe']);
open(cmd_grid);
pause(3);
fprintf('������������ ����ת����ɣ�\n');
fprintf(fid_info,'������������ ����ת����ɣ�\n');

fprintf('\n');
fprintf(fid_info,'\n');

%% ����������ʼ�����ϳ�Ա

% ---------------- ����������������ϵ������ START ---------------- %

% ����������ģ��Դ�������� Beta ͬ��
if (flag_da_beta==1)
    
    % ��ȡ��������ͬ����Ϣ
    DA_area_info=load('.\Program_CFD_Base\Input\InputDAgrid.dat');
    % ��һ�б�ʾͬ������������ڵ�
    NDAV=DA_area_info(1);
    % �ڶ��б�ʾͬ�������ڵ�Ԫ����
    NDAG=DA_area_info(2);
    % ��ȡͬ������������ڵ��뵥Ԫ�����Ϣ
    DA_VN=DA_area_info(3:(NDAV+2));
    DA_GN=DA_area_info((NDAV+3):(NDAV+NDAG+2));
    
    % ��ȡ��ǰ������������ڵ㡢��Ԫ����������Ϣ
    VN_grid=load('.\Program_CFD_Base\Input\Input_NNP_Coordinate.dat');
    GN_grid=load('.\Program_CFD_Base\Input\Input_NG_Coordinate.dat');
    % ��ȡ��������ڵ�����
    VN_grid_x=VN_grid(:,1);
    VN_grid_y=VN_grid(:,2);
    NNP_tot=length(VN_grid);
    % ��ȡ��������Ԫ��������
    GN_grid_x=GN_grid(:,1);
    GN_grid_y=GN_grid(:,2);
    NG_tot=length(GN_grid);
    
    % ��ȡѡȡ����������ڵ㡢��Ԫ����������Ϣ
    DA_area_coor_info=load('.\Program_CFD_Base\Input\InputDAgrid_Coordinate.dat');
    % ��ȡ����������ڵ����
    NDAV_coor=DA_area_coor_info(1,1);
    % ��ȡ����������Ԫ����
    NDAG_coor=DA_area_coor_info(1,2);
    % ��ȡѡȡ����������ڵ��뵥Ԫ�����Ϣ
    DA_VN_x=DA_area_coor_info((2:(NDAV_coor+1)),1); % ѡȡ����������ڵ� x ����
    DA_VN_y=DA_area_coor_info((2:(NDAV_coor+1)),2); % ѡȡ����������ڵ� y ����
    DA_GN_x=DA_area_coor_info(((NDAV_coor+2):(NDAV_coor+NDAG_coor+1)),1); % ѡȡ����������Ԫ x ����
    DA_GN_y=DA_area_coor_info(((NDAV_coor+2):(NDAV_coor+NDAG_coor+1)),2); % ѡȡ����������Ԫ x ����
    
    % ��ǳ�����ͬ�������ڵ�����Ԫ��ţ��� - 1���� - 0
    DA_GN_flag_tot=zeros(NG_tot,1);
    for i=1:NDAG
        DA_GN_flag_tot(DA_GN(i))=1;
    end
    
    % ��ʼ�� ������Beta ����
    Va_Beta=zeros(NG_tot,N_max);
    Va_NNP_Beta=zeros(NNP_tot,N_max);
    
    % ��������ȡ�������ÿ���Ŷ������ or �̶�����
    %     if (flag_sel_dot==1)
    %         % �ⶨ����٣�ֱ��ѡȡ��ر�ż��������ѡȡ������ÿ��һ����Ž��е��������ѡȡɢ��
    %         DA_NG_gap=200;
    %     elseif (flag_sel_dot==2)
    %         % �ⶨ����ڣ�Բ������
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
    
    % ��ȡ�߽���һ������Ԫ�������ͷ���ʱ�ڴ�
    % Ŀ�ģ���ѡ������֮�⣬�Լ��߽���һ������Ԫ�� Beta�Ŷ�ֵ
    % ����Ϊ��׼ֵ 1.0���ң�����У����߽���һ������Ԫ�ų���ѡ������֮��
    Grid_info=textread('.\Program_CFD_Base\Input\US2D.dat');
    BD_num=Grid_info(1,1); % ��ȡ�߽���һ�������Ԫ����
    clear Grid_info; % ���������Ϣ
    
    for pp=1:N_max

        % �����ⶨ����ѡȡ����ɢ�㣬�������Ŷ�Betaֵ�ĸ�ֵ
        XY_NG_Beta=ones(NG_tot,1); % ��ʼ�� Beta ����
        XY_NNP_Beta=ones(NNP_tot,1); % ��ʼ�� Beta ����
        
        if (flag_sel_dot==1)
            % �ⶨ����٣�ֱ��ѡȡ��ر�ż��������ѡȡ������ÿ��һ����Ž��е��������ѡȡɢ��
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
            % ����MATLAB�Դ��� scatteredInterpolant() ������ѡȡ��ɢ����������ֵ������ V = F(x,y) ���������
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
            
            % �ⶨ����ڣ���ϸ������ + �Ŷ�
            % �����Ŷ���������
            XY_dots_num_max=round(0.1*NDAG_coor);
            XY_dots_dis_lim=dis_ratio*max(((max(XDA)-min(XDA))/XY_dots_num_max),((max(YDA)-min(YDA))/XY_dots_num_max));
            
            % �� flag_beta_turb_type==1���������Ŷ�����Ԫ���أ���ÿ�ξ���������� Beta��ʼ�Ŷ��㣬������� XY_cnt �� XY_sel_flag ֵ
            % �� flag_beta_turb_type==0����ر�����Ŷ�����Ԫ���أ��������һ�� Beta��ʼ�Ŷ��㣬֮��̶�
            if ((flag_beta_turb_type==1) || ((flag_beta_turb_type==0) && (pp==1)))
                [XY_cnt_pre,Xp_sca,Yp_sca,XY_sel_flag_pre]=Region_Beta_Dots_Selected(NDAG_coor,DA_GN_x,DA_GN_y,XY_dots_num_max,XY_dots_dis_lim);
            end
            
            XY_sel_flag_real=zeros(XY_cnt_pre,1);
            for i=1:XY_cnt_pre
                XY_sel_flag_real(i)=DA_GN(XY_sel_flag_pre(i));
            end
            
            % ���� Beta �Ŷ�ֵ��ʼ��
            % ��ȫ�� Beta ֵ����Ϊ��׼ֵ 1.0���̶��߽�ֵ��
            Vp_NG_std_all=ones(NG_tot,1);
            % ��ָ������Ԫ���� Beta �Ŷ�
            for i=1:XY_cnt_pre
                Vp_NG_std_all(XY_sel_flag_real(i))=normrnd(std_beta,sigma_beta);
            end
            % ���߽���һ�㵥ԪBetaֵ����Ϊ��׼ֵ 1.0 ������
            for i=1:BD_num
                Vp_NG_std_all(i)=1.0;
            end
            % ��ȡ��ƽ���ֵ������Ϣ������ȡ���Ŷ��㣩
            cnt_non_fit=0;
            for i=1:NG_tot
                if ((i<=BD_num) || (DA_GN_flag_tot(i)==0) || (sum(find(XY_sel_flag_real==i))~=0))
                    cnt_non_fit=cnt_non_fit+1;
                    Xp_non_all(cnt_non_fit)=GN_grid_x(i);
                    Yp_non_all(cnt_non_fit)=GN_grid_y(i);
                    Vp_non_all(cnt_non_fit)=Vp_NG_std_all(i);
                end
            end
            
            % ����MATLAB�Դ��� scatteredInterpolant() ������ѡȡ��ɢ����������ֵ������ V = F(x,y) ���������
            F=scatteredInterpolant(Xp_non_all',Yp_non_all',Vp_non_all');
            F.Method='natural';
            F.ExtrapolationMethod='nearest';
            
            % ��ֵ�õ���������Ԫ Beta ֵ
            Vp_NG_fit=F(DA_GN_x,DA_GN_y); % ָ������Ԫ��ֵ
            Xp_NG_all=GN_grid_x;
            Yp_NG_all=GN_grid_y;
            XY_NG_Beta=F(Xp_NG_all,Yp_NG_all); % ָ������Ԫ��ֵ
            
            % ������������Ԫ Beta ֵ��ָ�������ڵĵ�Ԫ�ڵ���в�ֵ������ϵõ�ȫ������ڵ�� Beta ֵ��Ȼ��ʵ�ֿ��ӻ�
            % ��ȫ�� Beta ֵ��ȫ������ڵ㣩����Ϊ��׼ֵ 1.0
            % F=scatteredInterpolant(Xp_non_all',Yp_non_all',XY_NG_Beta');
            Vp_NNP_std_all=ones(NNP_tot,1);
            Vp_NNP_fit=F(DA_VN_x,DA_VN_y); % ָ����������ڵ��ֵ
            Xp_NNP_all=VN_grid_x;
            Yp_NNP_all=VN_grid_y;
            XY_NNP_Beta=F(Xp_NNP_all,Yp_NNP_all); % ��ֵ�õ���������ڵ� Beta ֵ����������Ԫֵ��ֵ��
            
        end
        
        % ��ɳ�ʼ Beta �ֲ�������������
        Va_Beta(:,pp)=XY_NG_Beta;
        Va_NNP_Beta(:,pp)=XY_NNP_Beta;
        
        % �����ڴ�
        clear Xp_non_all;
        clear Yp_non_all;
        clear Vp_non_all;
    
    end
    
    % ��ѡ�Ŷ�������ӻ�����������Ԫ Beta ���ݣ�
    figure;    
    scatter(Xp_sca,Yp_sca);
    xlabel('x'),ylabel('y');
    % ����ɱ༭ͼƬ
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
    % ����ɱ༭ͼƬ
    saveas(gcf,['.\Output_Figures\',savename,'\Init_Final_Beta_SampleContour.fig']);
    
end

% ---------------- ����������������ϵ������ END ---------------- %

% -------------- �޸Ĳ���Ԥ�� + ���ϳ�Ա���� START -------------- %

% [S_Va_init_Combined,NNP,NG,MDA,M]=State_Ensem_Gen_Parellel_Common(S_Va_init,N_max,MStepNS_cov,targetPath);
[S_Va_init_Combined_eff,target_width,Val_Sen,N_new,NNP,NG,MDA,M]=State_Ensem_Gen_Parellel_Common_Beta_BL(S_Va_init,N_max,MStepNS_cov,targetPath); % ����֧�� �Ŷ���Beta(x) ��ͬ��

fprintf('���ϳ�Ա���ɹ��� ��������󼯺ϳ�Ա�� = %d����Ч���ϳ�Ա�� = %d��������ָ�� TMu_re = %.6f\n',N_max,N_new,Val_Sen);
fprintf('\n');
fprintf(fid_info,'���ϳ�Ա���ɹ��� ��������󼯺ϳ�Ա�� = %d����Ч���ϳ�Ա�� = %d��������ָ�� TMu_re = %.6f\n',N_max,N_new,Val_Sen);
fprintf(fid_info,'\n');

% ����ȫ״̬�������󣬳�ʼ���ϳ�Ա����
S_Va_ana=S_Va_init_Combined_eff;

% --------------- �޸Ĳ���Ԥ�� + ���ϳ�Ա���� END --------------- %

% ----------------- ��������ͬ����Ϣ���� START ----------------- %

% ��������������ͬ��
if (flag_da_part==1)

    % �����壺���涨ͬ������������Ԫ��Ӧ����ճϵ������ͬ����
    %           ͬ���������ʣ������Ԫ��ճϵ������ԭֵ��

    % ��ȡ��������ͬ����Ϣ
    DA_area_info=load('.\Program_CFD_Base\Input\InputDAgrid.dat');
    % ��һ�б�ʾͬ������������ڵ�
    NDAV=DA_area_info(1);
    % �ڶ��б�ʾͬ�������ڵ�Ԫ����
    NDAG=DA_area_info(2);
    % ��ȡͬ������������ڵ��뵥Ԫ�����Ϣ
    DA_VN=DA_area_info(3:(NDAV+2));
    DA_GN=DA_area_info((NDAV+3):(NDAV+NDAG+2));
    
    % ��ǳ�����ͬ�������ڵ�����Ԫ���
    DA_GN_flag=zeros(NG,1);
    for i=1:NDAG
        DA_GN_flag(DA_GN(i))=1;
    end

    % ����ָ��ͬ��������� MDA_new �� M_new ֵ����С�ڻ����ԭʼֵ
    MDA_new=NDAG+PT;
    M_new=MDA_new+n_com;
else
    MDA_new=MDA;
    M_new=M;
end

% ------------------ ��������ͬ����Ϣ���� END ------------------ %

% �����ʼ������Ϣ
% ����ɹ�����ֱ�ӵ��û�����Ĺ���������
save([new_output_folder,savename,'_Pre.mat']);  % ���漯�ϳ�Ա���ɹ������н��
save([new_output_folder,savename,'_S_Va_init.mat'],'S_Va_ana'); % ������³�Ա�����������
save([new_output_folder,'S_Va_init.dat'],'S_Va_ana','-ascii');
save([new_output_folder,'S_Va_init.mat'],'S_Va_ana');

%% �����ģ�����ʵ���������������ӳ���ϵ
if ((strcmp(char(S_Va_ref_Name(1)),'Cp')==1) || (strcmp(char(S_Va_ref_Name(1)),'Cf')==1) || (strcmp(char(S_Va_ref_Name(1)),'yplus')==1))
    
    % �����۲����Ӿ���H
    H=zeros(k,M_new);
    for i=1:k
        H(i,(M_new-n_com)+num_da_arr(i))=1;  % !!! ͨ�ó����У��۲����� H ���޸�
    end
    H_Cp=H(:,((M_new-n_com+1):M_new));

    % % !!! ��Ҫ�޸ģ�����λ�ø��� !!!
    % ����ԭʼ�۲����ݼ��۲���������
    D_obs_e_ori=zeros(k,N_max);
    Y_err_ori=normrnd(0,sigma_obs_ini,[k,N_max]);
    for i=1:k
        % D_obs_e(i,:)=(ones(1,N).*exp_data(i,2))+Y_err(i,:);
        % D_obs_e(i,:)=(ones(1,N).*shape_combined_arr_sor(num_da_arr(i),4))+Y_err(i,:);
        % D_obs_e_ori(i,:)=(ones(1,N_max).*shape_combined_arr_init(num_da_arr(i),7))+Y_err_ori(i,:);   % ��ʵ�۲�ֵ�����ò�ֵ�õ���
        D_obs_e_ori(i,:)=(ones(1,N_max).*shape_combined_arr_init(num_da_arr(i),7));   % ��ʵ�۲�ֵ�����ò�ֵ�õ���
    end

    % ����ģʽ���Э�����ʵ���пɼ����޹���������
    std_q=0.0;
    sigma_2_q=std_q^2;
    C_qq_e=ones(M_new,M_new)*sigma_2_q;
    
elseif (strcmp(char(S_Va_ref_Name(1)),'BL')==1)
    
    % !!!!! %
    n_com=kb;
    % num_da_arr=ones(1,kb);
    num_da_arr=(1:kb); % !!! ע��ԭʼ����BUG !!!
    
    % �����۲����Ӿ���H
    H=zeros(kb,M_new);
    for i=1:kb
        H(i,(M_new-n_com)+num_da_arr(i))=1;  % !!! ͨ�ó����У��۲����� H ���޸�
    end
    H_Cp=H(:,((M_new-n_com+1):M_new));

    % % !!! ��Ҫ�޸ģ�����λ�ø��� !!!
    % ����ԭʼ�۲����ݼ��۲���������
    D_obs_e_ori=zeros(kb,N_max);
    Y_err_ori=normrnd(0,sigma_obs_ini,[kb,N_max]);
    for i=1:kb
        % D_obs_e(i,:)=(ones(1,N).*exp_data(i,2))+Y_err(i,:);
        % D_obs_e(i,:)=(ones(1,N).*shape_combined_arr_sor(num_da_arr(i),4))+Y_err(i,:);
        % D_obs_e_ori(i,:)=(ones(1,N_max).*shape_combined_arr_init(num_da_arr(i),7))+Y_err_ori(i,:);   % ��ʵ�۲�ֵ�����ò�ֵ�õ���
        D_obs_e_ori(i,:)=ones(1,N_max).*Exp_Obs_Arr(i);   % ��ʵ�۲�ֵ�����ò�ֵ�õ���
    end

    % ����ģʽ���Э�����ʵ���пɼ����޹���������
    std_q=0.0;
    sigma_2_q=std_q^2;
    C_qq_e=ones(M_new,M_new)*sigma_2_q;
    
end

%% �����壺��ʼ�ڵ������̣�����ͬ���㷨ʵ��
% �ؼ��㣺ʵ�ֶ�CFD�������㹦�ܵĿ���
% Case 1�����ر�����ģʽ JX_NS = 0���� JX_Turb = 1��������� TMu(I) ����Ԫ����
%                                      JX_Turb = 0���򲻴����κ��Ѿ��������������Ԫ����
% Case 2������������ģʽ JX_NS = 1���� JX_Turb = 1���򽫴��� RU(I),U(I),V(I),P(I),TMu(I)���������������� ����Ԫ����
%                                      JX_Turb = 0���򽫴��� RU(I),U(I),V(I),P(I) ����Ԫ����
% �ڵ���ʵ�ֹ��̣�
% Step 1��Ԥ�ⲽ ���� �����ʼ��Ma,Alfa��TMu���ݣ�������� RU(I),U(I),V(I),P(I) �ȱ�����Ma,Alfa,TMu(I) ͬ����������
% Step 2���˲��� ���� �����˲�ͬ���㷨���������ͱ���ѹ���ֲ�ʵ��ֵ�����ֵ�Ĳ��죬����ͬ������ Ma,Alfa,TMu(I)
% Step 3������ ���� ���������������еĲ�̬���ݽ�������
% ע��ʵ������ʱ����������һ��������Ԥ�ⲽ����ÿ��һ��ʱ������˲����ͺ���ֱ��ͬ�������������L2Norm����

% ���ݽ���ģ���д
% д��������״̬��������Ϊ��һ�ֵ�������

% ----------------------- �ڵ������� START ------------------------ %

% p=0; % ���븳ֵ�������ʼ���׶��Ѷ��壻�����ֶ����ԣ�ȷ�������������������
while (p<=Max_inner_iter)

    % ------------- �޸Ĳ��� + ���ϳ�Ա���㣨Ԥ�ⲽ�� START ------------- %
    
    if (p>0)
        N_pre=N_new;
        [S_Va_pre_eff,target_width,Val_Sen,N_new]=State_Ensem_Renew_Parellel_Common_Beta_BL(p,S_Va_ana,N_pre,MDA,MStepNS_pre,targetPath);
        S_Va_pre=S_Va_pre_eff;
        % �������¼��ϳ�Ա��
        fprintf(['���ϳ�Ա�ڵ����� ',num2str(p,'%03d'),' �ֽ�������󼯺ϳ�Ա�� = %d����Ч���ϳ�Ա�� = %d��������ָ�� TMu_re = %.6f\n'],N_max,N_new,Val_Sen);
        fprintf('\n');
        fprintf(fid_info,['���ϳ�Ա�ڵ����� ',num2str(p,'%03d'),' �ֽ�������󼯺ϳ�Ա�� = %d����Ч���ϳ�Ա�� = %d��������ָ�� TMu_re = %.6f\n'],N_max,N_new,Val_Sen);
        fprintf(fid_info,'\n');
    else
        S_Va_pre=S_Va_ana;
    end
    
    % ------------ ���ɹ۲����ݼ�����Ӧģ�ͱ�׼����� START ------------- %
    
    % !!! ��Ҫ�޸ģ�����λ�ø��� !!!
    % ���ɹ۲����ݾ���
    if ((strcmp(char(S_Va_ref_Name(1)),'Cp')==1) || (strcmp(char(S_Va_ref_Name(1)),'Cf')==1) || (strcmp(char(S_Va_ref_Name(1)),'yplus')==1))
        D_obs_e=zeros(k,N_new);
        D_obs_e_z=zeros(k,N_new); % �µĹ۲���������Ϊ�����Ŀ��ͬ��Ϊ��С�۲�ֵ��ʵ��ֵ����ֵ��֮������������Ϊ�㣨��ҪBUG��
        for i=1:k
            D_obs_e(i,:)=(ones(1,N_new).*shape_combined_arr_init(num_da_arr(i),7));   % ��ʵ�۲�ֵ�����ò�ֵ�õ���
        end
    elseif (strcmp(char(S_Va_ref_Name(1)),'BL')==1)
        D_obs_e=zeros(kb,N_new);
        D_obs_e_z=zeros(kb,N_new); % �µĹ۲���������Ϊ�����Ŀ��ͬ��Ϊ��С�۲�ֵ��ʵ��ֵ����ֵ��֮������������Ϊ�㣨��ҪBUG��
        for i=1:kb
            D_obs_e(i,:)=ones(1,N_new).*Exp_Obs_Arr(i);   % ��ʵ�۲�ֵ�����ò�ֵ�õ���
        end       
    end
    
    if ((strcmp(char(S_Va_ref_Name(1)),'Cp')==1) || (strcmp(char(S_Va_ref_Name(1)),'Cf')==1) || (strcmp(char(S_Va_ref_Name(1)),'yplus')==1))
        % 2021/04/27 ֲ������Ӧ�۲�ֵ��׼��ģ�� - AOSDM�������ԣ�
        % ��ʼ��ģ�Ͳ���
        sigma_obs_sum=0;
        % Case 0���ر�ģ�Ϳ��أ��̶��۲�ֵ��׼��Ϊԭʼֵ
        if (flag_sigma_obs_adapt==0)
            sigma_obs_new=ones(k,1).*sigma_obs_ini;

        % Case 1����������Ӧģ�ͣ�1 - �̶�����ģʽ������ÿ��ͬ��ǰ���й۲�ֵ��׼������Ϊ��һ�ֱ�׼����Է�������
        elseif (flag_sigma_obs_adapt==1)
            if (p==0)
                sigma_obs_new=ones(k,1).*sigma_obs_ini;
                sigma_obs_val=sigma_obs_ini;
            else
                sigma_obs_val=flag_sigma_obs_coeff*sigma_obs_val;
                sigma_obs_new=ones(k,1).*sigma_obs_val;
            end
            % �������ˮƽ�߽�
            for i=1:k
                if (sigma_obs_new(i)<=sigma_obs_min)
                    sigma_obs_new(i)=sigma_obs_min;
                end
                if (sigma_obs_new(i)>=sigma_obs_max)
                    sigma_obs_new(i)=sigma_obs_max;
                end
            end

        % Case 2����������Ӧģ�ͣ�2 - ����ƽ��ģʽ������ÿ��ͬ��ǰ������һ�����й۲�ֵ��Ӧ�ļ���ֵ��������ƽ����ȣ������й۲�ֵ�ı�׼�����Ϊƽ��������׼��ˮƽ�����Է������ӣ�
        elseif (flag_sigma_obs_adapt==2)
            for i=1:k
                sigma_obs_sum=sigma_obs_sum+target_width(num_da_arr(i));
            end
            sigma_obs_ave=sigma_obs_sum/k;
            % �� 95% ���������Ӧ �� �� 1.96 �ң�����趨Ϊ 3.92 ��
            % ����õ��µĹ۲�ֵ��׼��ȡֵ
            sigma_obs_new=ones(k,1).*(flag_sigma_obs_coeff*(sigma_obs_ave/3.92));
            % �������ˮƽ�߽�
            for i=1:k
                if (sigma_obs_new(i)<=sigma_obs_min)
                    sigma_obs_new(i)=sigma_obs_min;
                end
                if (sigma_obs_new(i)>=sigma_obs_max)
                    sigma_obs_new(i)=sigma_obs_max;
                end
            end

        % Case 3������������Ӧģ�ͣ�3 - ������Ӧģʽ������ÿ��ͬ��ǰ������һ�ָ��۲�ֵ��Ӧ�ļ���ֵ���������ȣ������۲�ֵ�ı�׼�����Ϊ��Ӧ������׼��ˮƽ�����Է������ӣ�
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

        % ��������Ӧ�����������������ɹ۲�����������
        Y_err=zeros(k,N_new);
        for i=1:k
            for j=1:N_new
                Y_err(i,j)=normrnd(0,sigma_obs_new(i));
            end
        end
        % Y_err=normrnd(0,sigma_obs_init,[k,N_new]);
        C_ss_e=cov(Y_err'); % �۲����Э�������ά��Ϊ�۲������k*k
        
    elseif (strcmp(char(S_Va_ref_Name(1)),'BL')==1)
        
        sigma_obs_new=ones(kb,1).*sigma_obs_ini;
        Y_err=zeros(kb,N_new);
        for i=1:kb
            for j=1:N_new
                Y_err(i,j)=normrnd(0,sigma_obs_new(i));
            end
        end
        % Y_err=normrnd(0,sigma_obs_init,[k,N_new]);
        C_ss_e=cov(Y_err'); % �۲����Э�������ά��Ϊ�۲������k*k
        
    end
    
    % ------------- ���ɹ۲����ݼ�����Ӧģ�ͱ�׼����� END -------------- %
        
    % -- ���ݲ��� -- %
    save([new_output_folder,savename,'_S_Va_pre_restart_',num2str(p,'%03d'),'.mat']); % �����ֹĿǰ���� p+1 ��ͬ����ʼǰ��ȫ��������������������ֱ�ӵ��������ļ���ʼ����
    % -- ���ݲ��� -- %
    
    % ��ͬ��Ԥ����̣�Prior����õ���ճϵ��ȡ����������δ���ã�
    % for i=(PT+1):MDA
    %     for j=1:N
    %         S_Va_pre(i,j)=log10(S_Va_pre(i,j));
    %     end
    % end
    
    % ------------ ������������ͬ�� [ǰ����]״̬���� START -------------- %
    
    % ��������������ͬ��������״̬����
    if (flag_da_part==1)
    
        % ����ȫ����״̬���� S_Va_pre
        % !!!!! �޸� ��Ӧͨ�ó��� !!!!!
        S_Va_pre_NDA=S_Va_pre; % �����������״̬�����з�ͬ����������
        S_Va_pre_FDA=S_Va_pre; % �����������״̬��������ͬ����������
        % S_Va_pre_NDA(1,:)=zeros(1,N_new);
        % S_Va_pre_NDA(2,:)=zeros(1,N_new);
        % S_Va_pre_FDA(1,:)=S_Va_pre(1,:);
        % S_Va_pre_FDA(2,:)=S_Va_pre(2,:);
        S_Va_pre_NDA((1:PT),:)=zeros(PT,N_new);  % ����ͬ������
        S_Va_pre_FDA((1:PT),:)=S_Va_pre((1:PT),:);
        for i=1:NG
            % ��������ڵ�λ��ָ��ͬ�����򣬼� DA_GN_flag(i) = 1
            % ��״̬������NDA���еĶ�Ӧ��ճϵ�������㣬ʣ��״̬������FDA����ԭʼ״̬������Ӧ�и���
            if (DA_GN_flag(i)==1)
                S_Va_pre_NDA((i+PT),:)=zeros(1,N_new);
                S_Va_pre_FDA((i+PT),:)=S_Va_pre((i+PT),:);
            else
                S_Va_pre_NDA((i+PT),:)=S_Va_pre((i+PT),:);
                S_Va_pre_FDA((i+PT),:)=zeros(1,N_new);
            end
        end

        % ��ȡ����ͬ����״̬��������
        S_Va_pre_FDA_sum_flag=sum(S_Va_pre_FDA,2); % ��� S_Va_pre_FDA λ�����飺�����мӺ�Ϊ�㣬��Ϊ��ͬ�������У�����Ϊͬ��������
        S_Va_pre_DA=zeros(M_new,N_new);
        DA_pre_cnt=0; % ��ʼ��������
        % ����ԭʼ״̬��������
        for i=1:M
            if (S_Va_pre_FDA_sum_flag(i)~=0)
                DA_pre_cnt=DA_pre_cnt+1;
                S_Va_pre_DA(DA_pre_cnt,:)=S_Va_pre_FDA(i,:); % ��ȡ����ͬ����״̬�������飬Ԥ������ = M_new
            end
        end
        
    else
        
        % ����������������ͬ��������ȡ�õ��Ľ���ͬ����״̬����������Ϊԭʼ����
        S_Va_pre_DA=S_Va_pre;
        
    end
    
    % ------------- ������������ͬ�� [ǰ����]״̬���� END --------------- %
    
    % -------------- �޸Ĳ��� + ���ϳ�Ա���㣨Ԥ�ⲽ�� END -------------- %
    
    if (flag_da_method==1)
        
        % -------------- �����㷨�٣�EnKF �˲��������¼��ϳ�Ա START ---------------- %
        
        % EnKF�㷨ʵ��
        % ע�⣺EnKF���й����У���Ҫ����ľ�������ά��Ϊ[M * M]���ʵ�״̬����ά�ȹ���ʱ��EnKF�޷�ʹ�ã�
        % ����һ��Ԥ�����
        % �����ʱ��Ԥ��ļ������Э���� C_pp_e_f -- dim[M * M]
        C_pp_e_f=cov(S_Va_pre_DA')+C_qq_e;
        
        % ����������²���������۲���Ϣ
        % �� ���㿨��������Ke
        Ke=C_pp_e_f*H'*pinv(H*C_pp_e_f*H'+C_ss_e);
        
        % �� ����۲���Ϣ�����¼���Ԥ��ֵΪ����ֵ
        % �Լ��ϳ�Ա���и���
        S_Va_ana_DA=S_Va_pre_DA+Ke*(D_obs_e-H*S_Va_pre_DA);
        
        % �� �������Э����
        C_pp_e_a=(ones(M_new,M_new)-Ke*H)*C_pp_e_f;
        % C_pp_e_a=cov(S_Va_ana_DA');
        
        % --------------- �����㷨�٣�EnKF �˲��������¼��ϳ�Ա END ----------------- %
            
    elseif (flag_da_method==2)
        
        % -------------- �����㷨�ڣ�ETKF �˲��������¼��ϳ�Ա START ---------------- %

        % ͬ���������ṩ��Ӧλ�õĹ۲���Ϣ����ʵģ�����Ҵ���
        % ETKF�㷨ʵ��
        % �����ʱ��Ԥ��ļ������Э����
        % C_pp_e_f=cov(S_Va_pre')+C_qq_e;

        % ����һ������ֵ�ֽ�
        S_Va_pre_DA_m=(S_Va_pre_DA*ones(N_new,N_new))/N_new;
        S_Va_pre_DA_det=S_Va_pre_DA-S_Va_pre_DA_m;
        R=eye(N_new)+(S_Va_pre_DA_det')*(H')*pinv((N_new-1).*C_ss_e)*H*S_Va_pre_DA_det; % ���� Moore-Penrose α�棬��֤���ȿ������ͬ�����ȶ���
        [ZL,G,ZR]=svd(R); % ����ֵ�ֽ�

        % ����������²���������ʵ��۲���Ϣ
        % �� ���㿨��������Ke��ά��ΪM*k
        Ke=S_Va_pre_DA_det*ZL*inv(G)*ZR'*(S_Va_pre_DA_det')*(H')*pinv((N_new-1).*C_ss_e);
        % �� ����任����T
        T=ZL*sqrtm(inv(G))*ZR';
        % �� ����۲���Ϣ�����¼���Ԥ��ֵΪ����ֵ
        % ��״̬���ϳ�Ա���и���
        if (flag_da_relative==1)
            S_Va_ana_DA_m=S_Va_pre_DA_m+Ke*(D_obs_e-H*S_Va_pre_DA_m);
        elseif (flag_da_relative==0)
            S_Va_ana_DA_m=S_Va_pre_DA_m+Ke*(D_obs_e_z-H*S_Va_pre_DA_m);
        end
        S_Va_ana_DA_det=S_Va_pre_DA_det*T;
        S_Va_ana_DA=S_Va_ana_DA_m+S_Va_ana_DA_det;
        % �������Э����
        % C_pp_e_a=(eye(M)-Ke*H)*C_pp_e_f;

        % --------------- �����㷨�ڣ�ETKF �˲��������¼��ϳ�Ա END ----------------- %
        
    end

    % ------------- ������������ͬ�� [��ԭ]״̬���� START --------------- %
    
    % ��������������ͬ�������뻹ԭ״̬����
    if (flag_da_part==1)
        
        % ��ԭͬ�����ȫ����״̬�������飬��ά��Ϊԭʼ״̬���������Ӧά��
        S_Va_ana_FDA=zeros(M,N_new);
        DA_ana_cnt=0; % ��ʼ��������
        for i=1:M
            if (S_Va_pre_FDA_sum_flag(i)~=0)
                DA_ana_cnt=DA_ana_cnt+1;
                S_Va_ana_FDA(i,:)=S_Va_ana_DA(DA_ana_cnt,:); % ��ȡ����ͬ����״̬��������
            end
        end
        
        % !!! ͬ��״̬�������� + ԭʼ��ͬ��״̬�������� = ���ӵõ�����ȫ����״̬������������ά��Ϊ M * N_new��!!!
        S_Va_ana=S_Va_ana_FDA+S_Va_pre_NDA;
        
    else
        
        % ����������������ͬ������ͬ����õ���״̬�������鼴Ϊȫ����״̬����
        S_Va_ana=S_Va_ana_DA;
        
    end
    
    % -------------- ������������ͬ�� [��ԭ]״̬���� END ---------------- %
    
    % -------------- ������ + ������̬��ճϵ�� START ---------------- %
    
    for kk=1:S_Va_target_Num
        
        % ע�⣺����ͬ����ճʱ����״̬�����д�����ճ���󣬶Բ�̬���ݽ������������������������
        if (strcmp(char(S_Va_target_Name(kk)),'TMu')==1)
            
            % ��ͬ���˲���analytical����������ճϵ���еĲ�̬���ݽ������������ָ�ԭʼϵ��ֵ
            for i=(PT+1):MDA
                for j=1:N_new

                    % ������̬����
                    if (S_Va_ana(i,j)<1E-10)
                        S_Va_ana(i,j)=1E-10;
                    end

                    % !!!! ȡ����������δ����
                    % if (S_Va_ana(i,j)<-10)
                    %     S_Va_ana(i,j)=-10;
                    % end
                    % S_Va_ana(i,j)=10^(S_Va_ana(i,j));  % ת�����в��־�����ʧ

                end
            end
            
        % elseif (strcmp(char(S_Va_target_Name(kk)),'Beta')==1)
             % ������ Beta�� ���ӻ���Visualization��
             % �ڵ���������ȡ S_Va_ana ÿ�ֽ���е���Ч���ݲ����п��ӻ�
        
        end
        
    end
    
    % --------------- ������ + ������̬��ճϵ�� END ----------------- %
    
    % ����������¼��ǰ���ŷ������� S_Va_ana_ave �� S_Va_ana_e
    S_Va_ana_ave=mean(S_Va_ana'); % ����õ����Ź���ֵ
    S_Va_ana_ave_da=S_Va_ana_ave((PT+1):MDA)'; % ��ȡͬ������۲ⲿ�ֵ�״̬��������TMu����Beta���ȣ�!!! 0418_�޸�BUG ����MDA !!!
    
    % ÿ��ͬ����� ������Beta�� ���ӻ�
    for kk=1:S_Va_target_Num
        if (strcmp(char(S_Va_target_Name(kk)),'Beta')==1)
            % ����MATLAB�Դ��� scatteredInterpolant() ������ͬ���õ�������Ԫֵ���������ֵ������ڵ㣬���� V = F(x,y) ���������
            V=scatteredInterpolant(GN_grid_x,GN_grid_y,S_Va_ana_ave_da);
            V.Method='natural';
            V.ExtrapolationMethod='nearest';
            Va_NNP_Beta_new=V(VN_grid_x,VN_grid_y); % ��ֵ�õ���������ڵ� Beta ֵ����������Ԫֵ��ֵ
        end
    end
    
    % ���浱ǰ������Ϣ
    target_renew_dat=['S_Va_round_',num2str(p+1,'%03d'),'.dat'];
    target_renew_mat=['S_Va_round_',num2str(p+1,'%03d'),'.mat'];
    target_renew_da_dat=['S_Va_round_da_',num2str(p+1,'%03d'),'.dat'];
    target_renew_da_mat=['S_Va_round_da_',num2str(p+1,'%03d'),'.mat'];
    save([new_output_folder,target_renew_dat],'S_Va_ana','-ascii'); % ������³�Ա�����������
    save([new_output_folder,target_renew_mat],'S_Va_ana'); % ������³�Ա�����������
    save([new_output_folder,target_renew_da_dat],'S_Va_ana_DA','-ascii'); % �������ͬ����Ա�����������
    save([new_output_folder,target_renew_da_mat],'S_Va_ana_DA'); % �������ͬ����Ա�����������
    
    % ����ͬ���˳���־
    % ���� L2Norm���� �� ��ѭ�������ﵽ�����������˳���ѭ������
    % ����дING������
    
    % �ڵ�����Ŀ����
    p=p+1;
    
    % ----------------------- �ؼ������ START ------------------------ %

    S_Va_ana_final=S_Va_ana_ave';
    DA_Results_targetPath=[targetPath,'DA_Final_Results'];
    State_Ensem_Recal_New_Common_Beta_BL(p,S_Va_ana_final,MStepNS_recal_ensemble,DA_Results_targetPath);
    save([savename,'_DA_RUN_All.mat']);  % ���浱ǰͬ���������н��

    % ----------------------- �ؼ������ END ------------------------ %
    
end

% ----------------------- �ڵ������� END ------------------------ %

%% ������������ͬ���������
save([savename,'_All.mat']);  % ����ͬ���������н��
fprintf('\n');
fprintf(fid_info,'\n');
fprintf('����ͬ���������ʼ\n');
fprintf(fid_info,'����ͬ���������ʼ\n');
fprintf('\n');
fprintf(fid_info,'\n');

% �����ڵ���������������������
new_output_folder_cov=[new_output_folder,'Results_Recal_Covergence_Analysis'];
mkdir(new_output_folder_cov);
targetPath=new_output_folder_cov;

fid_cal_force_all=fopen([new_output_folder_cov,'\Recal_Covergence_Analysis_Ensemble_Force_ALL.PLT'],'w+');

fprintf('�ڵ��������������������߽�����ڵ���\n');
fprintf(fid_info,'�ڵ��������������������߽�����ڵ���\n');

for i=1:N_new
    Cal_force_init=load([savename,'\Pre_Ensemble_Generation\Init_Ensemble_Member_',num2str(i,'%03d'),'\Force.PLT']);
    Cal_force_inner_iter=Cal_force_init;
    for p=1:Max_inner_iter
        Cal_force_new=load([savename,'\Inner_Iteration_Round_',num2str(p,'%03d'),'\Renew_Ensemble_Member_',num2str(i,'%03d'),'\Force.PLT']);
        Cal_force_temp=[Cal_force_inner_iter;Cal_force_new];
        Cal_force_inner_iter=Cal_force_temp;
    end
    % �����ڵ���������������������
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
    
    fprintf('���ϳ�Ա %03d �������������߽���������\n',i);
    fprintf(fid_info,'���ϳ�Ա %03d �������������߽���������\n',i);
    pause(2);

end
fprintf('�ڵ��������������������߽����ȫ������\n');
fprintf(fid_info,'�ڵ��������������������߽����ȫ������\n');
fprintf('\n');
fprintf(fid_info,'\n');

fclose(fid_cal_force_all);

% ��ճϵ������ϵ���Ա�
% ���巽��ϵ������ͬ�������ڸ�����ڵ㣩
% ��ȡ����ͬ�������Ա��ճϵ������
fprintf('��ճϵ������ϵ��������ڵ���\n');
fprintf(fid_info,'��ճϵ������ϵ��������ڵ���\n');
% ����ͬ��ǰ��ճ����
TMuSP=zeros(NNP,N_new);
TMu_coeff=zeros(NNP,1); % ������ճϵ������ϵ��
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
    % ������ճϵ������ϵ��
    TMu_coeff(m)=(sqrt(TMu_err_sum/(N_new-1)))/TMu_mean(m);
end
% �滻ԭʼ����ڵ���ճϵ�������ļ�
save([ddir_post,'\Input\TMuPData.PLT'],'TMu_coeff','-ascii');
% ���� CFD_2D_Post_TMu ���򣬲��������������ͼ
cmd_post=[ddir_post,'\CFD_2D_Post_TMu.exe'];
open(cmd_post);
pause(5);
pre_post_TMu_Path=[ddir_post,'\Output\TMuPData_Post.PLT'];
copyfile(pre_post_TMu_Path,[new_output_folder,'Pre_Ensemble_Generation\TMuPData_Post.PLT']);
fprintf('ͬ��ǰ״̬ ��ճϵ������ϵ������������\n');
fprintf(fid_info,'ͬ��ǰ״̬ ��ճϵ������ϵ������������\n');
pause(2);
% ������ڵ���������ճ����
for p=1:Max_inner_iter
    TMuSP=zeros(NNP,N_new);
    TMu_coeff=zeros(NNP,1); % ������ճϵ������ϵ��
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
        % ������ճϵ������ϵ��
        TMu_coeff(m)=(sqrt(TMu_err_sum/(N_new-1)))/TMu_mean(m);
    end
    % �滻ԭʼ����ڵ���ճϵ�������ļ�
    save([ddir_post,'\Input\TMuPData.PLT'],'TMu_coeff','-ascii');
    % ���� CFD_2D_Post_TMu ���򣬲��������������ͼ
    cmd_post=[ddir_post,'\CFD_2D_Post_TMu.exe'];
    open(cmd_post);
    % system([ddir_post,'\CFD_2D_Post_TMu.exe']);
    pause(5);
    pre_post_TMu_Path=[ddir_post,'\Output\TMuPData_Post.PLT'];
    copyfile(pre_post_TMu_Path,[new_output_folder,'Inner_Iteration_Round_',num2str(p,'%03d'),'\TMuPData_Post.PLT']);
    fprintf('�ڵ���״̬ %03d ��ճϵ������ϵ������������\n',p);
    fprintf(fid_info,'�ڵ���״̬ %03d ��ճϵ������ϵ������������\n',p);
    pause(2);
end
fprintf('��ճϵ������ϵ�������ȫ������\n');
fprintf(fid_info,'��ճϵ������ϵ�������ȫ������\n');
fprintf('\n');
fprintf(fid_info,'\n');

fclose(fid_info);
fclose(fid_all);
