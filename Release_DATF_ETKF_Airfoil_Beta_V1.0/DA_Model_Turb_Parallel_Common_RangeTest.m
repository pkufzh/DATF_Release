%% ����ͬ���������㸨������Data Assimilation & Turbulent Flows��
%% ���ϱ任�������˲�ʵ�飨Ensemble Transform Kalman Filter��ETKF����
%% ����ETKF�������ʵ�����ݣ�ʵ�ֶ���������������������ͬ��
%% ������Ҫ���ܣ��ܹ����������������ָ����Χ��ָ���������Ŷ����������׼����ֵ��������
%%              �����ڶ�ָ�������ĳ�ʼ�����Ŷ��������

%% ��ʼ������
clear all
clc
close all

% �����ʼ��Ա����
% load('Pre_Cal_Path.mat');

format long;

% ����ȫ�ֱ���
global savename fid_all fid_info PPNS initPath CFL_Setting;
global S_Va_turb_Name S_Va_turb_Num S_Va_turb_FilePath S_Va_turb_ModiLine S_Va_turb_ModiFormat S_Va_turb_Std S_Va_turb_Det_base_mark;
global S_Va_target_Name S_Va_target_Num S_Va_target_TypeTurb S_Va_target_TypeDA S_Va_turb_target S_Va_Flow_length PT PS;
global S_Va_ref_Name S_Va_ref_Num;
global shape_upper shape_lower shape_combined_arr_init;

%% ������Ҫ����ͬ�����̲���

% ------------  ������Ҫ�ļ������벢����ز��� START ------------ %

% ���ø����������ݽ���·��
savename='Results_Model_Turb_Parallel_NACA4412_RangeTest_Cb1_Sigma_Cv1_CKS';
new_output_folder=['.\',savename,'\'];
mkdir(new_output_folder);
targetPath=new_output_folder;
initPath='Init_Ensem_Member';

% �򿪼�¼������Ϣ�ļ�
fid_info=fopen([new_output_folder,'DA_CompInfo_Brief.dat'],'w+');
fid_all=fopen([new_output_folder,'DA_CompInfo_Detailed.dat'],'w+');

% ���ò����߳�/���������ɸ��ݹ������߳����ã�
PPNS=7;

% ��������ԭʼ�������ݣ��������������Ͳ�ֵ
shape_upper=load('.\Input_Data\Input_Data_Shape_NACA4412_upper.dat');
shape_lower=load('.\Input_Data\Input_Data_Shape_NACA4412_lower.dat');
% �������������ʵ������
exp_upper_data=load('.\Input_Data\Input_Data_Exp_NACA4412_upper.dat');
exp_lower_data=load('.\Input_Data\Input_Data_Exp_NACA4412_lower.dat');

% -------------  ������Ҫ�ļ������벢����ز��� END ------------- %

% -------------  �������ɳ�ʼ���ϳ�Ա��ز��� START ------------- %

% ���弯�ϳ�Ա����
N=100;

% ----------  ָ����ʼ�Ŷ�������Ŀ��ͬ��������ͬ���۲����ģ�� START ---------- %

% �� ָ����ʼ�Ŷ�����ģ�飨ĿǰĬ��Ϊ�������Ŷ���
% ������׼״̬����
% �����ʼ�Ŷ��������ơ��ļ�·�����޸��������޸ĸ�ʽ�����������ʾ��һһ��Ӧ��
% ��ʼ�Ŷ���������
S_Va_turb_Name=["Cb1","Sigma","Cv1","CKS"];
S_Va_turb_CharName=char(S_Va_turb_Name);
% ��ʼ�Ŷ���������
S_Va_turb_Num=length(S_Va_turb_Name);
% ��ʼ�Ŷ������ļ�·����ʡ�����г����Ŀ¼��
% S_Va_turb_FilePath_raw=char('\Input\InputNS.dat','\Input\InputNS.dat','\Input\InputTurb.dat');
% S_Va_turb_FilePath=strip(S_Va_turb_FilePath_raw,'right');
S_Va_turb_FilePath=["\Input\InputSAconst.dat","\Input\InputSAconst.dat","\Input\InputSAconst.dat","\Input\InputSAconst.dat"];
% ��ʼ�Ŷ������޸���������Ӧ�ļ��ڣ�
S_Va_turb_ModiLine=[3,5,6,7];
% ��ʼ�Ŷ������޸ĸ�ʽ����Ӧ�ļ��ڣ�
S_Va_turb_ModiFormat=["%-16.5f","%-16.5f","%-16.5f","%-16.5f"];
% ��ʼ�Ŷ������Ŷ�����
% �����Ŷ���ʽ��0 - ����ֵ�����Ŷ���1 - ���������Ŷ������磺����50�����Ŷ�����Ϊ��׼ֵ��[50%,150%]��
S_Va_turb_way=[1,1,1,1];
S_Va_turb_Std=[0.1355,(2/3),7.1,0.41];  % ��ʼ�Ŷ�������׼ֵ
S_Va_turb_Det_base_lower=[10,40,80,10];  % ��ʼ�Ŷ�����������
S_Va_turb_Det_base_upper=[80,80,80,150];  % ��ʼ�Ŷ�����������

% S_Gap=8; % �����������״̬�Ŷ��������ٷֱȣ�����ֵ 5 �����������Ŷ�����ı����� 5% �����仯��
% S_Va_turb_Det_base_range=S_Va_turb_Det_base_lower+S_Va_turb_Det_base_upper;
% N=round(max(S_Va_turb_Det_base_range)/S_Gap)+1;
% S_Va_turb_Det_base_ratio=S_Va_turb_Det_base_range./(N-1); % ����ÿһ������ʵ�ʵĵȲ����仯�ʣ��ٷֱȣ�
% S_Va_turb_Det_base_mark=zeros(N,S_Va_turb_Num);
% for i=1:S_Va_turb_Num
%     S_Va_turb_Det_base_mark(:,i)=linspace(((-1)*S_Va_turb_Det_base_lower(i)),S_Va_turb_Det_base_upper(i),N);
% end

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

% S_Va_init=zeros(N,length(S_Va_turb_min));
% % �Ȳ�����ʼ�Ŷ���������
% for i=1:length(S_Va_turb_min)
%     Va_min=S_Va_turb_min(i);
%     Va_max=S_Va_turb_max(i);
%     Va_turb=linspace(Va_min,Va_max,N)';
%     S_Va_init(:,i)=Va_turb;
% end
% S_Va_init_ave=mean(S_Va_init);

% �����������������г�ʼ�Ŷ���������
Va_min=S_Va_turb_min;
Va_max=S_Va_turb_max;
% ��ʼͬ����������
PN=length(Va_min);
S_Va_coef=lhsdesign(N,PN);
% ����N*P���Ŷ���������
S_Va_init=(ones(N,1)*Va_min)+S_Va_coef.*(ones(N,1)*(Va_max-Va_min));
S_Va_init_ave=mean(S_Va_init);

S_Va_turb_Det_base_mark=zeros(N,S_Va_turb_Num);
for i=1:N
    for j=1:S_Va_turb_Num
        S_Va_turb_Det_base_mark(i,j)=((S_Va_init(i,j)-S_Va_turb_Std(j))/S_Va_turb_Std(j))*100;
    end
end

% �� ָ��Ŀ��ͬ������ģ��
% ȫ��Ŀ��ͬ����������
S_Va_target_Name=["TMu"];
S_Va_target_CharName=char(S_Va_target_Name);
% ȫ��Ŀ��ͬ����������
S_Va_target_Num=length(S_Va_target_Name);
% �ж�Ŀ��ͬ���������͢٣�TypeTurb����k - �ǵ�k����ʼ�Ŷ�������k~=0����0 - �ǳ�ʼ�Ŷ�����
% ע����Ŀ��ͬ�������ǳ�ʼ�Ŷ�����������Ҫ�ڸ��²���ֵ��Ϊԭʼ STD ֵ
S_Va_target_TypeTurb=[0];
% �ж�Ŀ��ͬ���������͢ڣ�TypeDA����1 - ����ֵ��ͬ��������2 - ��������ֵ��ͬ������
S_Va_target_TypeDA=[2];
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
S_Va_ref_CharName=char(S_Va_ref_Name);
% ȫ��ͬ���۲��������
S_Va_ref_Num=length(S_Va_ref_Name);

% -----------  ָ����ʼ�Ŷ�������Ŀ��ͬ��������ͬ���۲����ģ�� END ----------- %

% --------------  �������ɳ�ʼ���ϳ�Ա��ز��� END -------------- %

% -----------------  ������Ҫ���������� START ----------------- %

MStepNS_cov=50000;   % ����NS��������������������������ɹ��̣�
CFL_Setting=5; % �ֱ����ü��ϳ�Ա���ɡ��ڵ����Լ��ؼ�����̵� CFL ��
    
% ------------------  ������Ҫ���������� END ------------------ %

% ---------------  ���������Ҫ������Ϣ START --------------- %

% �����������ʾ
fprintf('******************************************************************\n');
fprintf('                      �����Ŷ�������Ҫ������Ϣ                    \n');
fprintf('******************************************************************\n');
fprintf('[���򱣴�·�� & ��������]\n');
fprintf(['�������·����',savename,'\n']);
fprintf(['������Ϣ�ļ���DA_CompInfo_Brief.dat\n']);
fprintf('�������߳�/������ PPNS = %d\n',PPNS);
fprintf('******************************************************************\n');
fprintf('[ͬ����Ҫ����]\n');
fprintf('���ϳ�Ա���� N = %d\n',N);
fprintf('******************************************************************\n');
fprintf('[��ʼ�Ŷ�������Ϣ]\n');
for i=1:S_Va_turb_Num
    fprintf(['��ʼ�Ŷ����� ',num2str(i,'%03d'),'��',S_Va_turb_CharName(:,:,i),' �� [ %.5f , %.5f ], STD = %.5f\n'],S_Va_turb_min(i),S_Va_turb_max(i),S_Va_turb_Std(i));
end
fprintf('******************************************************************\n');
fprintf('[ͬ�����̿���]\n');
fprintf('���ü������ɹ���������������%d, CFL = %d\n',MStepNS_cov,CFL_Setting(1));
fprintf('******************************************************************\n');
fprintf('\n');

% �ļ������ʾ
fprintf(fid_info,'******************************************************************\n');
fprintf(fid_info,'                      �����Ŷ�������Ҫ������Ϣ                    \n');
fprintf(fid_info,'******************************************************************\n');
fprintf(fid_info,'[���򱣴�·�� & ��������]\n');
fprintf(fid_info,['�������·����',savename,'\n']);
fprintf(fid_info,['������Ϣ�ļ���DA_CompInfo_Brief.dat\n']);
fprintf(fid_info,'�������߳�/������ PPNS = %d\n',PPNS);
fprintf(fid_info,'******************************************************************\n');
fprintf(fid_info,'[ͬ����Ҫ����]\n');
fprintf(fid_info,'���ϳ�Ա���� N = %d\n',N);
fprintf('******************************************************************\n');
fprintf(fid_info,'[��ʼ�Ŷ�������Ϣ]\n');
for i=1:S_Va_turb_Num
    fprintf(fid_info,['��ʼ�Ŷ����� ',num2str(i,'%03d'),'��',S_Va_turb_CharName(:,:,i),' �� [ %.5f , %.5f ], STD = %.5f\n'],S_Va_turb_min(i),S_Va_turb_max(i),S_Va_turb_Std(i));
end
fprintf(fid_info,'******************************************************************\n');
fprintf(fid_info,'[ͬ�����̿���]\n');
fprintf(fid_info,'���ü������ɹ���������������%d, CFL = %d\n',MStepNS_cov,CFL_Setting(1));
fprintf(fid_info,'******************************************************************\n');
fprintf(fid_info,'\n');

% ----------------  ���������Ҫ������Ϣ END ---------------- %

%% ����һ������ο�ʵ�������ּ������ݣ������в�ֵ����
% �������ݲ���
% �����������
com_data=load('.\Program_CFD_Base\Output\CP_CF_Comp.PLT');
% ָ��������������ݱ���·��
fid_com_adj_path='.\Program_CFD_Base\Output\CP_CF_Comp_New.PLT';
fid_tar_new_path='.\Program_CFD_Base\Output\CP_New.PLT';

% �Լ���ѹ��ϵ�����߽ڵ��Ž���������ʵʱ�����ǰͬ���ִ�ѹ���ֲ�����
% �������ɹ��̡��ڵ���������ȫ����Ա��ɼ���󱣴����߰�����
[shape_upper_num,shape_lower_num,shape_upper_arr_sor,shape_lower_arr_sor,n_com]=DA_Va_Ref_Adjustment_and_Visulization(com_data,fid_com_adj_path,fid_tar_new_path);

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

% % ������ʵ��ֵ��������ӳ���ϵ
% % ��¼����ͬ����ʵ�����ݵ���Ϣ���ɵ����ļ�
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
% % �������ͬ���Ĺ۲�ʵ��ֵ�ֲ�
% fid_exp_da=fopen('.\Program_CFD_Base\Output\InputExpData_for_DA.PLT','w+');
% fprintf(fid_exp_da,'VARIABLES="X","CP"\n');
% fprintf(fid_exp_da,'ZONE T="ZONE EXP for DA"\n');
% for i=1:k
%     fprintf(fid_exp_da,'%20.8f %20.8f\n',exp_for_DA(i,1),exp_for_DA(i,2));
% end
% fclose(fid_exp_da);

% ���Ʊ�׼����ģ�ͼ�����������ģ�廯��ʾ
copyfile('.\Program_CFD_Base\Output\CP_New.PLT','.\Program_CFD_Base\Output\CP_SA.PLT');

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
saveas(gcf,['.\Output_Figures\',savename,'_RangeTest_Init_State_Output.fig']);

%% ����һ�����г����鴴��

% ��ȡ�����������
% ��������ת��������
fprintf('������������ͬ���������򣨼��ϳ�Ա��ƣ�������\n');
fprintf(fid_info,'������������ͬ���������򣨼��ϳ�Ա��ƣ�������\n');
fprintf('\n');
fprintf(fid_info,'\n');
cmd_grid=('.\Program_CFD_Base\GridTreat_.exe');
open(cmd_grid);
pause(3);
fprintf('������������ת����ɣ�\n');
fprintf(fid_info,'������������ת����ɣ�\n');

% �������������ļ��У�����ΪPPNS��
for p=1:PPNS
    ddir=['Program_CFD_Parallel_',num2str(p,'%03d')];
    copyfile('Program_CFD_Base',ddir);
    % �״���������ת��
    cmd_grid=([ddir,'\GridTreat_.exe']);
    open(cmd_grid);
    pause(3);
    fprintf('���г����� %03d ����ת����ɣ�\n',p);
    fprintf(fid_info,'���г����� %03d ����ת����ɣ�\n',p);
end

fprintf('\n');
fprintf(fid_info,'\n');

%% ����������ʼ�����ϳ�Ա

% -------------- �޸Ĳ���Ԥ�� + ���ϳ�Ա���� START -------------- %

S_Va_init_Combined=State_Ensem_Gen_Parellel_Common_Assist_RangeTest(S_Va_init,N,MStepNS_cov,targetPath);

% ����ȫ״̬�������󣬳�ʼ���ϳ�Ա����
S_Va_ana=S_Va_init_Combined;

% --------------- �޸Ĳ���Ԥ�� + ���ϳ�Ա���� END --------------- %

% �����ʼ������Ϣ
% ����ɹ�����ֱ�ӵ��û�����Ĺ���������
save([new_output_folder,savename,'_Pre.mat']);  % ���漯�ϳ�Ա���ɹ������н��
save([new_output_folder,savename,'_S_Va_init.mat'],'S_Va_ana'); % ������³�Ա�����������
save([new_output_folder,'S_Va_init.dat'],'S_Va_ana','-ascii');
save([new_output_folder,'S_Va_init.mat'],'S_Va_ana');

fclose(fid_info);
fclose(fid_all);
