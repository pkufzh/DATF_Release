
% �������ڣ�2021/05/27
% �������ܣ�ʵ������ͬ�����ϳ�Ա�ĳ�ʼ�������������У�
% �����������˵��
% ���������S_Va_init ���� ��ʼ�Ŷ�״̬���������������ӭ�ǣ����ų�����
%           N ���� ���ϳ�Ա��Ŀ
%           MStepNS_cov ����  ���ö���NS�����������������
%           targetPath ���� ����Ŀ��·��
% ���������S_Va_init_reg ���� �������ϳ�Ա��ʼ�����״̬�������󣬴��������ڵ�������
%           NNP ���� ͬ����������ڵ�ά��
%           NG ���� ͬ����������ά��
%           MDA ���� ͬ����������ά��
%           M ���� ���ϳ�Աά�ȣ�״̬������

% ��Ҫ���£��޶����������������ļ��ϳ�Ա����ͬ�����̣��޸���Ч���ϳ�Ա N_max �� N_new

function [S_Va_init_Combined_eff,target_width,Val_Sen,N_new,NNP,NG,MDA,M]=State_Ensem_Gen_Parellel_Common_Beta_BL(S_Va_init,N_max,MStepNS_cov,targetPath)

% ����ȫ�ֱ���
global savename new_output_folder fid_all fid_info PPNS initPath CFL_Setting Resi_target Resi_target_lev JX_NS STD_Comp JX_Turb L2NORM_Setting flag_Va_init;
global S_Va_turb_Name S_Va_turb_Type S_Va_turb_Num S_Va_turb_FilePath S_Va_turb_ModiLine S_Va_turb_ModiFormat S_Va_turb_Std;
global S_Va_target_Name S_Va_target_Num S_Va_target_TypeTurb S_Va_target_TypeDA S_Va_turb_target S_Va_Flow_length PT PS;
global S_Va_ref_Name S_Va_ref_Num;
global shape_upper shape_lower shape_outline shape_combined_arr_init Va_Beta Va_NNP_Beta Va_NNP_Beta_new;
global flag_da_relative flag_sigma_obs_adapt flag_sigma_obs_coeff D_obs_e_ori D_obs_e flag_flow_type;
global BL_setting N_BL L_BL BL_num BL_lab BL_dot BL_obs_num kb;
global Resi_flag Resi_flag_label;

% --------------- �޸Ĳ���Ԥ��+���ϳ�Ա���� START --------------- %

% ��ʼ��target_width
target_width=-1;

% ��ʼ��������ָ��
TMu_RMSE_sum=0;

% ���嶨�����ϳ�Ա������������
% Resi_target_gen_level=Resi_target_log(1)+1; % �涨�������ɽضϲв�ˮƽΪ�����в�ֵ + 1���� 10 Ϊ�ף��������󼯺ϳ�Ա��׽��Χ
Resi_flag=zeros(1,N_max);  % Resi_flag_resi = 0 �� 1����ʾ ��Ӧ��ŵ��Ŷ�״̬�����Ƿ��������� - 1���� - 2����
Resi_label_pre=0;
Resi_flag_label=zeros(1,N_max);  % Resi_flag_label = k����ʾ��Ӧ�����ŵ��Ŷ�״̬�����ԭʼ��Ա��� k��

% �޸�NS����������������MStepNS -- NS��������������
% ָ���滻�У���Ӧ�����ļ���InputNS.dat����9��
rep_line_MStepNS=9;
rep_format_MStepNS='%-16d';

% �޸����㿪�أ� JX_NS -- �Ƿ���������������㣨0��1��
%                JX_Turb -- �Ƿ���������������㣨0��1��
% ָ���滻�У��ֱ��Ӧ�����ļ���InputNS.dat����10�У���InputTurb.dat����7��
rep_line_JX_NS=10;
rep_format_JX_NS='%-12d';
rep_line_JX_Turb=7;
rep_format_JX_Turb='%-12d';

% �޸ļ���CFL����CFL_Setting
% ָ���滻�У���Ӧ�����ļ���InputNS.dat����8��
rep_line_CFL=8;
rep_format_CFL='%-12d';
CFL_Gen=CFL_Setting(1);

% �޸ı�׼�����������ȽϿ��أ�STD_Comp
% ָ���滻�У���Ӧ�����ļ���InputNS.dat����22��
rep_line_STD_Comp=22;
rep_format_STD_Comp='%-12d';
STD_Comp_Gen=STD_Comp(1);

% �޸������в Resi_target
rep_line_Resi=15;
rep_format_Resi='%-10.2e';
Resi_Gen=Resi_target(1);

% �޸����� L2NORM ˮƽ�� L2NORM_Setting
rep_line_L2NORM=8;
rep_format_L2NORM='%-12d';
L2NORM_Gen=L2NORM_Setting(1);

addPath_pre=[targetPath,'Pre_Ensemble_Generation'];

% ѭ����������������ʵ���޸Ĳ���+���ϳ�Ա����
% ��ʼ�������k�����ϳ�Ա
k_flag=0;
% run_flag(i)   �����i���������ǰ�����״̬���
%           =0  ����ǰ��������ڿ���״̬
run_flag=zeros(1,PPNS);
cnt_flag=zeros(1,PPNS);
% stop_flag(i) = 0 �����i�����������
%              = 1 �����i����������ڹ���
stop_flag=zeros(1,PPNS);
% exit_flag = 0 �����������δ����
%           = 1 ������������ѽ���
exit_flag=0;

for kp=1:PPNS
    % ����٣���� target_run_Path Ŀ¼���Ƿ���� ProgramExitFlag.txt
    %         �����ڣ�������ɾ��
    % ����������н�����־�ļ�ģ��
    % ���ö�д·��
    deletePath=['Program_CFD_Parallel_',num2str(kp,'%03d')];
    if (exist([deletePath,'\ProgramExitFlag.txt'],'file'))~=0
        delete([deletePath,'\ProgramExitFlag.txt']);
    end
end

while (exit_flag~=1)
    
    % ���ñ�־�����Ҳ������е����������MATLAB������ͣɨ�����룬��������Ч��
    % ��ʼ״̬�������������
    % for i=1:PPNS
    %    stop_flag(i)=0;
    % end
    
    for kp=1:PPNS
    
        % ���õ�ǰ״̬��д·��
        target_runPath=['Program_CFD_Parallel_',num2str(kp,'%03d')];
        sourcePath=[target_runPath,'\Output'];

        % ����ڣ�ִ�в��д���ģ��
        if ((run_flag(kp)==0) && (k_flag<N_max+1))
            
            k_flag=k_flag+1;
            if (k_flag==N_max+1)
                break;
            end
            run_flag(kp)=k_flag;  % ��ǰ����״̬ʵ�ʼ����ڱ�� run_flag(kp)���� k_flag
            
            % �������޸�
            % 1. ���£��رգ�������ݿ���
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
            
            % 2. ����NS����������е�������
            fid2=fopen([target_runPath,'\Input\InputNS.dat'],'r+');
            for i=1:(rep_line_MStepNS-1)
               fgetl(fid2);
            end
            fseek(fid2,0,'cof');
            fprintf(fid2,rep_format_MStepNS,MStepNS_cov);
            fclose(fid2);
            
            % 3. ���ó�ʼ�Ŷ�����
            if (flag_Va_init==0)
                % �������������� Beta ͬ��
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
                % ������������ Beta ͬ������д����Ӧ��Ϣ��ȫ������Ԫ��ڵ� Beta ֵ��ע�⣺����ܶ�����Ԫ�� Beta ����ͬ����
                Beta_init=Va_Beta(:,run_flag(kp));
                save([sourcePath,'\Beta_plus.PLT'],'Beta_init','-ascii');
                Beta_init_NNP=Va_NNP_Beta(:,run_flag(kp));
                save([sourcePath,'\Beta_plus_NNP.PLT'],'Beta_init_NNP','-ascii');
            end
            
            % 4. ���ü���CFL��
            fid4=fopen([target_runPath,'\Input\InputNS.dat'],'r+');
            for i=1:(rep_line_CFL-1)
               fgetl(fid4);
            end
            fseek(fid4,0,'cof');
            fprintf(fid4,rep_format_CFL,CFL_Gen);
            fclose(fid4);
            
            % 5. �������׼�����ȽϿ���
            fid5=fopen([target_runPath,'\Input\InputNS.dat'],'r+');
            for i=1:(rep_line_STD_Comp-1)
               fgetl(fid5);
            end
            fseek(fid5,0,'cof');
            fprintf(fid5,rep_format_STD_Comp,STD_Comp_Gen);
            fclose(fid5);
            
            % 6. ���ü�����������
            fid6=fopen([target_runPath,'\Input\InputNS.dat'],'r+');
            for i=1:(rep_line_Resi-1)
               fgetl(fid6);
            end
            fseek(fid6,0,'cof');
            fprintf(fid6,rep_format_Resi,Resi_Gen);
            fclose(fid6);
            
            % 7. ����L2NORM����ˮƽ
            fid7=fopen([target_runPath,'\Input\InputTurb.dat'],'r+');
            for i=1:(rep_line_L2NORM-1)
               fgetl(fid7);
            end
            fseek(fid7,0,'cof');
            fprintf(fid7,rep_format_L2NORM,L2NORM_Gen);
            fclose(fid7);

            % 8. ���� CFD �������ɳ�ʼ���ϳ�Ա����������
            % ����CFD���������
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
                
            fprintf(['���ϳ�Ա���ɹ��̣�״̬ ',num2str(run_flag(kp),'%03d')]);
            fprintf(fid_info,['���ϳ�Ա���ɹ��̣�״̬ ',num2str(run_flag(kp),'%03d')]);
            if (flag_Va_init==0)
                % �������������� Beta ͬ��
                for ns=1:S_Va_turb_Num
                    fprintf([' ',char(S_Va_turb_Name(ns)),' = %.5f'],S_Va_init(run_flag(kp),ns));
                    fprintf(fid_info,[' ',char(S_Va_turb_Name(ns)),' = %.5f'],S_Va_init(run_flag(kp),ns));
                end
            else
                % ������������ Beta ͬ��
                fprintf(' ��ʼ�Ŷ� Beta ��');
                fprintf(fid_info,' ��ʼ�Ŷ� Beta ��');
            end
            fprintf(' ����������\n');
            fprintf(fid_info,' ����������\n');

            % ����CFD�����������д��ڣ�������С�������������г���
            % cmd=[target_runPath,'\CFD_2D.exe'];
            % open(cmd);
            [s,e]=dos(['cd ',target_runPath,' && run_CFD_2D.bat&']);
            stop_flag(kp)=1;
            pause(3); % ��ͣ����
            
        end
        
        % ����ۣ��ж��Ƿ�����Ѿ��������
        % ���ҵ� ProgramExitFlag.txt��������Ѿ�������ϣ����б���
        if (exist([target_runPath,'\ProgramExitFlag.txt'],'file'))~=0
            
            % �����ʼ���ϳ�Ա��Ϣ
            % �����Աָ����Ԫ������ճϵ��

            % -------------------- ���н������ START ------------------- %
            
            % ��ȡ��ǰ״̬�������ղв�ֵ Resi(dual)_Final
            if (JX_Turb(1)==0)
                Resi_All=load([target_runPath,'\Output\Residual_1.PLT']);
            elseif ((JX_Turb(1)==1) || (JX_Turb(1)==2))
                Resi_All=load([target_runPath,'\Output\Residual_1_Turb_Const.PLT']);
            elseif (JX_Turb(1)==3)
                Resi_All=load([target_runPath,'\Output\Residual_1_Turb_Beta.PLT']);
            end
            Resi_Final=Resi_All(size(Resi_All,1),2);
            % �����������
            com_data=load([sourcePath,'\CP_CF_Comp.PLT']);
            % ָ��������������ݱ���·��
            fid_com_adj_path=[sourcePath,'\CP_CF_Comp_New.PLT'];
            fid_tar_new_path=[sourcePath,'\CP_New.PLT'];
        
            % �Լ���ѹ��ϵ�����߽ڵ��Ž���������ʵʱ�����ǰͬ���ִ�ѹ���ֲ�����
            % �������ɹ��̡��ڵ���������ȫ����Ա��ɼ���󱣴����߰�����
            if (flag_flow_type==0)
                [shape_upper_num,shape_lower_num,shape_upper_arr_sor,shape_lower_arr_sor,n_com]=DA_Va_Airfoil_Ref_Adjustment_and_Visulization(com_data,fid_com_adj_path,fid_tar_new_path);
            elseif (flag_flow_type==1)
                [shape_num,shape_arr_sor,n_com]=DA_Va_NonAirfoil_Ref_Adjustment_and_Visulization(com_data,fid_com_adj_path,fid_tar_new_path);
            end
            
            % �����������ʾ
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
            
            fprintf(['���ϳ�Ա���ɹ��̣�״̬ ',num2str(run_flag(kp),'%03d'),' �ѱ�����ϣ�Total Step = %d��Final Residual = %.5f\n'],(size(Resi_All,1)-1),Resi_Final);
            fprintf(fid_info,['���ϳ�Ա���ɹ��̣�״̬ ',num2str(run_flag(kp),'%03d'),' �ѱ�����ϣ�Total Step = %d��Final Residual = %.5f\n'],(size(Resi_All,1)-1),Resi_Final);
            
            % ��Ǽ���ﵽ������������� Resi_flag = 1
            if (Resi_Final<=Resi_target_lev(1))
                Resi_flag(run_flag(kp))=1; % ����״̬��Ϊ 1��������״̬��Ϊ 0
                Resi_label_pre=Resi_label_pre+1; % ����״̬�� +1
                Resi_flag_label(Resi_label_pre)=run_flag(kp); % ���浱ǰ����״̬��� Resi_flag_label���ۼ��� = Resi_label_pre������ = N_new
            end
            
            if (strcmp(char(S_Va_ref_Name(1)),'BL')==1)
                
                % ������ȡ�߽����Ϣ�������ļ�
                fprintf(['���ϳ�Ա���ɹ��̣�״̬ ',num2str(run_flag(kp),'%03d'),' �߽��������ȡ��ʼ\n']);
                fprintf(fid_info,['���ϳ�Ա���ɹ��̣�״̬ ',num2str(run_flag(kp),'%03d'),' �߽��������ȡ��ʼ\n']);
                [s,e]=dos(['cd ',target_runPath,' && run_BL_Hump_Extract_Tecplot.bat']);
                fprintf(['���ϳ�Ա���ɹ��̣�״̬ ',num2str(run_flag(kp),'%03d'),' �߽��������ȡ���\n']);
                fprintf(fid_info,['���ϳ�Ա���ɹ��̣�״̬ ',num2str(run_flag(kp),'%03d'),' �߽��������ȡ���\n']);
                pause(3); % ��ͣ����
                
                % ��ȡ�߽���ٶ�������
                cnt_BL_cal=0;
                BL_Cal_Vel=zeros(kb,1);
                for i=1:N_BL
                    pre_BL=[sourcePath,'\BL_Output\BL_Hump_xc_',num2str(BL_lab(i),'%03d'),'.dat'];
                    BL_info=load(pre_BL);
                    BL_Vel=BL_info(:,3); % ��ȡ�����߽���ٶ���Ϣ
                    for j=1:BL_num(i)
                        cnt_BL_cal=cnt_BL_cal+1;
                        BL_Cal_Vel(cnt_BL_cal)=(BL_Vel(BL_dot(i,j)))/0.1; % ��ȡ����ͬ���ı߽���ٶ���Ϣ
                    end
                end
                if (sum(cnt_flag)==0)
                    BLSC=zeros(length(BL_Cal_Vel),N_max);
                    BLSC(:,run_flag(kp))=BL_Cal_Vel;
                else
                    BLSC(:,run_flag(kp))=BL_Cal_Vel;
                end
                
            end
             
            % ����ļ���棨��ֹ���ǣ�
            pre_Path=[addPath_pre,'\',initPath,'_',num2str(run_flag(kp),'%03d')];
            fileList=dir(sourcePath);
            mkdir(pre_Path); % ���� mkdir() �����������ļ���
            filename_s=cell(length(fileList),1);
            filename_t=cell(length(fileList),1);

            for i=3:length(fileList)
                  filename_s{i}=[sourcePath,'\',fileList(i).name];
                  filename_t{i}=pre_Path;
                  copyfile(filename_s{i},filename_t{i});
            end
            
%             if (strcmp(char(S_Va_ref_Name(1)),'BL')==1)
%                 % ��ȡ�߽���ٶ�������
%                 cnt_BL_cal=0;
%                 BL_Cal_Vel=zeros(kb,1);
%                 for i=1:N_BL
%                     pre_BL=[pre_Path,'\BL_Output\BL_Hump_xc_',num2str(BL_lab(i),'%03d'),'.dat'];
%                     BL_info=load(pre_BL);
%                     BL_Vel=BL_info(:,3); % ��ȡ�����߽���ٶ���Ϣ
%                     for j=1:BL_num(i)
%                         cnt_BL_cal=cnt_BL_cal+1;
%                         BL_Cal_Vel(cnt_BL_cal)=(BL_Vel(BL_dot(i,j)))/0.1; % ��ȡ����ͬ���ı߽���ٶ���Ϣ
%                     end
%                 end
%                 if (sum(cnt_flag)==0)
%                     BLSC=zeros(length(BL_Cal_Vel),N_max);
%                     BLSC(:,run_flag(kp))=BL_Cal_Vel;
%                 else
%                     BLSC(:,run_flag(kp))=BL_Cal_Vel;
%                 end
%             end
            
            % ��¼��ǰ�����ļ���·������������ļ�����
            pre_VcPath=[pre_Path,'\VcData.PLT'];
            pre_TMuPath=[pre_Path,'\TMuData.PLT'];
            pre_TMuPPath=[pre_Path,'\TMuPData.PLT'];
            pre_Cp_Cf_yplus_Path=[pre_Path,'\CP_CF_Comp_New.PLT']; % ע�⣺�Ѿ��޸�Ϊ����˳���ļ������ļ�
            pre_Beta=[pre_Path,'\Beta_plus.PLT'];
            pre_ErrInfo=[pre_Path,'\FlowField_Err_Info.txt'];

            % ��ȡ����ģ�����ݱ��� Vc �� TMu
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
            % �����������ż���������ָ��
            if (Resi_flag(run_flag(kp))==1)
                TMu_RMSE_sum=TMu_RMSE_sum+TMu_RMSE;
            end

            % ����˳����Ҫ�ر�ע�⣡�������޸�
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

            % -------------------- ���н������ END -------------------- %
            
            run_flag(kp)=0;  % �����������Ϊ����״̬
            delete([target_runPath,'\ProgramExitFlag.txt']);
            cnt_flag(kp)=cnt_flag(kp)+1; % �ۼƸ���������
            stop_flag(kp)=0;
            pause(3); % �����ļ�����ͣ����
            
        end
        
        if (sum(stop_flag)==PPNS)
            % �Ҳ������е��������MATLAB������ͣɨ������
            pause(10);
        end
        
        if ((sum(stop_flag)==0) && (sum(cnt_flag)==N_max))
            % ������������У��������˳�����������ɲ��г���������
            exit_flag=1;
            pause(10);
            break;
        end
    
    end 
    
end

% ����״̬���������������󣬼����µļ��ϳ�Ա���� N_new
N_new=sum(Resi_flag);
% ������������� Resi_flag_label ��β���޹���ɾȥ
Resi_flag_label((N_new+1):end)=[];
% ������������������ָ��
Val_Sen=TMu_RMSE_sum/N_new;

% !!!!!!!!!!!!!!!!!! ���޸�ģ�飬��ÿ����״̬�������ݼ���������������޸�
% % �����ʼ���º�״̬����
% S_Va_init_reg=[(S_Va_init(:,(1:PT))),TMuSC',CpSC']';
% NNP=length(TMuPS);   % ����ȫ��������ڵ�ά��
% NG=length(TMuS);   % ����ȫ��������ά��
% MDA=length(TMuS)+PT;   % ����ȫ������������ά��
% M=length(S_Va_init_reg);  % ���㼯�ϳ�Աά�ȣ�״̬������

% ͬ��״̬��������������ɣ��� ��ͬ���������� ����ͬ���������� �۲����
% ���·ֱ�����Ϊ S_Va_init_Part_One(Two,Three)
% ���յõ��ĺϲ����ͬ��״̬��������Ϊ S_Va_init_combined

% �� ��ȡ��ͬ������
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

% �� ��ȡ����ͬ������
S_Va_cnt_Two=0;
% ��˳�򱣴������ͬ������ά�ȣ����ȣ�
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

% �� ��ȡͬ���۲����
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

% �������������Ϣ
save([targetPath,'Com_Ensem_Gen.dat'],'S_Va_init_Part_Three','-ascii');

if ((strcmp(char(S_Va_ref_Name(1)),'Cp')==1) || (strcmp(char(S_Va_ref_Name(1)),'Cf')==1) || (strcmp(char(S_Va_ref_Name(1)),'yplus')==1))

    % !!!!! �޸Ĳ��֣���ȡ��Ч���ֵĹ۲���Ϣ !!!!!
    S_Va_init_Part_Three_eff=zeros(N_new,size(S_Va_init_Part_Three,2));
    n_combined_eff=0;
    for i=1:N_max
        if (Resi_flag(i)==1)
            n_combined_eff=n_combined_eff+1; % �����жϼ�������һ
            S_Va_init_Part_Three_eff(n_combined_eff,:)=S_Va_init_Part_Three(i,:);
        end
    end

    % �����м���״̬����õ���Ŀ��ͬ���������ֲ����İ����߽��л���
    S_Va_init_Part_Three_min=min(S_Va_init_Part_Three_eff); % �±߽�
    S_Va_init_Part_Three_max=max(S_Va_init_Part_Three_eff); % �ϱ߽�
    S_Va_init_Part_Three_mid=(S_Va_init_Part_Three_min+S_Va_init_Part_Three_max)/2; % ���°���������
    target_width=S_Va_init_Part_Three_max-S_Va_init_Part_Three_min; % ���°���߽���

    if (flag_flow_type==0)

        y1=S_Va_init_Part_Three_max(1:shape_lower_num);
        y2=S_Va_init_Part_Three_min(1:shape_lower_num);
        y3=S_Va_init_Part_Three_max(shape_lower_num+1:end);
        y4=S_Va_init_Part_Three_min(shape_lower_num+1:end);
        % �½����ӻ�ͼ��
        % !!! ���Ӱ����������� !!! Middle_line
        figure;
        hold on
        [hha1,hhb1,hhc1]=shadedplot(Xp(1:shape_lower_num),y1,y2,'[0.8706 0.9216 0.9804]');
        [hha2,hhb2,hhc2]=shadedplot(Xp(shape_lower_num+1:end),y3,y4,'[0.8706 0.9216 0.9804]');
        hh1=plot(shape_combined_arr_init(:,2),shape_combined_arr_init(:,7),'ro','MarkerSize',4.0); % ���Ʊ�׼ʵ��ֵ
        hh2=plot(shape_combined_arr_init(:,2),S_Va_init_Part_Three_mid,'-','LineWidth',0.8,'Color',[0.7451,0.7451,0.7451]); % ���ư�����������
        % plot(Xp,S_Va_init_Part_Three_max); % ���������������
        % plot(Xp,S_Va_init_Part_Three_min); % ���������������
        legend([hhb1,hhc1,hha1(2),hh2,hh1],{'���������������','���������������','���±߽���Χ����','���������������','ʵ��ֵ����ֵ'}); % ��Ӳ���ͼ��
        xlabel('\itx'),ylabel('\itCp'),title('�������ɼ�����');
        grid on
        box on
        hold off
        % ����ɱ༭ͼƬ
        saveas(gcf,['.\Output_Figures\',savename,'\Init_Ensem_Gen_Output.fig']);

    elseif (flag_flow_type==1)

        y1=S_Va_init_Part_Three_max;
        y2=S_Va_init_Part_Three_min;
        % �½����ӻ�ͼ��
        % !!! ���Ӱ����������� !!! Middle_line
        figure;
        hold on
        [hha1,hhb1,hhc1]=shadedplot(Xp,y1,y2,'[0.8706 0.9216 0.9804]');
        hh1=plot(shape_combined_arr_init(:,2),shape_combined_arr_init(:,7),'ro','MarkerSize',4.0); % ���Ʊ�׼ʵ��ֵ
        hh2=plot(shape_combined_arr_init(:,2),S_Va_init_Part_Three_mid,'-','LineWidth',0.8,'Color',[0.7451,0.7451,0.7451]); % ���ư�����������
        legend([hhb1,hhc1,hha1(2),hh2,hh1],{'���������������','���������������','���±߽���Χ����','���������������','ʵ��ֵ����ֵ'}); % ��Ӳ���ͼ��
        xlabel('\itx'),ylabel('\itCp'),title('�������ɼ�����');
        grid on
        box on
        hold off
        % ����ɱ༭ͼƬ
        saveas(gcf,['.\Output_Figures\',savename,'\Init_Ensem_Gen_Output.fig']);

    end
    
end

% �����ʼ���º�ȫ״̬�����������ֵ��ӣ�
if (flag_da_relative==1)
    S_Va_init_Combined=[S_Va_init_Part_One,S_Va_init_Part_Two,S_Va_init_Part_Three]';
elseif (flag_da_relative==0)
    % !!!! 2021/03/03��˼ģ����� !!!!
    % �Թ۲������ֵ����ת����Ŀ���ǵõ���ʵ��ֵ����ֵ����������ֵ
    % �ɽ�ʵ��ֵ����ֵ������Ϊ��
    % ͬ��Ŀ�꼴ת��Ϊ��С�۲�ֵ��ʵ��ֵ����ֵ���������������
    S_Va_init_Part_Three_err=(S_Va_init_Part_Three./D_obs_e_ori')-ones(size(D_obs_e_ori'));
    S_Va_init_Combined=[S_Va_init_Part_One,S_Va_init_Part_Two,S_Va_init_Part_Three_err]';
end

% -- ���ݲ��� -- %
save([new_output_folder,savename,'_S_Va_init_Combined.mat'],'S_Va_init_Combined'); % ��������ԭʼ��Ա�����������
% -- ���ݲ��� -- %

% % !!!!! �޸Ĳ��֣���ȡ��Ч���ֵ�״̬���� !!!!!
% S_Va_init_Combined_eff=zeros(N_new,size(S_Va_init_Combined,1));
% n_combined_eff=0;
% for i=1:N_max
%     if (Resi_flag(i)==1)
%         n_combined_eff=n_combined_eff+1;
%         S_Va_init_Combined_eff(n_combined_eff,:)=S_Va_init_Combined(i,:);
%     end
% end

% !!!!! �޸Ĳ��֣���ȡ��Ч���ֵ�״̬���� !!!!!
S_Va_init_Combined_eff=zeros(size(S_Va_init_Combined,1),N_new);
n_combined_eff=0;
for i=1:N_max
    if (Resi_flag(i)==1)
        n_combined_eff=n_combined_eff+1;
        S_Va_init_Combined_eff(:,n_combined_eff)=S_Va_init_Combined(:,i);
    end
end

% -- ���ݲ��� -- %
save([new_output_folder,savename,'_S_Va_init_Combined_eff.mat'],'S_Va_init_Combined_eff'); % ��������ԭʼ��Ա�����������
% -- ���ݲ��� -- %

% �����ʼ���º�ȫ״̬���� - ������ֵ
NNP=length(TMuPS);   % ����ȫ��������ڵ�ά��
NG=length(TMuS);   % ����ȫ��������ά��
% ͬ�����ά�ȼ��㣬�����ݳ�����
MDA=PT+size(S_Va_init_Part_Two,2);   % ����ȫ��ͬ����������ά�ȣ���ֵ���� + ������ֵ������
M=size(S_Va_init_Combined_eff,1);  % ���㼯�ϳ�Աά�ȣ�״̬������

% --------------- �޸Ĳ���Ԥ��+���ϳ�Ա���� END --------------- %

fprintf('\n');
fprintf(fid_info,'\n');

end
