
% �������ڣ�2021/05/27
% �������ܣ�ʵ������ͬ���ڵ������̵�Ԥ�ⲽ�����������У�
% �����������˵��
% ���������p ���� ��ǰ�ڵ����ִ�
%           S_Va_ana ���� ��ǰԤ��õ���״̬��������
%           N ���� ���ϳ�Ա��Ŀ
%           MDA ���� ͬ����������ά��
%           MStepNS_pre ����  �����ڵ���Ԥ�ⲽ��
%           targetPath ���� �ڵ����������Ŀ��·��
% ���������S_Va_pre ���� �����ڵ���Ԥ�ⲽ������״̬�������󣬴��������˲����»���

function [S_Va_pre_eff,target_width,Val_Sen,N_new]=State_Ensem_Renew_Parellel_Common_Beta_BL(p,S_Va_ana,N_pre,MDA,MStepNS_pre,targetPath)

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

% ��ʼ�� target_width
target_width=-1;

% ��ʼ��������ָ��
TMu_RMSE_sum=0;

% ���嶨���ڵ�����Ա������������
Resi_renew_flag=zeros(1,N_pre);  % Resi_flag_resi = 0 �� 1����ʾ ��Ӧ��ŵ��Ŷ�״̬�����Ƿ��������� - 1���� - 2����
Resi_renew_label_pre=0;
Resi_renew_flag_label=zeros(1,N_pre);  % Resi_flag_label = k����ʾ��Ӧ�����ŵ��Ŷ�״̬�����ԭʼ��Ա��� k��

% �޸Ĳ���Ԥ�ò�����Ԥ�ⲽ�����¹۲����Cp��
% CpSC=[];

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
CFL_Renew=CFL_Setting(2);

% �޸ı�׼�����������ȽϿ��أ�STD_Comp
% ָ���滻�У���Ӧ�����ļ���InputNS.dat����22��
rep_line_STD_Comp=22;
rep_format_STD_Comp='%-12d';
STD_Comp_Renew=STD_Comp(2);

% �޸������в Resi_target
rep_line_Resi=15;
rep_format_Resi='%-10.2e';
Resi_Renew=Resi_target(2);

% �޸����� L2NORM ˮƽ�� L2NORM_Setting
rep_line_L2NORM=8;
rep_format_L2NORM='%-12d';
L2NORM_Renew=L2NORM_Setting(2);

addPath_pre=[targetPath,'Inner_Iteration_Round_',num2str(p,'%03d')];

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
    % stop_flag=0;
    
    for kp=1:PPNS
        
        % ���ö�д·��
        target_runPath=['Program_CFD_Parallel_',num2str(kp,'%03d')];
        sourcePath=[target_runPath,'\Output'];
        
        % ����ڣ�ִ�в��д���ģ��
        if ((run_flag(kp)==0) && (k_flag<N_pre+1))
                        
            k_flag=k_flag+1;
            if (k_flag==N_pre+1)
                break;
            end
            run_flag(kp)=k_flag;
            
            % �������޸�
            % 1. ���£��򿪣�������ݿ���
            % !!! ���㿪���������⣿
            % JX_NS=1;
            % JX_Turb=1;
            fid1=fopen([target_runPath,'\Input\InputNS.dat'],'r+');
            for i=1:(rep_line_JX_NS-1)
               fgetl(fid1);
            end
            fseek(fid1,0,'cof');
            fprintf(fid1,rep_format_JX_NS,JX_NS(2));
            fclose(fid1);
            fid1=fopen([target_runPath,'\Input\InputTurb.dat'],'r+');
            for i=1:(rep_line_JX_Turb-1)
               fgetl(fid1);
            end
            fseek(fid1,0,'cof');
            fprintf(fid1,rep_format_JX_Turb,JX_Turb(2));
            fclose(fid1);

            % 2.����NS����������е�������
            fid2=fopen([target_runPath,'\Input\InputNS.dat'],'r+');
            for i=1:(rep_line_MStepNS-1)
               fgetl(fid2);
            end
            fseek(fid2,0,'cof');
            fprintf(fid2,rep_format_MStepNS,MStepNS_pre);
            fclose(fid2);

            % 3. ���ø���ͬ������
            % ����Ŀ��ͬ������(S_Va_ana)����Ŀ��ͬ�������ǳ�ʼ�Ŷ�����������Ҫ�ڸ��²���ֵ��Ϊԭʼ STD ֵ
            S_Va_cnt_One=0;
            for ns=1:S_Va_target_Num
                % ��ȡĿ��ͬ���������ͣ�TypeTurb����k - �ǵ�k����ʼ�Ŷ�������k~=0����0 - �ǳ�ʼ�Ŷ�����
                % ��Ŀ��ͬ������Ϊ��ʼ�Ŷ���������Ϊ��ͬ��������e.g. "Ma" "AoA"��
                if ((S_Va_target_TypeTurb(ns)>0) && (S_Va_target_TypeDA(ns)==1))
                    
                    k=S_Va_target_TypeTurb(ns);
                    S_Va_cnt_One=S_Va_cnt_One+1;
                    
                    fid3=fopen([target_runPath,char(S_Va_turb_FilePath(k))],'r+');
                    for i=1:(S_Va_turb_ModiLine(k)-1)
                       fgetl(fid3);
                    end
                    fseek(fid3,0,'cof');
                    fprintf(fid3,char(S_Va_turb_ModiFormat(k)),S_Va_ana(S_Va_cnt_One,run_flag(kp)));                
                    fclose(fid3);
                    
                end
            end
            
            % 4. ���ü���CFL��
            fid4=fopen([target_runPath,'\Input\InputNS.dat'],'r+');
            for i=1:(rep_line_CFL-1)
               fgetl(fid4);
            end
            fseek(fid4,0,'cof');
            fprintf(fid4,rep_format_CFL,CFL_Renew);
            fclose(fid4);
            
            % 5. �������׼�����ȽϿ���
            fid5=fopen([target_runPath,'\Input\InputNS.dat'],'r+');
            for i=1:(rep_line_STD_Comp-1)
               fgetl(fid5);
            end
            fseek(fid5,0,'cof');
            fprintf(fid5,rep_format_STD_Comp,STD_Comp_Renew);
            fclose(fid5);
            
            % 6. ���ü�����������
            fid6=fopen([target_runPath,'\Input\InputNS.dat'],'r+');
            for i=1:(rep_line_Resi-1)
               fgetl(fid6);
            end
            fseek(fid6,0,'cof');
            fprintf(fid6,rep_format_Resi,Resi_Renew);
            fclose(fid6);
            
            % 7. ����L2NORM����ˮƽ
            fid7=fopen([target_runPath,'\Input\InputTurb.dat'],'r+');
            for i=1:(rep_line_L2NORM-1)
               fgetl(fid7);
            end
            fseek(fid7,0,'cof');
            fprintf(fid7,rep_format_L2NORM,L2NORM_Renew);
            fclose(fid7);
            
            % 8. �ָ�ָ���Ŷ�����Ϊ STD ֵ
            if (flag_Va_init==0)
                for ns=1:S_Va_turb_Num
                    % ����ʼ�Ŷ�������Ŀ��ͬ����������Ϊ��������e.g. "CKS"��
                    % �򽫸ñ����ָ�ԭʼ STD ֵ
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
            
            % !!! ��Ҫ�޸� !!! �����ļ�������Դ���� Resi_flag_label
            % 9. ���������������� NSData �� TMuData �ļ��������� Output �ļ����µ������ļ�
            % sourcePath������ CFD �����ļ�����·��
            % targetPath�������ϳ�Ա������·��
            % ��ȡ��ǰ��������ֵ�����ļ�·��
            if (p==1)
                pre_Path=[targetPath,'Pre_Ensemble_Generation','\',initPath,'_',num2str(Resi_flag_label(run_flag(kp)),'%03d')];
            else
                addPath=[targetPath,'Inner_Iteration_Round_',num2str(p-1,'%03d')];
                pre_Path=[addPath,'\Renew_Ensemble_Member_',num2str(Resi_flag_label(run_flag(kp)),'%03d')];
            end
            % ���¹������� NSData �ļ�
            pre_FieldPath=[pre_Path,'\NSData.PLT'];
            field_val=load(pre_FieldPath);
            save([sourcePath,'\NSData.PLT'],'field_val','-ascii');
            
            % % ���µ�ǰ�������ݣ���������ͬ������������
            % TMu_new=S_Va_ana(((PT+1):MDA),run_flag(kp));
            % save([sourcePath,'\TMuData.PLT'],'TMu_new','-ascii');
            
            S_Va_cnt_Two=0;
            SS=PT+1; % �����������
            ST=PT; % ���������յ�
            for i=1:S_Va_target_Num
                if (S_Va_target_TypeDA(i)==2)
                    
                    S_Va_cnt_Two=S_Va_cnt_Two+1;

                    % �������������յ㣬��������
                    SL=S_Va_Flow_length(S_Va_cnt_Two);
                    ST=ST+SL;
                    
                    if (strcmp(char(S_Va_target_Name(i)),'Vc')==1)
                        Vc_new=S_Va_ana((SS:ST),run_flag(kp));
                        save([sourcePath,'\VcData.PLT'],'Vc_new','-ascii');
                        SS=ST+1;
                    elseif (strcmp(char(S_Va_target_Name(i)),'TMu')==1)
                        TMu_new=S_Va_ana((SS:ST),run_flag(kp));
                        save([sourcePath,'\TMuData.PLT'],'TMu_new','-ascii');
                        SS=ST+1;
                    elseif (strcmp(char(S_Va_target_Name(i)),'TMuP')==1)
                        TMuP_new=S_Va_ana((SS:ST),run_flag(kp));
                        save([sourcePath,'\TMuPData.PLT'],'TMuP_new','-ascii');
                        SS=ST+1;
                    elseif (strcmp(char(S_Va_target_Name(i)),'Beta')==1)
                        Beta_new=S_Va_ana((SS:ST),run_flag(kp));
                        save([sourcePath,'\Beta_plus.PLT'],'Beta_new','-ascii');
                        SS=ST+1;
                    end
                    
                end
            end

            % 10. ���� CFD �������ɳ�ʼ���ϳ�Ա����������
            % ����CFD���������
            fprintf(fid_all,'**************************************************\n');
            fprintf(fid_all,'|                 PROGRAM STATUS                  \n');
            fprintf(fid_all,'**************************************************\n');
            fprintf(fid_all,'|      INNER ITERATION PROCEDURE Round %03d       \n',p);
            fprintf(fid_all,'**************************************************\n');
            fprintf(fid_all,'|                                                 \n');
            fprintf(fid_all,'|    The calculation of %03d sample is RUNNING!   \n',run_flag(kp));
            fprintf(fid_all,'|                                                 \n');
            fprintf(fid_all,'**************************************************\n');
            fprintf(fid_all,'\n');
                
            fprintf(['���ϳ�Ա�ڵ����� ',num2str(p,'%03d'),' �֣�״̬ ',num2str(run_flag(kp),'%03d')]);
            fprintf(fid_info,['���ϳ�Ա�ڵ����� ',num2str(p,'%03d'),' �֣�״̬ ',num2str(run_flag(kp),'%03d')]);
            % ���ʱ����ʾ Ŀ��ͬ�������еĵ�����
            S_Va_cnt_One=0;
            for ns=1:S_Va_target_Num
                % ��Ŀ��ͬ������Ϊ��ʼ�Ŷ���������Ϊ��ͬ��������e.g. "Ma" "AoA"��
                if ((S_Va_target_TypeTurb(ns)>0) && (S_Va_target_TypeDA(ns)==1))
                    S_Va_cnt_One=S_Va_cnt_One+1;
                    fprintf([' ',char(S_Va_target_Name(ns)),' = %.5f'],S_Va_ana(S_Va_cnt_One,run_flag(kp)));
                    fprintf(fid_info,[' ',char(S_Va_target_Name(ns)),' = %.5f'],S_Va_ana(S_Va_cnt_One,run_flag(kp)));
                end
            end
            fprintf(' ����������\n');
            fprintf(fid_info,' ����������\n');
            
            % ����CFD�����������д��ڣ��������г���
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
            if (JX_Turb(2)==0)
                Resi_All=load([target_runPath,'\Output\Residual_1.PLT']);
            elseif ((JX_Turb(2)==1) || (JX_Turb(2)==2))
                Resi_All=load([target_runPath,'\Output\Residual_1_Turb_Const.PLT']);
            elseif (JX_Turb(2)==3)
                Resi_All=load([target_runPath,'\Output\Residual_1_Turb_Beta.PLT']);
            end
            % Resi_All=load([target_runPath,'\Output\Residual_1_Turb_Const.PLT']);
            Resi_Final=Resi_All(size(Resi_All,1),2);
            L2Norm_All=load([target_runPath,'\Output\L2NORM.PLT']);
            L2Norm_Init=L2Norm_All(1,(6:9));
            L2Norm_Final=L2Norm_All((size(L2Norm_All,1)-1),(6:9));
            L2Norm_Ratio=100*abs(L2Norm_Final./L2Norm_Init);
            L2Norm_Max=max(L2Norm_Ratio);
            L2Norm_Min=min(L2Norm_Ratio);
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
            fprintf(fid_all,'|      INNER ITERATION PROCEDURE Round %03d       \n',p);
            fprintf(fid_all,'**************************************************\n');
            fprintf(fid_all,'|                                                 \n');
            fprintf(fid_all,'|   The calculation of %03d sample is FINISHED!   \n',run_flag(kp));
            fprintf(fid_all,'|                                                 \n');
            fprintf(fid_all,'**************************************************\n');
            fprintf(fid_all,'\n');
            fprintf(fid_all,'\n');
            
            fprintf(['���ϳ�Ա�ڵ����� ',num2str(p,'%03d'),' �֣�״̬ ',num2str(run_flag(kp),'%03d'),' �ѱ�����ϣ�Total Step = %d��Final Residual = %.5f��Min L2Norm = %.2f%%��Max L2Norm = %.2f%%\n'],(size(Resi_All,1)-1),Resi_Final,L2Norm_Min,L2Norm_Max);
            fprintf(fid_info,['���ϳ�Ա�ڵ����� ',num2str(p,'%03d'),' �֣�״̬ ',num2str(run_flag(kp),'%03d'),' �ѱ�����ϣ�Total Step = %d��Final Residual = %.5f��Min L2Norm = %.2f%%��Max L2Norm = %.2f%%\n'],(size(Resi_All,1)-1),Resi_Final,L2Norm_Min,L2Norm_Max);
            
            % ��Ǽ���ﵽ������������� Resi_flag = 1
            if (Resi_Final<=Resi_target_lev(2))
                Resi_renew_flag(run_flag(kp))=1; % ����״̬��Ϊ 1��������״̬��Ϊ 0
                Resi_renew_label_pre=Resi_renew_label_pre+1; % ����״̬�� +1
                Resi_renew_flag_label(Resi_renew_label_pre)=run_flag(kp); % ���浱ǰ����״̬��� Resi_flag_label���ۼ��� = Resi_label_pre������ = N_new
            end
            
            if (strcmp(char(S_Va_ref_Name(1)),'BL')==1)
                
                % ������ȡ�߽����Ϣ�������ļ�
                fprintf(['���ϳ�Ա�ڵ����� ',num2str(p,'%03d'),' �֣�״̬ ',num2str(run_flag(kp),'%03d'),' �߽��������ȡ��ʼ\n']);
                fprintf(fid_info,['���ϳ�Ա�ڵ����� ',num2str(p,'%03d'),' �֣�״̬ ',num2str(run_flag(kp),'%03d'),' �߽��������ȡ��ʼ\n']);
                [s,e]=dos(['cd ',target_runPath,' && run_BL_Hump_Extract_Tecplot.bat']);
                fprintf(['���ϳ�Ա�ڵ����� ',num2str(p,'%03d'),' �֣�״̬ ',num2str(run_flag(kp),'%03d'),' �߽��������ȡ���\n']);
                fprintf(fid_info,['���ϳ�Ա�ڵ����� ',num2str(p,'%03d'),' �֣�״̬ ',num2str(run_flag(kp),'%03d'),' �߽��������ȡ���\n']);
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
                    BLSC=zeros(length(BL_Cal_Vel),N_pre); % !!!
                    BLSC(:,run_flag(kp))=BL_Cal_Vel;
                else
                    BLSC(:,run_flag(kp))=BL_Cal_Vel;
                end
                
            end
            
            % ���ӣ�����ļ�ÿһ�ڵ����������һ�ļ���
            savePath=[addPath_pre,'\Renew_Ensemble_Member_',num2str(run_flag(kp),'%03d')];
            fileList=dir(sourcePath);
            mkdir(savePath);
            filename_s=cell(length(fileList),1);
            filename_t=cell(length(fileList),1);

            for i=3:length(fileList)
                  filename_s{i}=[sourcePath,'\',fileList(i).name];
                  filename_t{i}=savePath;
                  copyfile(filename_s{i},filename_t{i});
            end

%             if (strcmp(char(S_Va_ref_Name(1)),'BL')==1)
%                 % ��ȡ�߽���ٶ�������
%                 cnt_BL_cal=0;
%                 BL_Cal_Vel=zeros(kb,1);
%                 for i=1:N_BL
%                     pre_BL=[savePath,'\BL_Output\BL_Hump_xc_',num2str(BL_lab(i),'%03d'),'.dat'];
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
            pre_Cp_Cf_yplus_Path=[savePath,'\CP_CF_Comp_New.PLT'];
            pre_ErrInfo=[savePath,'\FlowField_Err_Info.txt'];

            % �����ڵ���Ԥ�ⲽ��ɺ�Ĺ۲����ݣ�e.g. "Cp"��
            Cp_Cf_yplus_S=load(pre_Cp_Cf_yplus_Path);
            Err_val_all=load(pre_ErrInfo);
            
            Xp=Cp_Cf_yplus_S(:,1);
            Yp=Cp_Cf_yplus_S(:,2);
            CpS=Cp_Cf_yplus_S(:,3);
            CfS=Cp_Cf_yplus_S(:,4);
            yplusS=Cp_Cf_yplus_S(:,5);
            
            TMu_MSE=Err_val_all(1,7);
            TMu_RMSE=Err_val_all(2,7);
            % �����������ż���������ָ��
            if (Resi_renew_flag(run_flag(kp))==1)
                TMu_RMSE_sum=TMu_RMSE_sum+TMu_RMSE;
            end
            
            if (sum(cnt_flag)==0)
                CpSC=zeros(length(CpS),N_pre);
                CfSC=zeros(length(CfS),N_pre);
                yplusSC=zeros(length(yplusS),N_pre);
                
                CpSC(:,run_flag(kp))=CpS;
                CfSC(:,run_flag(kp))=CfS;
                yplusSC(:,run_flag(kp))=yplusS;
            else            
                CpSC(:,run_flag(kp))=CpS;
                CfSC(:,run_flag(kp))=CfS;
                yplusSC(:,run_flag(kp))=yplusS;
            end

            % -------------------- ���н������ END -------------------- %
            
            run_flag(kp)=0;
            delete([target_runPath,'\ProgramExitFlag.txt']);
            cnt_flag(kp)=cnt_flag(kp)+1; % �ۼƸ���������
            stop_flag(kp)=0;
            pause(3); % �����ļ�����ͣ����
            
        end
        
        if (sum(stop_flag)==PPNS)
            % �Ҳ������е��������MATLAB������ͣɨ������
            pause(10);
        end
        
        if (sum(stop_flag)==0) && (sum(cnt_flag)==N_pre)
            % ������������У��������˳�����������ɲ��г���������
            exit_flag=1;
            pause(10);
            break;
        end
        
    end
    
end

% ����״̬���������������󣬼����µļ��ϳ�Ա���� N_new
N_new=sum(Resi_renew_flag);
Resi_renew_flag_label((N_new+1):end)=[];
Resi_flag_label=Resi_renew_flag_label; % ���²�����ȫ�������������� Resi_flag_label�����������һ���ڵ�������
% ������������������ָ��
Val_Sen=TMu_RMSE_sum/N_new;

% ��ȡ���º��ͬ���۲����
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
        elseif (strcmp(char(S_Va_ref_Name(i)),'BL')==1)
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
        elseif (strcmp(char(S_Va_ref_Name(i)),'BL')==1)
            S_Va_init_Part_Three_F=[S_Va_init_Part_Three,BLSC'];
            S_Va_init_Part_Three=S_Va_init_Part_Three_F;
        end
    end
end
if (S_Va_cnt_Three==0)
    S_Va_init_Part_Three=[];
end

% �������������Ϣ
save([addPath_pre,'\Com_Ensem_Gen.dat'],'S_Va_init_Part_Three','-ascii');

if ((strcmp(char(S_Va_ref_Name(1)),'Cp')==1) || (strcmp(char(S_Va_ref_Name(1)),'Cf')==1) || (strcmp(char(S_Va_ref_Name(1)),'yplus')==1))

    % !!!!! �޸Ĳ��֣���ȡ��Ч���ֵĹ۲���Ϣ !!!!!
    S_Va_init_Part_Three_eff=zeros(N_new,size(S_Va_init_Part_Three,2));
    n_combined_eff=0;
    for i=1:N_pre
        if (Resi_renew_flag(i)==1)
            n_combined_eff=n_combined_eff+1;
            S_Va_init_Part_Three_eff(n_combined_eff,:)=S_Va_init_Part_Three(i,:);
        end
    end

    % !!!!! �޸Ĳ��� !!!!!
    % �����м���״̬����õ���Ŀ��ͬ���������ֲ����İ����߽��л���
    S_Va_init_Part_Three_min=min(S_Va_init_Part_Three_eff);
    S_Va_init_Part_Three_max=max(S_Va_init_Part_Three_eff);
    S_Va_init_Part_Three_mid=(S_Va_init_Part_Three_min+S_Va_init_Part_Three_max)/2; % ���°���������
    target_width=S_Va_init_Part_Three_max-S_Va_init_Part_Three_min; % ���°���߽���

    if (flag_flow_type==0)

        y1=S_Va_init_Part_Three_max(1:shape_lower_num);
        y2=S_Va_init_Part_Three_min(1:shape_lower_num);
        y3=S_Va_init_Part_Three_max(shape_lower_num+1:end);
        y4=S_Va_init_Part_Three_min(shape_lower_num+1:end);
        % �½����ӻ�ͼ��
        % !!! ���Ӱ����������� !!! Middle_line ������
        figure;
        hold on
        [hha1,hhb1,hhc1]=shadedplot(Xp(1:shape_lower_num),y1,y2,'[0.8706 0.9216 0.9804]');
        [hha2,hhb2,hhc2]=shadedplot(Xp(shape_lower_num+1:end),y3,y4,'[0.8706 0.9216 0.9804]');
        hh1=plot(shape_combined_arr_init(:,2),shape_combined_arr_init(:,7),'ro','MarkerSize',4.0); % ���Ʊ�׼ʵ��ֵ
        hh2=plot(shape_combined_arr_init(:,2),S_Va_init_Part_Three_mid,'-','LineWidth',0.8,'Color',[0.7451,0.7451,0.7451]); % ���ư�����������
        % plot(Xp,S_Va_init_Part_Three_max); % ���������������
        % plot(Xp,S_Va_init_Part_Three_min); % ���������������
        legend([hhb1,hhc1,hha1(2),hh2,hh1],{'���������������','���������������','���±߽���Χ����','���������������','ʵ��ֵ����ֵ'}); % ��Ӳ���ͼ��
        xlabel('\itx'),ylabel('\itCp'),title(['�ڵ����� ',num2str(p,'%03d'),' �ּ�����']);
        grid on
        box on
        hold off
        % ����ɱ༭ͼƬ
        saveas(gcf,['.\Output_Figures\',savename,'\Renew_Round_',num2str(p,'%03d'),'_Output.fig']);

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
        xlabel('\itx'),ylabel('\itCp'),title(['�ڵ����� ',num2str(p,'%03d'),' �ּ�����']);
        grid on
        box on
        hold off
        % ����ɱ༭ͼƬ
        saveas(gcf,['.\Output_Figures\',savename,'\Renew_Round_',num2str(p,'%03d'),'_Output.fig']);

    end
    
end

% �����ʼ���º�״̬���������ɺ���������������
if (flag_da_relative==1)
    S_Va_pre=[(S_Va_ana((1:MDA),:))',S_Va_init_Part_Three]';
elseif (flag_da_relative==0)
    % !!!! 2021/03/03ģ����� !!!!
    % �Թ۲������ֵ����ת����Ŀ���ǵõ���ʵ��ֵ����ֵ����������ֵ
    % �ɽ�ʵ��ֵ����ֵ������Ϊ��
    % ͬ��Ŀ�꼴ת��Ϊ��С�۲�ֵ��ʵ��ֵ����ֵ���������������
    S_Va_init_Part_Three_err=(S_Va_init_Part_Three./D_obs_e')-ones(size(D_obs_e'));
    S_Va_pre=[(S_Va_ana((1:MDA),:))',S_Va_init_Part_Three_err]';
end

% !!!!! �޸Ĳ��֣���ȡ��Ч���ֵ�״̬���� !!!!!
S_Va_pre_eff=zeros(size(S_Va_pre,1),N_new);
n_combined_eff=0;
for i=1:N_pre
    if (Resi_renew_flag(i)==1)
        n_combined_eff=n_combined_eff+1;
        S_Va_pre_eff(:,n_combined_eff)=S_Va_pre(:,i);
    end
end

% -- ���ݲ��� -- %
save([new_output_folder,savename,'_S_Va_pre_renew_',num2str(p,'%03d'),'.mat'],'S_Va_pre'); % ��������ԭʼ��Ա�����������
% -- ���ݲ��� -- %

% --------------- �޸Ĳ���Ԥ��+���ϳ�Ա���� END --------------- %

fprintf('\n');
fprintf(fid_info,'\n');

end
