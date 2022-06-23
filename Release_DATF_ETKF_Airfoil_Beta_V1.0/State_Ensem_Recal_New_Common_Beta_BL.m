
% �������ڣ�2021/05/20
% �������ܣ�ʵ������ͬ�������ؼ���
% �����������˵��
% ���������p ���� ��ǰ�ڵ����ִ�
%           S_Va_ana_ave ���� ��ǰԤ��õ���״̬��������
%           MDA ���� ͬ����������ά��
%           MStepNS_recal ����  �����ؼ�����㲽��
%           targetPath ���� �ڵ����������Ŀ��·��
% �����������

function State_Ensem_Recal_New_Common_Beta_BL(p,S_Va_ana_ave,MStepNS_recal,targetPath)

% ����ȫ�ֱ���
global savename new_output_folder fid_all fid_info PPNS initPath CFL_Setting Resi_target Resi_target_lev JX_NS STD_Comp JX_Turb L2NORM_Setting flag_Va_init;
global S_Va_turb_Name S_Va_turb_Type S_Va_turb_Num S_Va_turb_FilePath S_Va_turb_ModiLine S_Va_turb_ModiFormat S_Va_turb_Std;
global S_Va_target_Name S_Va_target_Num S_Va_target_TypeTurb S_Va_target_TypeDA S_Va_turb_target S_Va_Flow_length PT PS;
global S_Va_ref_Name S_Va_ref_Num;
global shape_upper shape_lower shape_outline shape_combined_arr_init Va_Beta Va_NNP_Beta Va_NNP_Beta_new;
global flag_da_relative flag_sigma_obs_adapt flag_sigma_obs_coeff D_obs_e_ori D_obs_e flag_flow_type;
global BL_setting N_BL L_BL BL_num BL_lab BL_dot BL_obs_num kb;
global Resi_flag Resi_flag_label;

% �޸�NS����������������MStepNS -- NS��������������
% ָ���滻�У���Ӧ�����ļ���InputNS.dat����9��
rep_line_MStepNS=9;
rep_format_MStepNS='%-16d';

% �޸����㿪�أ� JX_NS -- �Ƿ���������������㣨0, 1��
%                JX_Turb -- �Ƿ���������������㣨0, 1, 3��
% ָ���滻�У��ֱ��Ӧ�����ļ���InputNS.dat����10�У���InputTurb.dat����7��
rep_line_JX_NS=10;
rep_format_JX_NS='%-12d';
rep_line_JX_Turb=7;
rep_format_JX_Turb='%-12d';

% �޸ļ���CFL����CFL_Setting
% ָ���滻�У���Ӧ�����ļ���InputNS.dat����8��
rep_line_CFL=8;
rep_format_CFL='%-12d';
CFL_Recal=CFL_Setting(3);

% �޸ı�׼�����������ȽϿ��أ�STD_Comp
% ָ���滻�У���Ӧ�����ļ���InputNS.dat����22��
rep_line_STD_Comp=22;
rep_format_STD_Comp='%-12d';
STD_Comp_Recal=STD_Comp(3);

% �޸������в Resi_target
rep_line_Resi=15;
rep_format_Resi='%-10.2e';
Resi_Recal=Resi_target(3);

% �޸����� L2NORM ˮƽ�� L2NORM_Setting
rep_line_L2NORM=8;
rep_format_L2NORM='%-12d';
L2NORM_Recal=L2NORM_Setting(3);

target_runPath='.\Program_CFD_Parallel_Recal';

% 1. ���£��رգ�������ݿ���
% JX_NS=0; % !!!!!!!!!!! ���ǿ���ǰһ���ؼ���������
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

% 2. ����NS����������е�������
fid3=fopen([target_runPath,'\Input\InputNS.dat'],'r+');
for i=1:(rep_line_MStepNS-1)
   fgetl(fid3);
end
fseek(fid3,0,'cof');
fprintf(fid3,rep_format_MStepNS,MStepNS_recal);
fclose(fid3);

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
        fprintf(fid3,char(S_Va_turb_ModiFormat(k)),S_Va_ana_ave(S_Va_cnt_One));                
        fclose(fid3);

    end
end

% 4. ���ü���CFL��
fid4=fopen([target_runPath,'\Input\InputNS.dat'],'r+');
for i=1:(rep_line_CFL-1)
   fgetl(fid4);
end
fseek(fid4,0,'cof');
fprintf(fid4,rep_format_CFL,CFL_Recal);
fclose(fid4);

% 5. �������׼�����ȽϿ���
fid5=fopen([target_runPath,'\Input\InputNS.dat'],'r+');
for i=1:(rep_line_STD_Comp-1)
   fgetl(fid5);
end
fseek(fid5,0,'cof');
fprintf(fid5,rep_format_STD_Comp,STD_Comp_Recal);
fclose(fid5);
            
% 6. ���ü�����������
fid6=fopen([target_runPath,'\Input\InputNS.dat'],'r+');
for i=1:(rep_line_Resi-1)
   fgetl(fid6);
end
fseek(fid6,0,'cof');
fprintf(fid6,rep_format_Resi,Resi_Recal);
fclose(fid6);

% 7. ����L2NORM����ˮƽ
fid7=fopen([target_runPath,'\Input\InputTurb.dat'],'r+');
for i=1:(rep_line_L2NORM-1)
   fgetl(fid7);
end
fseek(fid7,0,'cof');
fprintf(fid7,rep_format_L2NORM,L2NORM_Recal);
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

% 9. �����������������ļ�

% % �������̹������� TMuData �ļ��������� Output �ļ����µ������ļ�
% % ���µ�ǰ�������ݣ���������������
% TMu_new=S_Va_ana_ave((PT+1):MDA);
% % ���¹������� TMuData �ļ�
% save('.\Program_CFD_Parallel_Recal\Output\TMuData.PLT','TMu_new','-ascii');

S_Va_cnt_Two=0;
SS=PT+1; % �����������
ST=PT; % ���������յ㣨��ֵ���������������£�
for i=1:S_Va_target_Num
    if (S_Va_target_TypeDA(i)==2)

        S_Va_cnt_Two=S_Va_cnt_Two+1;

        % �������������յ㣬��������
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
            save('.\Program_CFD_Parallel_Recal\Output\Beta_plus.PLT','Beta_new','-ascii'); % ��ͬ���õ��ĵ�Ԫ������ֵ��Ϣд��ָ���ļ�
            save('.\Program_CFD_Parallel_Recal\Output\Beta_plus_NNP.PLT','Va_NNP_Beta_new','-ascii'); % ��ͬ������ֵ�õ��Ľڵ�������ֵ��Ϣд��ָ���ļ�
            SS=ST+1;
        end

    end
end

% ɾ�����г�����ɱ���ļ�
if (exist([target_runPath,'\ProgramExitFlag.txt'],'file'))~=0
    delete([target_runPath,'\ProgramExitFlag.txt']);
end

% 10. ���� CFD �������ɳ�ʼ���ϳ�Ա����������
% ����CFD���������
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

fprintf(['���ϳ�Ա�ؼ��� ',num2str(p,'%03d'),' �֣�״̬']);
fprintf(fid_info,['���ϳ�Ա�ؼ��� ',num2str(p,'%03d'),' �֣�״̬']);
% ���ʱ����ʾ Ŀ��ͬ�������еĵ�����
S_Va_cnt_One=0;
for ns=1:S_Va_target_Num
    % ��Ŀ��ͬ������Ϊ��ʼ�Ŷ���������Ϊ��ͬ��������e.g. "Ma" "AoA"��
    if ((S_Va_target_TypeTurb(ns)>0) && (S_Va_target_TypeDA(ns)==1))
        S_Va_cnt_One=S_Va_cnt_One+1;
        fprintf([' ',char(S_Va_target_Name(ns)),' = %.5f'],S_Va_ana_ave(S_Va_cnt_One));
        fprintf(fid_info,[' ',char(S_Va_target_Name(ns)),' = %.5f'],S_Va_ana_ave(S_Va_cnt_One));
    end
end
fprintf(' ����������\n');
fprintf(fid_info,' ����������\n');

% cmd=('.\Program_CFD_Parallel_Recal\CFD_2D.exe');
% open(cmd);
[s,e]=dos('cd Program_CFD_Parallel_Recal && run_CFD_2D.bat&');
pause(5); % ��ͣ����

% �����ص��������˳���־
exit_recal_flag=0;

while (exit_recal_flag~=1)

    if (exist([target_runPath,'\ProgramExitFlag.txt'],'file'))~=0
        
        % -------------------- ���н������ START ------------------- %

        % ��ȡ��ǰ״̬�������ղв�ֵ Resi(dual)_Final
        if (JX_Turb(3)==0)
            Resi_All=load([target_runPath,'\Output\Residual_1.PLT']);
        elseif ((JX_Turb(3)==1) || (JX_Turb(3)==2))
            Resi_All=load([target_runPath,'\Output\Residual_1_Turb_Const.PLT']);
        elseif (JX_Turb(3)==3)
            Resi_All=load([target_runPath,'\Output\Residual_1_Turb_Beta.PLT']);
        end
        % Resi_All=load([target_runPath,'\Output\Residual_1_Turb_Const.PLT']);
        Resi_Final=Resi_All(size(Resi_All,1),2);
        % �����������
        com_data=load('.\Program_CFD_Parallel_Recal\Output\CP_CF_Comp.PLT');
        % ָ��������������ݱ���·��
        fid_com_adj_path='.\Program_CFD_Parallel_Recal\Output\CP_CF_Comp_New.PLT';
        fid_tar_new_path='.\Program_CFD_Parallel_Recal\Output\CP_New.PLT';

        % �Լ���ѹ��ϵ�����߽ڵ��Ž���������ʵʱ�����ǰͬ���ִ�ѹ���ֲ�����
        % �������ɹ��̡��ڵ���������ȫ����Ա��ɼ���󱣴����߰�����
        if (flag_flow_type==0)
            [shape_upper_num,shape_lower_num,shape_upper_arr_sor,shape_lower_arr_sor,n_com]=DA_Va_Airfoil_Ref_Adjustment_and_Visulization(com_data,fid_com_adj_path,fid_tar_new_path);
            % �ϲ��������������ʵ��ֵ����ֵ��������
            shape_combined_arr_recal=[shape_lower_arr_sor;shape_upper_arr_sor]; % �������棬��������
        elseif (flag_flow_type==1)
            [shape_num,shape_arr_sor,n_com]=DA_Va_NonAirfoil_Ref_Adjustment_and_Visulization(com_data,fid_com_adj_path,fid_tar_new_path);
            shape_combined_arr_recal=shape_arr_sor;
        end
        
        % ��ʼ��׼״̬���������ӻ�
        figure;
        hold on
        plot(shape_combined_arr_init(:,2),shape_combined_arr_init(:,7),'ro','MarkerSize',4.0);
        plot(shape_combined_arr_recal(:,2),shape_combined_arr_recal(:,4),'b-','LineWidth',2.0);
        xlabel('\itx'),ylabel('\itCp'),title(['�ؼ��� ',num2str(p,'%03d'),' �ּ�����']);
        legend('ʵ�����ֵ','�ؼ�������ѹ��ϵ���ֲ�');
        grid on
        box on
        hold off
        % ����ɱ༭ͼƬ
        saveas(gcf,['.\Output_Figures\',savename,'\Recal_Round_',num2str(p,'%03d'),'_Output.fig']);
        
        % �����������ʾ
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
        
        fprintf(['���ϳ�Ա�ؼ��� ',num2str(p,'%03d'),' �֣�״̬�ѱ������, Total Step = %d��Final Residual = %.5f\n'],(size(Resi_All,1)-1),Resi_Final);
        fprintf(fid_info,['���ϳ�Ա�ؼ��� ',num2str(p,'%03d'),' �֣�״̬�ѱ������, Total Step = %d��Final Residual = %.5f\n'],(size(Resi_All,1)-1),Resi_Final);

        if (strcmp(char(S_Va_ref_Name(1)),'BL')==1)
            % ������ȡ�߽����Ϣ�������ļ�
            fprintf(['���ϳ�Ա�ؼ��� ',num2str(p,'%03d'),' �߽��������ȡ��ʼ\n']);
            fprintf(fid_info,['���ϳ�Ա�ؼ��� ',num2str(p,'%03d'),' �߽��������ȡ��ʼ\n']);
            % [s,e]=dos(['cd ',target_runPath,' && run_BL_Hump_Extract_Tecplot.bat']);
            [s,e]=dos('cd Program_CFD_Parallel_Recal && run_BL_Hump_Extract_Tecplot.bat');
            fprintf(['���ϳ�Ա�ؼ��� ',num2str(p,'%03d'),' �߽��������ȡ���\n']);
            fprintf(fid_info,['���ϳ�Ա�ؼ��� ',num2str(p,'%03d'),' �߽��������ȡ���\n']);
            pause(3); % ��ͣ����
        end
            
        % ����ļ���棬����Ԥ�ⲽ֮ǰ�ļ�����
        save_recal='.\Program_CFD_Parallel_Recal\Output';
        fileList=dir(save_recal);
        mkdir([targetPath,'\Recalculation_Round_',num2str(p,'%03d')]); % ���� mkdir() �����������ļ���
        filename_s=cell(length(fileList),1);
        filename_t=cell(length(fileList),1);

        for i=3:length(fileList)
              filename_s{i}=[save_recal,'\',fileList(i).name];
              filename_t{i}=[targetPath,'\Recalculation_Round_',num2str(p,'%03d')];
              copyfile(filename_s{i},filename_t{i});
        end

        fprintf('\n');
        fprintf(fid_info,'\n');

        % -------------------- ���н������ END -------------------- %
        
        exit_recal_flag=1;  % �˳���־��һ
        
    end
    
    pause(10);
    
end

end
