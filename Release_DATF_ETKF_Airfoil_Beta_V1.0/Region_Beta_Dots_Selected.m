
%% ����ʱ�䣺2021/04/02
%% �������ܣ�����ѡ�������ڵ�����Ԫ����ƽ���������
%% ������򣺵������ΪXY_dots_num_max������������֮����벻С��XY_dots_dis_lim
function [XY_cnt,Xp_sca,Yp_sca,XY_sel_flag]=Region_Beta_Dots_Selected(NDAG_coor,DA_GN_x,DA_GN_y,XY_dots_num_max,XY_dots_dis_lim)

    % ��ʼ���Ѳ�������
    XY_cnt=1;
    % ���ñ���Ƿ�ѡ�еļ��ϣ� 1 - ѡ�У�0 - δѡ��
    % XY_sel_flag=zeros(NDAG_coor,1);
    % ѡȡ����������һ��Ϊ��һ����
    kf=randi(NDAG_coor);
    XY_sel_flag(XY_cnt)=kf;
    % �����׸��������������Ϣ
    Xp_sca(XY_cnt)=DA_GN_x(kf);
    Yp_sca(XY_cnt)=DA_GN_y(kf);

    % ���ò���ֹͣ��־����ʼΪ0������Ϊ1�����˳��������̣������ظ���һ����Ч����
    stop_flag=0;

    % ��ʼ��ÿ��ɨ�������ֱ���ɹ�Ѱ�ҵ�Ŀ��㣩
    cnt_scan=0;
    % ����ÿ�����ɨ���������������ֵ���� stop_flag ��Ϊ 1���˳���������
    cnt_scan_lim=round(NDAG_coor/2);

    % ������������
    while ((XY_cnt<XY_dots_num_max) && (stop_flag==0))

        % ��ǰ���ѡȡ��ţ���������ΧΪ[1,NDAG_coor]��
        k=randi(NDAG_coor);
        % ɨ�������һ
        cnt_scan=cnt_scan+1;

        % �����ǰ���û�б�ѡ����������жϾ�������
        if (sum(find(XY_sel_flag==k))==0)
            % ��¼��ǰѡȡ��Ŷ�Ӧ����Ԫ������
            temp_x=DA_GN_x(k);
            temp_y=DA_GN_y(k);
            % ��ʼ������������������ж�
            temp_dis=zeros(XY_cnt,1);
            for i=1:XY_cnt
                temp_dis(i)=sqrt(((temp_x-Xp_sca(i))^2)+((temp_y-Yp_sca(i))^2));
            end
            % �ж�������������Ƿ�������������С��XY_dots_dis_lim��
            % ����������
            if (sum(find(temp_dis<=XY_dots_dis_lim))==0)
                % �����������һ
                XY_cnt=XY_cnt+1;
                % �����¼���Ĳ�����
                Xp_sca(XY_cnt)=DA_GN_x(k);
                Yp_sca(XY_cnt)=DA_GN_y(k);
                % �������ɲ������ţ������ظ�����
                XY_sel_flag(XY_cnt)=k;
                % ��ʼ��ɨ�����
                cnt_scan=0;
            end
        end
        
        % ����ǰɨ�����������ֵ���� stop_flag ��Ϊ 1���˳���������
        if (cnt_scan>cnt_scan_lim)
            stop_flag=1;
        end

    end
    
    Xp_sca=Xp_sca';
    Yp_sca=Yp_sca';
    XY_sel_flag=XY_sel_flag';
    
end
