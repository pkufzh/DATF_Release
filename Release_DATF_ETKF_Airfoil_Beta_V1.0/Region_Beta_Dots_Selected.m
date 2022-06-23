
%% 更新时间：2021/04/02
%% 函数功能：对所选定区域内的网格单元进行平面随机采样
%% 满足规则：点数最大为XY_dots_num_max，且任意两点之间距离不小于XY_dots_dis_lim
function [XY_cnt,Xp_sca,Yp_sca,XY_sel_flag]=Region_Beta_Dots_Selected(NDAG_coor,DA_GN_x,DA_GN_y,XY_dots_num_max,XY_dots_dis_lim)

    % 初始化已采样点数
    XY_cnt=1;
    % 设置标记是否选中的集合： 1 - 选中；0 - 未选中
    % XY_sel_flag=zeros(NDAG_coor,1);
    % 选取集合中任意一点为第一个点
    kf=randi(NDAG_coor);
    XY_sel_flag(XY_cnt)=kf;
    % 保存首个采样点的坐标信息
    Xp_sca(XY_cnt)=DA_GN_x(kf);
    Yp_sca(XY_cnt)=DA_GN_y(kf);

    % 设置采样停止标志，初始为0，若变为1，则退出采样过程，否则重复下一次有效采样
    stop_flag=0;

    % 初始化每轮扫描次数（直至成功寻找到目标点）
    cnt_scan=0;
    % 设置每轮最大扫描次数，若超过阈值，则将 stop_flag 设为 1，退出采样过程
    cnt_scan_lim=round(NDAG_coor/2);

    % 采样过程启动
    while ((XY_cnt<XY_dots_num_max) && (stop_flag==0))

        % 当前随机选取编号（正整数范围为[1,NDAG_coor]）
        k=randi(NDAG_coor);
        % 扫描次数加一
        cnt_scan=cnt_scan+1;

        % 如果当前编号没有被选过，则继续判断距离条件
        if (sum(find(XY_sel_flag==k))==0)
            % 记录当前选取编号对应网格单元的坐标
            temp_x=DA_GN_x(k);
            temp_y=DA_GN_y(k);
            % 初始化距离矩阵，用以条件判断
            temp_dis=zeros(XY_cnt,1);
            for i=1:XY_cnt
                temp_dis(i)=sqrt(((temp_x-Xp_sca(i))^2)+((temp_y-Yp_sca(i))^2));
            end
            % 判断任意两点距离是否满足条件（不小于XY_dots_dis_lim）
            % 若满足条件
            if (sum(find(temp_dis<=XY_dots_dis_lim))==0)
                % 满足条件点加一
                XY_cnt=XY_cnt+1;
                % 保存新加入的采样点
                Xp_sca(XY_cnt)=DA_GN_x(k);
                Yp_sca(XY_cnt)=DA_GN_y(k);
                % 标记已完成采样点编号，避免重复采样
                XY_sel_flag(XY_cnt)=k;
                % 初始化扫描次数
                cnt_scan=0;
            end
        end
        
        % 若当前扫描次数超过阈值，则将 stop_flag 设为 1，退出采样过程
        if (cnt_scan>cnt_scan_lim)
            stop_flag=1;
        end

    end
    
    Xp_sca=Xp_sca';
    Yp_sca=Yp_sca';
    XY_sel_flag=XY_sel_flag';
    
end
