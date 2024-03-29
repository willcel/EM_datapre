%% 存储电流的转折点

% clc
% clear
% close all

% ns = 17;

point1set = [];
point2set = [];
point3set = [];
point4set = [];

% load('data_avg.mat')
% load('current_avg.mat') 
load('data_avg2.mat')  % 处理零点偏移问题
load('current_avg2.mat') 


% fs = 10^6;
dt = 1/fs;

time = (1:size(data_avg_all, 2)).*dt;

% figure
for k=1:ns

    data_receive_coil = data_avg_all(k, :);
    current = current_avg_all(k, :);
    
    % peak point
    [max_val, ind] = max(current);
     index_peak = min(find(abs(current-max_val)<0.1));

    % end point
    ind1 = find(current<0);
    set2 = ind:length(current);
    point2 = min(intersect(set2, ind1));
%     point2 = find(time>3.0508e-3,1);

%     if(k==4 || k==5 ||  k ==6 || k==8 || k==9 || k==17)
%         point2 = 2128;
%     end


    % 衰减起始点
    diff_current = diff(current);
    ind2 = find(diff_current<-4);
    point3 = min(intersect(set2, ind2));


    % start point
    set1 = 1:ind;   % 电流到达峰值的前半段
    ind11 = find(diff_current>1.5);
    point1 = min(intersect(set1,ind11));
    if(isempty(point1))
        point1 = 1;
    end
    if(0)

        kk = k; %7;
        data_receive_coil = data_avg_all(kk, :);
        current = current_avg_all(kk, :);

%         figure(Position=[242.33333333333	151	1106.00000000000	643.333333333333])
%         subplot(2,1,1)
%         semilogy(time, data_receive_coil)
%         plot(time, data_receive_coil, 'Color', [1 0.5 0], 'LineWidth',1.5)
%         xlabel('time(s)')
%         ylabel('voltage (V)')
%         xlim([0.0,0.005])
% %         xlim([0.02,0.03])
% %         ylim([-1,1]*1e-5)
%         title(['测点',num2str(kk)])
%         grid on
%         set(gca,'FontSize',16,'FontWeight','bold')
%         text(-0.1,1.05,'(a)','Units','normalized','FontSize',22)
% 
%         figure
% %         subplot(2,1,2)
%         plot(time, current, 'Color', [0 1 1], 'LineWidth',1.5)
%         hold on
%         scatter(time(index_peak), current(index_peak),'ro')
%         hold on
%         scatter(time(point1), current(point1),'bo')
%         hold on
%         scatter(time(point2), current(point2),'go')
%         hold on
%         scatter(time(point3), current(point3),'ko')
%         xlim([0,0.005])
% %         ylim([-20 200])
% 
%         xlabel('time(s)')
%         ylabel('current (A)')
%         grid on
%         set(gca,'FontSize',16,'FontWeight','bold')
%         text(-0.1,1.05,'(b)','Units','normalized','FontSize',22)
        
        
        figure(Position=[242.33333333333	151	1106.00000000000	643.333333333333])
        sig1 = data_receive_coil;
       
        sig2 = -1*ones(1, length(sig1));
        idx = find(sig1 < 0);
        sig2(idx) = - sig1(idx);

        semilogy(time, (sig1))
        hold on
        semilogy(time, sig2, '-r')
        legend('正值', '负值')
                xlabel('time(s)')
        ylabel('voltage (V)')
                title(['测点',num2str(kk)])
        grid on
        set(gca,'FontSize',16,'FontWeight','bold')
        xlim([0.0,0.03])
        1;
        ylim([1e-10 1])
        %%
        name1 = ['ns',num2str(k),'_measure',num2str(k)];
        svfig(name1, '.\rawVolt\step2')
        close all
        
    end
    
    

    tt1 = time(index_peak)-time(point1);
    tt2 = time(point3)-time(index_peak);    % 中间段
    tt3 = time(point2)-time(point3);
    
    pulse_width = time(point2)-time(point1);
    
%     pulse_width*1000
    
    point1set = [point1set; time(point1); current(point1)];
    point2set = [point2set; time(index_peak); current(index_peak)];
    point3set = [point3set; time(point3); current(point3)];
    point4set = [point4set; time(point2); current(point2)];
    
end

%%
% point4set( 27 * 2 - 1)
% point4set(1:2:2*ns-1) = 2.18e-3;
%%
save('point1set.txt','point1set','-ascii')
save('point2set.txt','point2set','-ascii')
save('point3set.txt','point3set','-ascii')
save('point4set.txt','point4set','-ascii')
