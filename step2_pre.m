%% 存储电流的转折点

% clc
% clear
% close all


point1set = [];
point2set = [];
point3set = [];
point4set = [];

load('data_avg.mat')
load('current_avg.mat') 
% current_avg_all = -current_avg_all; save('current_avg.mat','current_avg_all') % 不知为什么，工程解决
%{
    data_avg_all(4,:) = 0.5*( data_avg_all(3,:)+data_avg_all(5,:));
    current_avg_all(4,:) = 0.5*( current_avg_all(3,:)+current_avg_all(5,:));

    save data_avg.mat data_avg_all
    save current_avg.mat current_avg_all
%}

% fs = 10^6;
dt = 1/fs;

% time = (1:size(data_avg_all, 2)).*dt;

point1bank = [];
% figure
for k=1:ns

    data_receive_coil = data_avg_all(k, :);
    current = current_avg_all(k, :);

    % peak point
    [max_val, ind] = max(current);

    diff_current = diff(current);
    % start point
    set1 = 1:ind;   % 电流到达峰值的前半段
    ind11 = find(diff_current>1.5);
    point1 = min(intersect(set1,ind11));

%     point1 = 1;

    point1bank = [point1bank; point1];

    
end

minpoint = min(point1bank)+1;

data_avg_all = data_avg_all(:,minpoint:end);
current_avg_all = current_avg_all(:,minpoint:end);

save('data_avg2.mat','data_avg_all')
save('current_avg2.mat','current_avg_all') 


%%

        if(0)
            %%
            signal = data_receive_coil;
            dt = 1/fs;
            time = (1:length(signal)).*dt;
            
            figure
            subplot(2,1,1)
            plot(time*1000, signal,'LineWidth',2.0)
            set(gca,'FontSize',12,'FontWeight','bold')
            xlabel('time(ms)')
            ylabel('voltage (V)')
            grid on
            xlim([0 5])
            title('叠加')
%             ylim([0,3]*1e-4)
                
            subplot(2,1,2)
            plot(time*1000, current,'LineWidth',2.0)                
            xlim([0,5])
            xlabel('time(ms)')
            ylabel('current (A)')
            grid on
            set(gca,'FontSize',12,'FontWeight','bold')
        end
