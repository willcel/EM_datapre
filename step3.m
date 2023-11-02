%% 均值滤波处理程序 2ms


load('point1set.txt')
load('point4set.txt')
% load('data_avg.mat')
load('data_avg2.mat')
% data_avg_all(4,:) = 0.5*( data_avg_all(3,:)+data_avg_all(5,:));

dt = 1/fs;



time = (1:size(data_avg_all, 2)).*dt;


start_point = mean(point1set(1:2:end));
disp(['start_point = ', num2str(start_point*1000)])

% 测点16，17的电流转折点point4set不太对
% end_point = mean(point4set(1:2:15*2-1));

end_point = mean(point4set(1:2:end));
disp(['end_point = ', num2str(end_point*1000)])

[~, start_index] =min(abs(time-start_point));  %
[~, end_index] =min(abs(time-end_point));    %

time = time-start_point;  % 时间从0开始

index = end_index+1;             % 必须大于end point % time(index) 彰显起始时刻


pure_sec_field = (data_avg_all(:,index:end));
% pure_sec_field = abs(data_avg_all(:,index:end));

ind_neg_final = size(data_avg_all, 2)-index+1;

filtered_sec_field = [];

%{ 
%没有滤波的二次场信号
close all
t = time(index:end);    % 这里t的起始值就是抽道的起始时刻
t = t(1:ind_neg_final);
t1 = t*1000;
figure('Position',[10	68.333333333333	1331	729.666666666667])
for i=1:1
    loglog(t1, pure_sec_field(i,:))
    hold on
%     legend_str1{i} = ['测点',num2str(i)];
end
xlim([0,20])
xlabel('time (ms)')
ylabel('Voltage (V)')
% legend(legend_str1,'NumColumns',4)
grid on
set(gca,'FontSize',16,'FontWeight','bold')
xlim([0 10])
%}

%%
for i=1:ns
    filter_signal = [];
    t = time(index:end);    % 这里t的起始值就是抽道的起始时刻
    t = t(1:ind_neg_final);
    
    t1 = t*1000;
    
    % low pass filter 30kHz
    sfx1 = lowpass(data_avg_all(i,:),30000,fs,'ImpulseResponse','fir','Steepness',0.95);
    sfx2 = sfx1(:,index:end);
    sfx3(i,:) = sfx2;
    ind_neg = find(sfx2<0); 
    ind_paint = find(sfx2>0);
%     figure
%     semilogy(ind_paint, sfx2(ind_paint),'.')
%     grid on
%     figure
%     semilogy(ind_neg, -sfx2(ind_neg),'.')
%     grid on

    sfx = (sfx2);
%     sfx = sfx2;

%     index_sfx = find(sfx<10^(-7));
%     sfx(index_sfx) = 10^(-7);
    

%     sizeNum = [50	100	100	200	200	400	400	800	800	800	1000    ];
%     sizeNum = 5*ones(1,17);
    j=1;
    
    % 时间段分割点
    timeline = t_st * 1e3 + [1  1.5   2   3   4   5   6   9       12   15     20   ];
    sizeNum = [50	80	100	150	200	300	300	500	    1000 1000	1000  ]; % 每个时间段对应的窗口宽度
    
    for j = 1:length(timeline)
        windowSize = sizeNum(j);
        b = (1/windowSize)*ones(1,windowSize); a = 1;
     
        filter_index =find(t*1e3 > timeline(j), 1);
        y = filter(b,a,sfx(:)');
        
        switch j
            case 1
                filter_signal = [filter_signal, pure_sec_field(i,1:150), y(151:filter_index-1)];
            case length(timeline)
                filter_signal = [filter_signal, y(lastidx:end)];
            otherwise
                filter_signal = [filter_signal, y(lastidx:filter_index-1)];
        end
        lastidx = filter_index;
    end
    
    
    %%
    %{ 
      % 滤波前后对比
        k=10;
        figure('Position',[511	255.666666666667	700	503.333333333333])
        loglog(t1, pure_sec_field(k,:))
        hold on
        loglog(t1, sfx)
    %     xlim([0,20])
        xlabel('time (ms)')
        ylabel('voltage (V)')
        hold on
        loglog(t1, filter_signal,'g')
        grid on
        legend('Input Data','30kHz Low pass filter','Average filter')
        xlim([0,20])
        title(['measurement point ',num2str(k)])
        set(gca,'FontName','Calibri','FontSize',16,'FontWeight','bold')

    %}
    filtered_sec_field = [filtered_sec_field; filter_signal];
    
end



%% 选择抽道时刻
delta_t = (log10(t_ed)-log10(t_st))/(nt-1);

time_log = zeros(1,nt);
ind_time = [];

for i=1:nt
    
    logt = log10(t_st)+(i-1)*delta_t;
    time_log(i) = 10^logt;
    
    tmp1 = t-time_log(i);
    
    tmp = find(tmp1>0, 1);

%     tmp2 = find(abs(tmp1) < 1e-6);
%     tmp = tmp2(1);
    
    if(i==1)
        if(isempty(tmp))
            ind_time = [ind_time; 1]; % 防止第一个抽道时刻，tmp为0的向量
        else
            ind_time = [ind_time; tmp];
        end
        
%         
    else
        ind_time = [ind_time; tmp];
    end
end

time_sample = t(ind_time);
signal_sample = filtered_sec_field(:,ind_time);

%%

%%
% {
% 单挑线画图保存，取点
k=1;
for i= 1:ns%
    %%
%     figure
%     loglog(t, sfx3(i,1:length(t)))
%     hold on
%     loglog(t, -sfx3(i,1:length(t)), '-r')
%     adjFig_1026(i)

    %%
%     figure
%     loglog(t, filtered_sec_field(i,1:length(t)))
%     hold on
%     loglog(t, -filtered_sec_field(i,1:length(t)), '-r')
%     adjFig_1026(i)
%     legend_str1 = ['抽道前测点',num2str(i)];
%     savefileName = fullfile(savefolder, [legend_str1,'.tif']);
%     saveas(gcf, savefileName)

    %%
    figure(Position=[10.333333333333	14.333333333333	1208.66666666667	618])
    loglog(time_sample*1000, signal_sample(i,:),'-o', 'LineWidth',1.0)
    hold on
    loglog(time_sample*1000, -signal_sample(i,:),'-ro', 'LineWidth',1.0)
    
    adjFig_1026(i); set(gca,'FontSize',24,'FontWeight','bold')
    for j = 1:nt
        text(time_sample(j)*1000, abs(signal_sample(i,j)), num2str(j), ...
        'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'bottom', 'FontSize', 10);
    end
    
    ylim([0.8e-5 0.5e-0])
    savefolder = '.\rawVolt';
    
    if ~exist(savefolder, 'dir')
        mkdir(savefolder);
        disp(['Folder "', savefolder, '" created.']);
    end
    legend_str1 = ['测点',num2str(i)];
    savefileName = fullfile(savefolder, [legend_str1,'.tif']);
    saveas(gcf, savefileName)
    close all
end
%}



%%
%{

figure('Position', [20	20	1260	628])
for i=1:nt
    semilogy(delta_pset.*(pset-min(pset)), 1e3 * signal_sample(:,i))
    hold on
    xlabel('Measurement Line / m')
    ylabel('Voltage / mV')
    legend_str2{i} = ['t = ',num2str(time_sample(i)*1000),' ms']; 
end
grid on
% legend(legend_str2,'NumColumns',4)
set(gca,'FontSize',14,'FontWeight','bold')
xlim([min(delta_pset.*(pset-min(pset))),max(delta_pset.*(pset-min(pset)))])

xlabel_pos = delta_pset.*(pset-min(pset));
for i = 1:ns
    text(xlabel_pos(i), 1e3 * signal_sample(i,2), ['',num2str(i)], ...
    'HorizontalAlignment', 'center', ...
    'VerticalAlignment', 'bottom', 'FontSize', 16);
end

% legend('')

savefileName = fullfile(savefolder, ['poumian.tif']);
saveas(gcf, savefileName)
close all


线性坐标系的剖面图

figure
for i=1:nt
    plot(delta_pset.*(pset-min(pset)), 1e3 * signal_sample(:,i))
    hold on
    xlabel('Measurement Line / m')
    ylabel('Voltage / mV')
    legend_str2{i} = ['t = ',num2str(time_sample(i)*1000),' ms'];
end
grid on
legend(legend_str2,'NumColumns',4)
set(gca,'FontSize',14,'FontWeight','bold')
xlim([min(delta_pset.*(pset-min(pset))),max(delta_pset.*(pset-min(pset)))])
%}
%%
save('signal_sample.mat','signal_sample')
step4_write_txt

function adjFig_1026(i)
    xlabel('time (s)');    ylabel('voltage (V)');    legend_str1 = ['测点',num2str(i)];
    grid on;legend(legend_str1,'NumColumns',4) ; set(gca,'FontSize',14,'FontWeight','bold')
end