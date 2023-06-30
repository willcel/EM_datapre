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


pure_sec_field = abs(data_avg_all(:,index:end));

ind_neg_final = size(data_avg_all, 2)-index+1;

filtered_sec_field = [];

%{ 
没有滤波的二次场信号


t = time(index:end);    % 这里t的起始值就是抽道的起始时刻
t = t(1:ind_neg_final);
t1 = t*1000;
figure('Position',[1000	608.333333333333	931	729.666666666667])
for i=1:ns
    loglog(t1, pure_sec_field(i,:))
    hold on
    legend_str1{i} = ['测点',num2str(i)];
end
xlim([0,20])
xlabel('time (ms)')
ylabel('Voltage (V)')
legend(legend_str1,'NumColumns',4)
grid on
set(gca,'FontSize',16,'FontWeight','bold')

%}

%%
for i=1:ns
    filter_signal = [];
    t = time(index:end);    % 这里t的起始值就是抽道的起始时刻
    t = t(1:ind_neg_final);
    
    t1 = t*1000;
    
    % low pass filter 30kHz
    sfx1 = lowpass(data_avg_all(i,:),30000,fs,'ImpulseResponse','fir');
    sfx2 = sfx1(:,index:end);
    ind_neg = find(sfx2<0); 
    ind_paint = find(sfx2>0);
%     figure
%     semilogy(ind_paint, sfx2(ind_paint),'.')
%     grid on
%     figure
%     semilogy(ind_neg, -sfx2(ind_neg),'.')
%     grid on

    sfx = abs(sfx2);
%     sfx = sfx2;

%     index_sfx = find(sfx<10^(-7));
%     sfx(index_sfx) = 10^(-7);
    

%     sizeNum = [50	100	100	200	200	400	400	800	800	800	1000    ];
    sizeNum = [50	100	200	300	300	600	600	1000	1000 1000	1000  ];
    if(ismember(i,[2 6 ]))
        sizeNum = [50	100	200	300	1000	1000	1000	1000	1000 1000	1000 ];
    elseif(ismember(i,[8 9]))
        sizeNum = [50	100	200	300	1000	2000	2000	2000	2000 1000	1000 ];
    elseif(ismember(i,[11]))
        sizeNum = [50	100	200	300	1000	1500	1000	1000	1000 1000	1000 ];
    elseif(ismember(i,[12]))
        sizeNum = [50	100	200	300	1000	2000	3000	1000	1000 1000	1000 ];
    elseif(ismember(i,[20]))
        sizeNum = [50	100	800	1000	1000	2000	3000	1000	1000 1000	1000 ];
    end
%     sizeNum = 5*ones(1,17);
    j=1;
    %% 2.3ms
    windowSize = sizeNum(j); j=j+1;
    b = (1/windowSize)*ones(1,windowSize);
    a = 1;
     
    filter_index =find(abs(t-2.3*10^(-3))<5*10^(-7));
    y = filter(b,a,sfx(:)');
    
    filter_signal = [filter_signal, pure_sec_field(i,1:150), y(151:filter_index-1)];
    
    %% 3ms
    windowSize = sizeNum(j); j=j+1;
    
    % if(i==4 || i==9 || i==16)
    %      windowSize = sizeNum(j); j=j+1;;
    % end
    b = (1/windowSize)*ones(1,windowSize);
    a = 1;
    
    filter_index1 =find(abs(t-3*10^(-3))<5*10^(-7));
    y = filter(b,a,sfx(:)');
    
    filter_signal = [filter_signal, y(filter_index:filter_index1-1)];
    
    %% 3.5ms
    
    windowSize = sizeNum(j); j=j+1;
    
    b = (1/windowSize)*ones(1,windowSize);
    a = 1;
    
    filter_index2 =find(abs(t-3.5*10^(-3))<5*10^(-7));
    y = filter(b,a,sfx(:)');
    
    filter_signal = [filter_signal, y(filter_index1:filter_index2-1)];
    
    
    %% 4ms
    windowSize = sizeNum(j); j=j+1;
    
    b = (1/windowSize)*ones(1,windowSize);
    a = 1;
    
    filter_index3 =find(abs(t-4*10^(-3))<5*10^(-7));
    y = filter(b,a,sfx(:)');
    
    filter_signal = [filter_signal, y(filter_index2:filter_index3-1)];
    
    %% 4.5ms
    windowSize = sizeNum(j); j=j+1;
    
    b = (1/windowSize)*ones(1,windowSize);
    a = 1;
    
    filter_index4 =find(abs(t-4.5*10^(-3))<5*10^(-7));
    y = filter(b,a,sfx(:)');
    
    filter_signal = [filter_signal, y(filter_index3:filter_index4-1)];
    
    %% 5ms -para 6
    windowSize = sizeNum(j); j=j+1;
    
    
    b = (1/windowSize)*ones(1,windowSize);
    a = 1;
    
    filter_index5 =find(abs(t-5*10^(-3))<5*10^(-7));
    y = filter(b,a,sfx(:)');
    
    filter_signal = [filter_signal, y(filter_index4:filter_index5-1)];
    
    %% 6.5ms
    
    windowSize = sizeNum(j); j=j+1;
    
    b = (1/windowSize)*ones(1,windowSize);
    a = 1;
    
    filter_index6 =find(abs(t-6.5*10^(-3))<5*10^(-7));
    y = filter(b,a,sfx(:)');
    
    filter_signal = [filter_signal, y(filter_index5:filter_index6-1)];
    
    
    %% 8.5ms
    windowSize = sizeNum(j); j=j+1;
    
    b = (1/windowSize)*ones(1,windowSize);
    a = 1;
    
    filter_index7 =find(abs(t-8.5*10^(-3))<5*10^(-7));
    y = filter(b,a,sfx(:)');
    
    filter_signal = [filter_signal, y(filter_index6:filter_index7-1)];
    
    %% 12ms
    windowSize = sizeNum(j); j=j+1;
    
    b = (1/windowSize)*ones(1,windowSize);
    a = 1;
    
    filter_index8 =find(abs(t-12*10^(-3))<5*10^(-7));
    y = filter(b,a,sfx(:)');
    
    filter_signal = [filter_signal, y(filter_index7:filter_index8-1)];
    
    %% 16ms
    windowSize = sizeNum(j); j=j+1;

%     if(i==4 || i==3 || i==14 || i==11)
% %          windowSize = sizeNum(j); j=j+1;;
%     end

    b = (1/windowSize)*ones(1,windowSize);
    a = 1;
    
    filter_index9 =find(abs(t-16*10^(-3))<5*10^(-7));
    y = filter(b,a,sfx(:)');
    
    filter_signal = [filter_signal, y(filter_index8:filter_index9-1)];
    
    %% 20ms
    
    windowSize = sizeNum(j); j=j+1;
    
    %  if(i==4 || i==3 || i==14 || i==11)
    %      windowSize = sizeNum(j); j=j+1;;
    %  end

    
    b = (1/windowSize)*ones(1,windowSize);
    a = 1;
    
    y = filter(b,a,sfx(:)');
    
    filter_signal = [filter_signal, y(filter_index9:end)];
    
    
    %%
    %{ 
      % 滤波前后对比

        figure('Position',[511	255.666666666667	700	503.333333333333])
        loglog(t1, pure_sec_field(i,:))
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
        title(['measurement point ',num2str(i)])
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
    
    tmp = find(abs(tmp1)<10^(-6));
    
    if(i==1)
        ind_time = [ind_time; 1]; % 防止第一个抽道时刻，tmp为0的向量
    else
        ind_time = [ind_time; tmp(1)];
    end
end

time_sample = t(ind_time);
signal_sample = filtered_sec_field(:,ind_time);

% signal_sample(5,:) = signal_sample(2,:);
% signal_sample(6,:) = signal_sample(1,:);
% signal_sample(21,:) = signal_sample(2,:);
% signal_sample(22,:) = signal_sample(2,:);
% signal_sample(23,:) = signal_sample(2,:);
%%
%{
figure(Position=[410.333333333333	174.333333333333	1808.66666666667	1118])
k=1;
for i= 1:ns%
    
    loglog(time_sample*1000, signal_sample(i,:))%,'-*')
    
    hold on
    xlabel('time (ms)')
    ylabel('voltage (V)')
    legend_str1{k} = ['测点',num2str(i)];k=k+1;
%     for j = 1:nt
%         text(time_sample(j)*1000, signal_sample(i,j), num2str(j), ...
%         'HorizontalAlignment', 'center', ...
%         'VerticalAlignment', 'bottom', 'FontSize', 6);
%     end
end
grid on
% legend('1','2','3')
% legend('5','6','7','8','9')
legend(legend_str1,'NumColumns',4)
set(gca,'FontSize',24,'FontWeight','bold')
%}

% signal_sample(17,[22 23 24]) = [1.4 1.3 1.3]*1E-5;
signal_sample(20,[6 7]) = [4.5E-5 2.8E-5];

% {
% 单挑线画图保存，取点
k=1;
for i= 8%1:ns%
    figure(Position=[410.333333333333	174.333333333333	1808.66666666667	1118])
    loglog(time_sample*1000, signal_sample(i,:),'-*', 'LineWidth',1.0)
    
    hold on
    xlabel('time (ms)')
    ylabel('voltage (V)')
    legend_str1 = ['测点',num2str(i)];
    for j = 1:nt
        text(time_sample(j)*1000, signal_sample(i,j), num2str(j), ...
        'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'bottom', 'FontSize', 10);
    end
    grid on
    legend(legend_str1,'NumColumns',4)
    set(gca,'FontSize',24,'FontWeight','bold')
    savefolder = '.\rawVolt';

    if ~exist(savefolder, 'dir')
        mkdir(savefolder);
        disp(['Folder "', savefolder, '" created.']);
    end

    savefileName = fullfile(savefolder, [legend_str1,'.tif']);
    saveas(gcf, savefileName)
%     close all
end
%}



%%
%{

figure('Position', [200	200	2060	1028])
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

legend('')

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