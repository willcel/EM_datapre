%% 均值滤波处理程序 2ms


load('point1set.txt')
load('point4set.txt')
% load('data_avg.mat')
load('data_all_ns7.mat')

% 有个偏移量，参考step2_pre的结果
data_avg_all=data_offset(:,173:end);





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

numMeasure = size(data_avg_all, 1);
%%
for i=1:numMeasure
    filter_signal = [];
    t = time(index:end);    % 这里t的起始值就是抽道的起始时刻
    t = t(1:ind_neg_final);
    
    t1 = t*1000;
    
    % low pass filter 30kHz
    sfx1 = lowpass(data_avg_all(i,:),30000,fs,'ImpulseResponse','fir');
    sfx2 = sfx1(:,index:end);
    ind_neg = find(sfx2<0); 
    ind_paint = find(sfx2>0);


    sfx = abs(sfx2);

%     sizeNum = [50	100	100	200	200	400	400	800	800	800	1000    ];
    sizeNum = [50	100	200	300	300	600	600	1000	1000 1000	1000  ];

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



% {
% 单挑线画图保存，取点

for i= 1:numMeasure%
    figure(Position=[410.333333333333	174.333333333333	1808.66666666667	1118])
    loglog(time_sample*1000, signal_sample(i,:),'-*', 'LineWidth',1.0)
    
    hold on
    xlabel('time (ms)')
    ylabel('voltage (V)')
    legend_str1 = ['测点7 测量次数',num2str(i)];
    for j = 1:nt
        text(time_sample(j)*1000, signal_sample(i,j), num2str(j), ...
        'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'bottom', 'FontSize', 10);
    end
    grid on
    legend(legend_str1,'NumColumns',4)
    set(gca,'FontSize',24,'FontWeight','bold')
    savefolder = '.\rawVolt\测点7全测量次数';

    if ~exist(savefolder, 'dir')
        mkdir(savefolder);
        disp(['Folder "', savefolder, '" created.']);
    end

    savefileName = fullfile(savefolder, [legend_str1,'.tif']);
    saveas(gcf, savefileName)
    close all
end
%}


