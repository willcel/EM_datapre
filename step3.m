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

end_point = t_st;  % mean(point4set(1:2:end));
disp(['end_point = ', num2str(end_point*1000)])

[~, start_index] =min(abs(time-start_point));  %
[~, end_index] =min(abs(time-end_point));    %

time = time-start_point;  % 时间从0开始

index = end_index+1;             % 必须大于end point % time(index) 彰显起始时刻


% pure_sec_field = (data_avg_all(:,index:end));
pure_sec_field = abs(data_avg_all(:,index:end));

ind_neg_final = size(data_avg_all, 2)-index+1;

filtered_sec_field = [];

%%
for i=1:ns
    
    t = time(index:end);    % 这里t的起始值就是抽道的起始时刻
    t = t(1:ind_neg_final);
    
    t1 = t*1000;
    
    % low pass filter 30kHz
    sfx1 = lowpass(data_avg_all(i,:),10000,fs,'ImpulseResponse','fir','Steepness',0.95);
%     sfx1 = lowpass(data_avg_all(i,:),30000,fs,'ImpulseResponse','fir');
    sfx2 = sfx1(:,index:end);

    ind_neg = find(sfx2<0); 
    ind_paint = find(sfx2>0);
    sfx = abs(sfx2);
    
    j=1;
    
    % 时间段分割点
           timeline = [0.5  1   2   3   4   5   6   9   12   15    20   ] + t_st * 1e3;
    switch i
        otherwise
            sizeNum = [50	50	100	150	200	300	400	500	800 1000	1000  ]; % 每个时间段对应的窗口宽度
    end
    
    %% 滤波50Hz
    % {
    
    idx_arr = 4.5e-3*fs : length(sfx2);
    sig_tar = sfx2(idx_arr);
% 
%{
    beipin_bank = [];
    for iter = 1:1
        tic
        f0 = 50;
        [peak1, peak2, phi0] = beipin(sig_tar,f0,fs);
        beipin_bank(iter,:) = [f0, peak1, peak2, phi0];
        toc
    end
    sig_tar2 = sig_tar; 

    N = length(sig_tar2);
    N2 = idx_arr(1) - 2 ;
    t2 = (-N2:N) / fs;
    for iter = 1:1
        f0 = 50;
        peak2 = beipin_bank(iter,3);
        phi0 = beipin_bank(iter,4);
        cfs = peak1*sin(2*pi*f0*t2 + phi0);
%         sig_tar2 = sig_tar2 - cfs;
        sfx2 = sfx2 - cfs;
        sfx3 = abs(sfx2);
    end
%}
    sig_tar2 = notchfilt(sig_tar, 50, fs);

%     figure
%     drawFFT(sig_tar, fs)
%     xlim([0 1e3])
%     hold on; 
%     drawFFT(sfx2(idx_arr), fs)
%     legend('去噪前','去噪后')



    filter_signal2 = meanFilt(t, sfx3, timeline, sizeNum, pure_sec_field(i,:));
    %}
    %%

    filter_signal = meanFilt(t, sfx, timeline, sizeNum, pure_sec_field(i,:));
    

    %{ 
      % 滤波前后对比
        %%
%             if(mod(i,20)==1)
            if(i==1)
                
%             tst = 13; ted = 23.5;
%             tst = 11; ted = 23.5;
%             filterNew = cal_expfit(filter_signal, t1, tst, ted);
        
            figure('Position',[511	255.666666666667	700	503.333333333333])
            loglog(t1, pure_sec_field(i,:))
            hold on
%             loglog(t1, sfx)
%             loglog(t1, sfx3, "LineWidth", 1)
            loglog(t1, filter_signal,'g', "LineWidth", 1.5)
            loglog(t1, filter_signal2,'r', "LineWidth", 1.5)
%             loglog(t1, filterNew, 'm', "LineWidth", 1.5)
            
            
            legend('Input Data','低通+均值','低通+陷波+均值')
            xlim([0,80])
            ylim([1e-10 1])
            title(['measurement point ',num2str(i)])
            set(gca,'FontSize',16,'FontWeight','bold')
            xlabel('time (ms)');  ylabel('voltage (V)') ; hold on; grid on

            svfig(num2str(i), '.\rawVolt\step3滤波对比')

            end
    
    %}
    close all

    filtered_sec_field = [filtered_sec_field; filter_signal2];
%     filtered_sec_field = [filtered_sec_field; filter_signal2];
    
end



%% 选择抽道时刻
[time_sample, signal_sample] = chouDao(t, filtered_sec_field, t_st, t_ed, nt);

% load Xhat
% signal_sample = Xhat;
%%
% figure
% loglog(t, filtered_sec_field(10,1:length(t)))
%%
savefolder = '.\rawVolt';
% {

%% 单挑线画图保存，取点
% signal_sample(160,:) = (signal_sample(159,:) +signal_sample(161,:) )/2;
% signal_sample(8,:) = cal_expfit(signal_sample(8,:), time_sample, time_sample(38), time_sample(42));

% for i= 1:ns%
%     figure(Position=[10.333333333333	14.333333333333	908.66666666667	518])
%     loglog(time_sample*1000, signal_sample(i,:),'-*', 'LineWidth',1.0)
%     
%     legend_str = ['测点',num2str(i)];
% 
%     hold on; grid on;
%     xlabel('time (ms)');  ylabel('voltage (V)')
%     set(gca,'FontSize',24,'FontWeight','bold')
%     legend(legend_str,'NumColumns',4)
%     for j = 1:nt
%         text(time_sample(j)*1000, signal_sample(i,j), num2str(j), ...
%         'HorizontalAlignment', 'center', ...
%         'VerticalAlignment', 'bottom', 'FontSize', 10);
%     end
%     ylim([1e-6 1e-1])
%     svfig(legend_str, savefolder)
%     close all
%     
% end



%% 把所有抽道的线画在同一张图中
% k=1;
% figure(Position=[10.333333333333	14.333333333333	908.66666666667	518])
% 
% for i= 1:ns%
%     
%     loglog(time_sample*1000, signal_sample(i,:),'-*', 'LineWidth',1.0)
%     
%     legend_str3{k} = ['测点',num2str(i)]; k= k+1;
% 
%     hold on; grid on;
%     xlabel('time (ms)');  ylabel('voltage (V)')
%     set(gca,'FontSize',24,'FontWeight','bold')
% 
% %     ylim([1e-6 1e-1])
% end
% legend(legend_str3,'NumColumns',4)
% svfig('抽道图', savefolder)




%%
% 剖面图

figure('Position', [20	20	1260	628])
for i=1:nt
    semilogy(delta_pset.*(pset-min(pset)), 1e3 * signal_sample(:,i))
    hold on
    xlabel('Measurement Line / m')
    ylabel('Voltage / mV')
    legend_str2{i} = ['t = ',num2str(time_sample(i)*1000),' ms']; 
end
grid on

set(gca,'FontSize',14,'FontWeight','bold')
xlim([min(delta_pset.*(pset-min(pset))),max(delta_pset.*(pset-min(pset)))])

xlabel_pos = delta_pset.*(pset-min(pset));
for i = 1:ns
    text(xlabel_pos(i), 1e3 * signal_sample(i,2), ['',num2str(i)], ...
    'HorizontalAlignment', 'center', ...
    'VerticalAlignment', 'bottom', 'FontSize', 16);
end

% legend('')

svfig('剖面图', savefolder)
legend(legend_str2,'NumColumns',4)
svfig('剖面图带legend', savefolder)
% close all
%}
%{

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

