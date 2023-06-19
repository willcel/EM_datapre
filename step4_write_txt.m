%% 将观测数据写入txt文件中

load('signal_sample.mat')
signal = signal_sample;



delta_t = (log10(t_ed)-log10(t_st))/(nt-1);

time_log = zeros(1,nt);

for i=1:nt
    logt = log10(t_st)+(i-1)*delta_t;
    time_log(i) = 10^logt;
    
end



%%
% figure
% id = 0;
% for i=1:ns
%     
%     loglog(time_log.*1000, signal(i,:))
%     hold on
%     xlabel('time (ms)')
%     ylabel('voltage (mV)')
%     id = id+1;
%     legend_str1{id} = ['测点',num2str(i)];
% end
% legend(legend_str1,'NumColumns',4)
% set(gca,'FontSize',12,'FontWeight','bold')
% grid on
% 
% measureline = 1:ns;
% 
% figure
% for i=1:nt
%     semilogy(measureline, signal(:,i))
%     hold on
%     xlabel('measurement line(m)')
%     ylabel('voltage (mV)')
%     legend_str2{i} = ['t = ',num2str(time_log(i)),' s'];
% end
% grid on
% set(gca,'FontSize',12,'FontWeight','bold')
% legend(legend_str2,'NumColumns',4)

%%
vobs = zeros(size(signal,1)*size(signal,2), 1);

for uu=1:ns
    for i=1:nt
        vobs(i+(uu-1)*nt) = signal(uu, i);
    end
end

namestep4 = sprintf('vobs_%.0fms.txt',t_ed*1e3);
save(namestep4,'vobs','-ascii')

