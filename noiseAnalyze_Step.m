load("data_test.mat")
% load("data_1.mat")
data_test;
%
figure
% time = 1:length(data_avg_all)
% subplot(2,1,1)
time = t1*1e-3;
loglog(time*1000, data_test(:),'LineWidth',2.0)
hold on
loglog(time*1000, -data_test(:),'LineWidth',2.0)
% loglog(time*1000, data_avg_all(selePoint,:).*factor,'LineWidth',2.0)
% xlim([2.2 4])
set(gca,'FontName','Calibri','FontSize',12,'FontWeight','bold')
xlabel('time(ms)')
ylabel('voltage (V)')
grid on
% xlim([0 5])
% title(['cedian',num2str(k)])



%%
s1 = data_test;
t1 = time;

% sfx1 = lowpass(s1,30000,fs,'ImpulseResponse','fir','Steepness',0.95);
% s1 = sfx1;


stpoint = 30e-3*fs; 
edpoint = 250e-3*fs; 

% stpoint = 16e-3*fs; 
% edpoint = 40e-3*fs; 

indarr = (stpoint:edpoint);

s1 = s1(indarr);
t1 = t1(indarr);

%%
figure
drawFFT(s1, fs)
% xlim([10 100])
% ylim([0 1.8e-6])
hold on

%%
% svfig('测点3频谱比较')
% legend("测点3频谱", "测点2频谱")
% legend("15ms前的频谱", "15ms后的频谱")