clc
clear
close all
%%
load('signal_sample.mat');
% load('hz1obs_tunnel_2w.mat');
bac = signal_sample([1:29 41:64],:);
tun = signal_sample(30:40,:);

n1 = 6; %隧道测点
n2 = 44; %大地测点
nt = 56;

ch = 1:56;
s = signal_sample;

% y=data_avg_all(1,:);
% x = (1:length(y))./256e3;
% x1 = (1:length(y))./1;
% semilogy(x1,y)
% 777/256e3
%%
% sa = [];
% saa = [];
% for j = 1:n1
% 
%     t_st = 3.01e-3;
%     t_ed = 11e-3;
%     delta_t = (log10(t_ed)-log10(t_st))/(nt-1);
%     time_log = zeros(1,nt);
%     ind_time = [];
%     t = (t_st*(256e3):length(tun(j,:)))./(1e6);
%     for i=1:nt
%         logt = log10(t_st)+(i-1)*delta_t;
%         time_log(i) = 10^logt;
%         tmp1 = t-time_log(i);
%         tmp = find(abs(tmp1)<10^(-6));
%         if(i==1)
%             ind_time = [ind_time; 1]; % 防止第一个抽道时刻，tmp为0的向量
%         else
%             ind_time = [ind_time; tmp(1)];
%         end
%     end
%     s11 = tun(j,round(ind_time));
%     saa = [saa;s11];
% end
% 
% sb  = [];
% sbb  = [];
% for j = 1:n2
% %     nbac=awgn(bac,20*log10(1/0.05),'measured');
% 
%     noise = randn(size(bac)) * epsilon * rms(bac);
%     nbac = bac + noise;
% 
%     t_st = 2.04e-3;
%     t_ed = 20e-3;
% 
%     delta_t = (log10(t_ed)-log10(t_st))/(nt-1);
%     time_log = zeros(1,nt);
%     ind_time = [];
%     t = (t_st*1e6:length(tun))./(1e6);
%     for i=1:nt
%         logt = log10(t_st)+(i-1)*delta_t;
%         time_log(i) = 10^logt;
%         tmp1 = t-time_log(i);
%         tmp = find(abs(tmp1)<10^(-6));
%         if(i==1)
%             ind_time = [ind_time; 1]; % 防止第一个抽道时刻，tmp为0的向量
%         else
%             ind_time = [ind_time; tmp(1)];
%         end
%     end
%     s1 = nbac(round(ind_time));
%     s11 = bac(round(ind_time));
%     sb = [sb;s1];
%     sbb = [sbb;s11];
% end    

% s = [sb(1:n2/2,:); sa; sb(n2/2+1:n2,:)];
% ss = [sbb(1:n2/2,:); saa; sbb(n2/2+1:n2,:)];

% h1 = semilogy(abs(s(:,[3  27 45 48 51])),'r-','linewidth',1.2);
% hold on
% h2=semilogy(abs(ss(:,[ 3  27  45 48 51])),'b-','linewidth',1.1);
% legend([h1(1) h2(1)],{'withnoise','noisefree'})
% xlabel('survey point')
% ylabel('dB/dt(nT/s)')
% set(gca,'fontsize',14)

%%
figure(1)
semilogy(s(:,:),'-','linewidth',1);
% hold on
% plot(s(:,ch(2)),'g-','linewidth',1);
% hold on
% plot(s(:,ch(3)),'b-','linewidth',1);
% hold on
% plot(ss(:,ch(1)),'r-*','linewidth',1.8);
% hold on
% plot(ss(:,ch(2)),'g-*','linewidth',1.8);
% hold on
% plot(ss(:,ch(3)),'b-*','linewidth',1.8);
% hold on
% legend('chan45-noise','chan48-noise','chan51-noise','chan45-noisefree','chan48-noisefree','chan51-noisefree')
title('原始')
xlabel('survey point')
ylabel('dB/dt(nT/s)')
set(gca,'fontsize',14)
set(gca,'linewidth',1.2)
grid on
%% step1 
% mx1 = mean(s([1:22 29:end],:),1);
% mx2 = mean(s(23:28,:),1);
% mx = [mx1(ones(22,1),:); mx2(ones(6,1),:); mx1(ones(22,1),:) ] ;
mx = mean(s,1);
X = s - mx;
X = X';
% [m,n] = size(X);
% Cn = zeros(m,m);
% for k = 1:m
%     for q = 1:m
%         for j = 1:n-1           
%             Zk(j) = X(k,j+1)-X(k,j);
%             Zq(j) = X(q,j+1)-X(q,j);    
%         end
%         Zk_ = mean(Zk);
%         Zq_ = mean(Zq);      
%        Cn(k,q) =  0.5*sum((Zk-Zk_).*(Zq-Zq_))/(n-1);   
%     end
% end

DET = [X(:,1:30)-X(:,2:31)  X(:,40:63)-X(:,41:64)];
% noi = ss-s;
% Cn = cov(noi);
Cn = cov(DET')/2;
%% step2
[U, DN, Ut] = svd(Cn);
P = Ut*inv(DN.^0.5);
C = cov(X'); 
Cd = P'*C*P;
[V, DD, Vt] = svd(Cd);
R = P*Vt;
ph = R'*X;

l = 1;
up = sum(sum(DD(1:l,1:l)));
dn = sum(sum(DD));
up/dn

B = diag([ones(l,1) ;zeros(nt-l,1)]);
Xhat = inv(R')*B*ph  + mx' ;      
Xhat = Xhat';
%%
% save Xhat Xhat
%%
% X1 = X
figure(2)
semilogy(Xhat(:,:),'-','linewidth',1.1);
title('MNF')
% hold on
% plot(Xhat(:,ch(2)),'r-','linewidth',1.1);
% hold on
% plot(Xhat(:,ch(3)),'r-','linewidth',1.1);
% hold on
% plot(ss(:,ch(1)),'b-','linewidth',1.1);
% hold on
% plot(ss(:,ch(2)),'b-','linewidth',1.1);
% hold on
% plot(ss(:,ch(3)),'b-','linewidth',1.1);
% legend('chan45-denoise','chan48-denoise','chan51-denoise')
xlabel('survey point')
ylabel('dB/dt(nT/s)')
set(gca,'fontsize',14)
set(gca,'linewidth',1.2)
grid on
1;




