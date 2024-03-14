clc
clear
close all
%%

nt = 56;
load("signal_sample.mat")

Lprior = 2;
[Xhat, energyRatio] = MNF(signal_sample, Lprior);
%%
% save Xhat Xhat
%%
% X1 = X
figure(2)
semilogy(Xhat(:,:),'-','linewidth',1.1);
title(sprintf('前%d个成分 PCA, 能量占比%.2f%%', Lprior, energyRatio*100))
adjustFig
xlabel('survey point')
% ylabel('dB/dt(nT/s)')
set(gca,'fontsize',14)
set(gca,'linewidth',1.2)
grid on
1;




