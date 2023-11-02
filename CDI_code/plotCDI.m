% clear
rho_a = load("rho1.dat");
hz1_a = load('log1.dat');
%%
% close all
% volt = load("hz1raw.dat");

tst = 3.1;
ted = 20;
nt = 56;
% time(1) = 2.04;
delta = (log10(ted)-log10(tst)) / (nt-1);
time = 10.^(log10(tst) + (0:nt-1) * delta);

figure; figidx=1;
% get(gcf, 'Position')
set(gcf, "Position", [1.0915    1.0610    0.6140    0.7215]*1e3)

% for i = 150:50:350 %400:10:460%1:100:500
for i = 400:10:460 
% i = 400;
rho_cur = rho_a(i);
volt2 = hz1_a(i, :);


subplot(5,1,figidx) ; figidx = figidx+1;
% semilogy(volt, '-ro')
loglog(time, (volt2), '-b')
hold on
loglog(time, -(volt2), '-r')
% loglog((volt2), '-bo')
grid on
legend('正值','负值')
title( ['电阻率为', num2str(rho_cur)])
end
% legend('原响应曲线', '极化模型曲线')
%%
figure
% semilogy
plot(time, (volt2), '-b')
grid on

%% 画一下查表的图

figure
for timeNt = 1:10:510
volt_perRho = hz1_a(:, timeNt);


loglog(rho_a, volt_perRho, '-')
hold on
end
grid on
title('正值部分')
xlabel('rho (欧米)')
ylabel('V');set(gca,'FontSize',12,'FontWeight','bold')
%%

figure
for timeNt = 1:10:510
volt_perRho = hz1_a(:, timeNt);


loglog(rho_a, -volt_perRho, '-')
hold on
end
xlabel('rho (欧米)')
ylabel('V')
grid on
title('负值部分');set(gca,'FontSize',12,'FontWeight','bold')

%% 画一下不同eta下的转折点
figure
rho_a = load("rho1.dat");
rho_seleIndx = 250;
rho_sele = rho_a(rho_seleIndx);
for eta = 0.3:0.3:0.9
    filename = sprintf("log%0.1f.dat", eta);
    volt = load(filename);
    volt2 = volt(rho_seleIndx, :);
    loglog(time, (volt2), '-bo')
    hold on
    loglog(time, -(volt2), '-r*')
    grid on
    
    1;
end
legend('正值','负值')
set(gca,'FontSize',14,'FontWeight','bold')
xlabel('time(ms)')
ylabel('V')
svfig('不同eta')