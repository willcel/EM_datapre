%%
clc
clear
close all
dbstop if error
%% ---------------------- 参数设置 ----------------------------------

%% 数据预处理的参数
% 0.5+(1:11)*0.5
pset = [1:25];  % 测点的坐标，文件夹的名称
delta_pset = 1;            % 测点之间的距离 （m）

ns = length(pset);                  % 测点的个数

factor = 1/0.65/10;                             % 回波的放大倍数
factor_current = 1/0.65*200;              % 电流的放大倍数

% no_of_measurement = 100;        % 每个测点的测量次数

fs = 1.25*10^6;

%% 数据采样
nt = 56;                        % 抽道时间
t_st = 2.2e-3; % 2.02e-3; %           % 起始时间        
t_ed = 20e-3;       % 结束时间 


%% 反演参数
total_depth = 30;           % 最大深度 m
nolayer = 5;

%% 发射参数
hr = 0.01;   % 接收线圈的高度

rt = 0.5;                 % 发射线圈的半径 m
nturn = 3;              % 线圈的匝数

rr = 0.25;               % 接收线圈的半径 m
nturn1 = 80;          % 接收线 圈的匝数

xr = 0.58;    % 中心距

%% ------------------- 数据预处理 ----------------------------------
% step1    % 数据去直流偏置, 叠加求平均
% step2_pre
% step2     % 存储电流的转折点

% 要设 t_ed
% step3                          % 滤波
% step4_write_txt                  % 将原始数据写入文件

% { 
% -----------
path_cdi = ['.\CDI_code\'];   % 视电阻率成像的程序所在文件夹
para = [nt; ns; nturn; nturn1; t_st; t_ed; xr; hr; rt; rr];
save('parameter_in.txt','para','-ascii')
path_cdi = ['.\CDI_code\'];   % 视电阻率成像的程序所在文件夹
% 移动文件
copyfile('point1set.txt', path_cdi)
copyfile('point2set.txt', path_cdi)
copyfile('point3set.txt', path_cdi)
copyfile('point4set.txt', path_cdi)
copyfile('parameter_in.txt', path_cdi)
%}
%{
% -------------------- 视电阻率成像 -------------------------------
delete([path_cdi,'flag.dat'])
winopen([path_cdi,'CDI_code.exe'])
% 
while 1
    pause(1);
    
    if(exist([path_cdi,'flag.dat'],'file'))
        a = load([path_cdi,'flag.dat']);
        if (size(a,1) == 1)
            break
        end
    end
end
%}
