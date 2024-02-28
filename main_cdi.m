%%
clc
clear
% close all
dbstop if error
%% ---------------------- 参数设置 ----------------------------------
addpath("D:\willcel\subfunc_mat")

%% 数据预处理的参数
% 0.5+(1:11)*0.5
pset = [1:163];  % 测点的坐标，文件夹的名称
delta_pset = 1;            % 测点之间的距离 （m）

ns = length(pset);                  % 测点的个数

% no_of_measurement = 100;        % 每个测点的测量次数

fs = 256e3;  % 1.25*10^6;

%% ------------------- 数据预处理 ----------------------------------



%% 数据采样
nt = 56;                        % 抽道时间
t_st = 2.04e-3; %            % 起始时间        
t_ed = 20e-3;       % 结束时间 

%% 反演参数
total_depth = 40;           % 最大深度 m
nolayer = 5;

%% 发射参数
factor = 1/0.22/10;                        % 回波的放大倍数
factor_current = 1/0.22*200;              % 电流的放大倍数

hr = 0.63;   % 接收线圈的高度

rt = 0.5;                 % 发射线圈的半径 m
nturn = 3;              % 线圈的匝数

rr = 0.25;               % 接收线圈的半径 m
nturn1 = 80;          % 接收线圈的匝数

xr = 0.59;    % 中心距

%%
step1    % 数据去直流偏置, 叠加求平均
%%
% 要设 t_ed
step3                          % 滤波
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

% -------------------- 视电阻率成像 -------------------------------
%%

cdi_2

%%

get_priori
% get_priori_depDivide

%%
cd ..
% git clone git@github.com:willcel/EM_singleBP.git
cd .\EM_singleBP
