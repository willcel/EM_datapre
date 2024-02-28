%%
clc
clear
% close all
dbstop if error
%% ---------------------- �������� ----------------------------------
addpath("D:\willcel\subfunc_mat")

%% ����Ԥ����Ĳ���
% 0.5+(1:11)*0.5
pset = [1:163];  % �������꣬�ļ��е�����
delta_pset = 1;            % ���֮��ľ��� ��m��

ns = length(pset);                  % ���ĸ���

% no_of_measurement = 100;        % ÿ�����Ĳ�������

fs = 256e3;  % 1.25*10^6;

%% ------------------- ����Ԥ���� ----------------------------------



%% ���ݲ���
nt = 56;                        % ���ʱ��
t_st = 2.04e-3; %            % ��ʼʱ��        
t_ed = 20e-3;       % ����ʱ�� 

%% ���ݲ���
total_depth = 40;           % ������ m
nolayer = 5;

%% �������
factor = 1/0.22/10;                        % �ز��ķŴ���
factor_current = 1/0.22*200;              % �����ķŴ���

hr = 0.63;   % ������Ȧ�ĸ߶�

rt = 0.5;                 % ������Ȧ�İ뾶 m
nturn = 3;              % ��Ȧ������

rr = 0.25;               % ������Ȧ�İ뾶 m
nturn1 = 80;          % ������Ȧ������

xr = 0.59;    % ���ľ�

%%
step1    % ����ȥֱ��ƫ��, ������ƽ��
%%
% Ҫ�� t_ed
step3                          % �˲�
% step4_write_txt                  % ��ԭʼ����д���ļ�

% { 
% -----------
path_cdi = ['.\CDI_code\'];   % �ӵ����ʳ���ĳ��������ļ���
para = [nt; ns; nturn; nturn1; t_st; t_ed; xr; hr; rt; rr];
save('parameter_in.txt','para','-ascii')
path_cdi = ['.\CDI_code\'];   % �ӵ����ʳ���ĳ��������ļ���
% �ƶ��ļ�
copyfile('point1set.txt', path_cdi)
copyfile('point2set.txt', path_cdi)
copyfile('point3set.txt', path_cdi)
copyfile('point4set.txt', path_cdi)
copyfile('parameter_in.txt', path_cdi)
%}

% -------------------- �ӵ����ʳ��� -------------------------------
%%

cdi_2

%%

get_priori
% get_priori_depDivide

%%
cd ..
% git clone git@github.com:willcel/EM_singleBP.git
cd .\EM_singleBP
