%%
clc
clear
close all
dbstop if error
%% ---------------------- �������� ----------------------------------

%% ����Ԥ����Ĳ���
% 0.5+(1:11)*0.5
pset = [1:25];  % �������꣬�ļ��е�����
delta_pset = 1;            % ���֮��ľ��� ��m��

ns = length(pset);                  % ���ĸ���

factor = 1/0.65/10;                             % �ز��ķŴ���
factor_current = 1/0.65*200;              % �����ķŴ���

% no_of_measurement = 100;        % ÿ�����Ĳ�������

fs = 1.25*10^6;

%% ���ݲ���
nt = 56;                        % ���ʱ��
t_st = 2.2e-3; % 2.02e-3; %           % ��ʼʱ��        
t_ed = 20e-3;       % ����ʱ�� 


%% ���ݲ���
total_depth = 30;           % ������ m
nolayer = 5;

%% �������
hr = 0.01;   % ������Ȧ�ĸ߶�

rt = 0.5;                 % ������Ȧ�İ뾶 m
nturn = 3;              % ��Ȧ������

rr = 0.25;               % ������Ȧ�İ뾶 m
nturn1 = 80;          % ������ Ȧ������

xr = 0.58;    % ���ľ�

%% ------------------- ����Ԥ���� ----------------------------------
% step1    % ����ȥֱ��ƫ��, ������ƽ��
% step2_pre
% step2     % �洢������ת�۵�

% Ҫ�� t_ed
% step3                          % �˲�
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
%{
% -------------------- �ӵ����ʳ��� -------------------------------
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
