% 2010.10.29
% 
% ���ܣ���ʮ��������ת��Ϊʮ���������з��ţ�
% ���ô˺�����  hex2decWithSign({'FFFE', '0002'}, 4);

function decval=hex2decWithSign(hexval, lengt)

test = 0;
if (test ==1)
    clc;close all;clear all;
    hexval={ '01451AC2', '005ACA6D', 'FF8237B5', '00855EFC', '00726593', 'FFA3CB58'};
    lengt = 8;    % λ��
end

decval = hex2dec(hexval);
sign = bitget(decval, 4*lengt);             % ���λΪ����λ��ȡ�������λ
negative_numbers = (sign == 1);              % ȡΪ�߼����桱�򡰼١�
pvivn_decval = decval(negative_numbers);     % ȡ��decval��Ϊ������

decval(negative_numbers) = pvivn_decval - bitshift(1, 4*lengt);  % ת��Ϊ�з�����  -2^32
