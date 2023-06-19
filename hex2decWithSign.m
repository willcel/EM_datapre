% 2010.10.29
% 
% 功能：将十六进制数转换为十进制数（有符号）
% 调用此函数用  hex2decWithSign({'FFFE', '0002'}, 4);

function decval=hex2decWithSign(hexval, lengt)

test = 0;
if (test ==1)
    clc;close all;clear all;
    hexval={ '01451AC2', '005ACA6D', 'FF8237B5', '00855EFC', '00726593', 'FFA3CB58'};
    lengt = 8;    % 位宽
end

decval = hex2dec(hexval);
sign = bitget(decval, 4*lengt);             % 最高位为符号位，取出其符号位
negative_numbers = (sign == 1);              % 取为逻辑“真”或“假”
pvivn_decval = decval(negative_numbers);     % 取出decval中为负的数

decval(negative_numbers) = pvivn_decval - bitshift(1, 4*lengt);  % 转换为有符号数  -2^32
