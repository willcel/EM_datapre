
load('D:\0602����\0602����1����Ԥ����\1M�Ŵ���_���ٵ���\current_avg.mat')

tmp = current_avg_all(1:15, :);
current_avg_all = [ tmp];

save('current_avg.mat','current_avg_all')

