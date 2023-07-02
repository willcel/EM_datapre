
load('D:\0602测量\0602测线1数据预处理\1M放大器_减少叠加\current_avg.mat')

tmp = current_avg_all(1:15, :);
current_avg_all = [ tmp];

save('current_avg.mat','current_avg_all')

