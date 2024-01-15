% 定义文件夹路径
% clear
dbstop if error
folder_path = 'D:\willcel\0109测试\1.9汤山\1.9汤山30m新阻尼';
needExaminDetail = 1;

data_avg_all = [];
current_avg_all = [];
for k = 1:ns
    k
    % 获取文件夹内所有的 txt 文件名
    volt_files = dir(fullfile(folder_path,num2str(k), '*_192.168.2.80_0_*.txt'));
    curr_files = dir(fullfile(folder_path,num2str(k), '*_192.168.2.80_1_*.txt'));

    % 对每个 CSV 文件进行处理
    data_offset = [];
    current_offset = [];

    for i = 1:numel(volt_files)
        i
        if ( (k==1 && i==1) || (k==3 && i==5) ||  (k==6 && i==3))
%             continue
        end
        % 读取文件     
        volt = read2(volt_files, i);
        curr = read2(curr_files, i);

        % 如果文件非空
        stidx = floor( length(volt) * 3 / 4 );
        offset = mean(volt(stidx:end));                % 接收线圈的偏置
        offset_current = mean(curr(stidx:end));    % 发射电流的偏置
        
        % 减掉偏置            
        current = curr-offset_current;
        signal = volt-offset;
        % current = datanew(:,2);
        % signal = datanew(:,3);
        
        if(0)
            %%
            dt = 1/fs;
            time = (1:length(signal)).*dt;
            
            figure
            subplot(2,1,1)
            plot(time*1000, signal.*factor,'LineWidth',2.0)
            set(gca,'FontName','Calibri','FontSize',12,'FontWeight','bold')
            xlabel('time(ms)')
            ylabel('voltage (V)')
            grid on
            xlim([0 5])
            title(['cedian',num2str(k)])
%             xlim([0 30])
%             ylim([-1e-3 1e-3])
            
%             svdata = signal.*factor;
%             save data_step svdata time
%             ylim([0,3]*1e-4)
                
            subplot(2,1,2)
            plot(time*1000, current.*factor_current,'LineWidth',2.0)                
            xlim([0,5])
            xlabel('time(ms)')
            ylabel('current (A)')
            grid on
            set(gca,'FontName','Calibri','FontSize',12,'FontWeight','bold')
            
%             name1 = ['ns',num2str(k),'_measure',num2str(i)];
%             svfig(name1)
%             close all

            % find(time > 0.58e-3, 1) = 149
        end
%             data_array = [data_array; file_data];
        data_offset = [data_offset; signal'];
        current_offset = [current_offset; current'];
    end
    %------------------- 取平均 ---------------------
    if(size(data_offset, 1) > 1)
        data_avg = mean(data_offset);
        current_avg = mean(current_offset);
    else
        data_avg = (data_offset);
        current_avg = (current_offset);
    end
    data_avg = data_avg.*factor;                            % 除以放大倍数 V
    current_avg = current_avg.*factor_current;      % 乘以衰减倍数  A
    
    data_avg_all = [data_avg_all; data_avg];
    current_avg_all = [current_avg_all; current_avg];
    
    if needExaminDetail && k == 7
%         data_offset = data_offset .* factor; 
%         save('data_all_ns7.mat', 'data_offset')
    end
end

save('data_avg.mat','data_avg_all')
% current_avg_all = -current_avg_all; % 不知为什么，工程解决
save('current_avg.mat','current_avg_all')

%%
step2_pre % 零点操作
%%
step2     % 存储电流的转折点

