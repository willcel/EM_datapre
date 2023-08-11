% �����ļ���·��
% clear
dbstop if error
folder_path = 'D:\willcel\0629�����\����4';
needExaminDetail = 1;


data_avg_all = [];
current_avg_all = [];
for k = 1:ns
    k
    % ��ȡ�ļ��������е� txt �ļ���
    txt_files = dir(fullfile(folder_path,['���',num2str(k)],'save_time_data', '*.save_data_time.txt'));
    
    % ��ÿ�� CSV �ļ����д���
    data_offset = [];
    current_offset = [];

    for i = 1:numel(txt_files)
        i
        if (k==5 && i==10) || (k==10 && i==38)
            continue
        end
        % ��ȡ CSV �ļ�
        % filename = fullfile(folder_path, ['���',num2str(k)],txt_files(i).name)
        % ��Ҫ�ϸ�ѡǰ20��
        filename = fullfile(folder_path, ['���',num2str(k)],'save_time_data',[num2str(i),'.save_data_time.txt']);

        [A,B,C] = textread(filename,'%s %s %s', 'headerlines', 1);   % ��ȡʮ����������
        
        datanew = zeros(length(B), 3);
        
        % datanew(:,1)=hex2decWithSign(A, 8)/2^23;  % 24λAD
        datanew(:,2)=hex2decWithSign(B, 8)/2^23;
        datanew(:,3)=hex2decWithSign(C, 8)/2^23;

        % ����ļ��ǿ�
        offset = mean(datanew(20000:end,2));                % ������Ȧ��ƫ��
        offset_current = mean(datanew(20000:end,3));    % ���������ƫ��
        
        % ����ƫ��            
        current = datanew(:,3)-offset_current;
        signal = datanew(:,2)-offset;
        % current = datanew(:,2);
        % signal = datanew(:,3);
        is_plot = 0;
        if(is_plot)
            dt = 1/fs;
            time = (1:length(datanew(:,2))).*dt;
            
            figure
            subplot(2,1,1)
            plot(time*1000, signal.*factor,'LineWidth',2.0)
            set(gca,'FontName','Calibri','FontSize',12,'FontWeight','bold')
            xlabel('time(ms)')
            ylabel('voltage (V)')
            grid on
            xlim([0 5])
%             ylim([0,3]*1e-4)
            
            subplot(2,1,2)
            plot(time*1000, current.*factor_current,'LineWidth',2.0)                
            xlim([0,5])
            xlabel('time(ms)')
            ylabel('current (A)')
            grid on
            set(gca,'FontName','Calibri','FontSize',12,'FontWeight','bold')
        end
%             data_array = [data_array; file_data];
        data_offset = [data_offset; signal'];
        current_offset = [current_offset; current'];
    end
    %------------------- ȡƽ�� ---------------------
    data_avg = mean(data_offset);
    current_avg = mean(current_offset);
    
    data_avg = data_avg.*factor;                            % ���ԷŴ��� V
    current_avg = current_avg.*factor_current;      % ����˥������  A
    
    data_avg_all = [data_avg_all; data_avg];
    current_avg_all = [current_avg_all; current_avg];
    
    if needExaminDetail && k == 7
        data_offset = data_offset .* factor; 
        save('data_all_ns7.mat', 'data_offset')
    end
end

save('data_avg.mat','data_avg_all')
% current_avg_all = -current_avg_all; % ��֪Ϊʲô�����̽��
save('current_avg.mat','current_avg_all')

