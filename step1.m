% �����ļ���·��
% clear
dbstop if error
folder_path = 'D:\willcel\0109����\1.9��ɽ\1.9��ɽ30m������';
needExaminDetail = 1;

data_avg_all = [];
current_avg_all = [];
for k = 1:ns
    k
    % ��ȡ�ļ��������е� txt �ļ���
    volt_files = dir(fullfile(folder_path,num2str(k), '*_192.168.2.80_0_*.txt'));
    curr_files = dir(fullfile(folder_path,num2str(k), '*_192.168.2.80_1_*.txt'));

    % ��ÿ�� CSV �ļ����д���
    data_offset = [];
    current_offset = [];

    for i = 1:numel(volt_files)
        i
        if ( (k==1 && i==1) || (k==3 && i==5) ||  (k==6 && i==3))
%             continue
        end
        % ��ȡ�ļ�     
        volt = read2(volt_files, i);
        curr = read2(curr_files, i);

        % ����ļ��ǿ�
        stidx = floor( length(volt) * 3 / 4 );
        offset = mean(volt(stidx:end));                % ������Ȧ��ƫ��
        offset_current = mean(curr(stidx:end));    % ���������ƫ��
        
        % ����ƫ��            
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
    %------------------- ȡƽ�� ---------------------
    if(size(data_offset, 1) > 1)
        data_avg = mean(data_offset);
        current_avg = mean(current_offset);
    else
        data_avg = (data_offset);
        current_avg = (current_offset);
    end
    data_avg = data_avg.*factor;                            % ���ԷŴ��� V
    current_avg = current_avg.*factor_current;      % ����˥������  A
    
    data_avg_all = [data_avg_all; data_avg];
    current_avg_all = [current_avg_all; current_avg];
    
    if needExaminDetail && k == 7
%         data_offset = data_offset .* factor; 
%         save('data_all_ns7.mat', 'data_offset')
    end
end

save('data_avg.mat','data_avg_all')
% current_avg_all = -current_avg_all; % ��֪Ϊʲô�����̽��
save('current_avg.mat','current_avg_all')

%%
step2_pre % ������
%%
step2     % �洢������ת�۵�

