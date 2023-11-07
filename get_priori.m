load('depth_.mat')
load('cdi_rho_.mat')
load('mat_cdi_time.mat')
% close all
mat = mat_cdi;

%{
    kk = 6;
    figure(Position=[200 200 500 500])
    semilogx(cdi_rho_(kk,:),depth_(kk,:),'-ro','LineWidth',1.2)
    set(gca,'ydir','reverse')
    set(gca,'xdir','reverse')
    title(['测点',num2str(kk)])
    grid on
    xlabel('rho(ohm*m)')
    ylabel('Depth(m)')
%}


nolayer = 5;
% ns=24;

rho_pro = zeros(ns, nolayer);
dep_pro = zeros(ns, nolayer);
thrbank = ones(ns, nolayer+1);
for i = 1:ns
    thrbank(i,:) = [99	2	0.1 1e-2 1e-96  1e-97 ];
%     thrbank(i,:) = [99	1	0.05	3e-3 1e-96  1e-97 ];
    % thrbank(2,:) = [10	1	0.1	0.01	0.001 1e-6];
    % 个性化定制
    if ismember(i,[3,6,7])
%         thrbank(i,:) = [99	1	0.1 1e-2 1e-96  1e-97 ];
    end

end



for i=1:ns
    thr = thrbank(i,:);

    count = zeros(nolayer,1);
    max_rec = zeros(nolayer,1);
    
    for j=1:size(mat,2)

        depth = j/100;
        
        flag=0;
        for k=1:nolayer

            if(mat(i,j)<=thr(k) && mat(i,j) > thr(k+1))
                count(k) = count(k)+1;
                rho_pro(i,k) = rho_pro(i,k) +mat(i,j);
                dep_pro(i,k) = depth;
                
                flag=1;
                break
            end
        end
        
        if(flag==0)
            if mat(i,j) ~= 100
                1;
            end
            count(nolayer) = count(nolayer)+1; % 最后一层
            rho_pro(i,nolayer) = rho_pro(i,nolayer) +mat(i,j);
            dep_pro(i,nolayer) = depth;
        end
    end

    rho_pro(i,:) = rho_pro(i,:)./count';
   

end


rho_pro(isnan(rho_pro)) = 100;
dep_pro(dep_pro(:,nolayer-1)==0,nolayer-1) = total_depth-0.1; % 防止dep_pro1最后一列减为0

% use interpolation is better
% tmp_dep_pro = dep_pro;
% tmp_rho_pro = rho_pro;
% for i = 1:ns
%     for j = 2:nolayer-1
%         if ~ (tmp_dep_pro(i,j) >= dep_pro(i,j-1) && tmp_dep_pro(i,j) <= dep_pro(i,j+1))
%             tmp_dep_pro(i,j) = (dep_pro(i,j-1) + dep_pro(i,j+1)) /2;
%         end
%         if isnan(tmp_rho_pro(i,j)) % ~(rho_pro(i,j) >= rho_pro(i,j-1) && rho_pro(i,j) <= rho_pro(i,j+1)) || isnan(rho_pro(i,j))
%             tmp_rho_pro(i,j) = (rho_pro(i,j-1) + rho_pro(i,j+1)) /2;
%         end
%     end
% end
% dep_pro = tmp_dep_pro;
% rho_pro = tmp_rho_pro;


% rho_pro可能第一第二层为100
for i=1:size(rho_pro,1)
    for j = 1:4
        if rho_pro(i,j)==100
            tmp = rho_pro(i,1:4);
            idx = find(tmp<100,1);
            rho_pro(i,j) = tmp(idx);
        end
    end
end

%% 初始厚度都设置为5m试试
% for k = 1:ns
% dep_pro(k,:) = 3*(1:nolayer);
% end

dep_pro(:,1) = dep_pro(:,2) / 2;

scale_factor = 100;
%%
mat = zeros(ns, total_depth*scale_factor);

for x = 1:ns

    for y = 1:total_depth*scale_factor
        y1 = y/scale_factor;
        for nn=1:nolayer
            if(nn==1)
                if(y1<=dep_pro(x, nn))
                    mat(x,y)=rho_pro(x,nn);
                end
            else
                if(y1<=dep_pro(x, nn) && y1>dep_pro(x, nn-1))
                    mat(x,y)=rho_pro(x,nn);
                end
            end
            
        end
    end
end


%%

y = 0:0.01:total_depth-0.01;
xdraw_range = [pset, pset(end)+1]; mat = [mat;zeros(1,total_depth*scale_factor)];
figure(Position=[137	142.333333333333	991.333333333333	650.666666666667])
pcolor(delta_pset*(xdraw_range - min(xdraw_range)),y,log10(mat'))
shading flat
colormap jet
xlabel('X axis (m)','FontName','Calibri','FontSize',15,'FontWeight','bold')
ylabel('Z axis (m)','FontName','Calibri','FontSize',15,'FontWeight','bold')
h=colorbar;
set(get(h,'title'),'string','log10(\rho)');
caxis([-4,2])
title('Apparent resistivity imaging')
set(gca,'FontName','Calibri','FontSize',12,'FontWeight','bold')
set(gca,'ydir','reverse')
for i = 1:ns
%     scatter(i-0.25,1,'^')
        text(xdraw_range(i)-delta_pset*0.5, 3, num2str(i), ...
        'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'bottom', 'FontSize', 12);
end



%% save to txt
rho_pro1 = zeros(ns*nolayer,1);
dep_pro1 = zeros(ns*nolayer,1);

for i=1:size(rho_pro,1)
    if(i==14)
        1;
    end
    for j=1:size(rho_pro,2)
        
%         if(rho_pro(i,j)<=0.001)
%             rho_pro1(j+(i-1)*size(rho_pro,2)) = 10.0;
%         else
            rho_pro1(j+(i-1)*size(rho_pro,2)) = rho_pro(i,j);
%         end
        
        if(j==1)
            dep_pro1(j+(i-1)*size(rho_pro,2)) = dep_pro(i,j);
        else
            dep_pro1(j+(i-1)*size(rho_pro,2)) = dep_pro(i,j)-dep_pro(i,j-1);
        end
        
    end
end

dep_pro1(dep_pro1<=0) = 1;

save('rho_pro_tunnel_20ms.txt','rho_pro1','-ascii')
save('dep_pro_tunnel_20ms.txt','dep_pro1','-ascii')

%{
figure
imagesc(delta_pset*(pset-min(pset)),y,log10(mat'))
colormap jet
xlabel('X axis (m)','FontName','Calibri','FontSize',15,'FontWeight','bold')
ylabel('Z axis (m)','FontName','Calibri','FontSize',15,'FontWeight','bold')
h=colorbar;
set(get(h,'title'),'string','log10(\rho)');
% caxis([-5,5])
title('Apparent resistivity imaging')
set(gca,'FontName','Calibri','FontSize',12,'FontWeight','bold')
caxis([-4,2])
%}
