
load('point1set.txt')
load('point4set.txt')
close all
t0 = point4set(1)-point1set(1);
% t0 = 0.002;
% t0 = t_st;

load('signal_sample.mat')

path_cdi = '.\CDI_code';

depth_ = [];
cdi_rho_ = [];


delta_t = (log10(t_ed)-log10(t_st))/(nt-1);
tlog = zeros(1,nt);
t = zeros(1,nt);

for i=1:nt
    tlog(i)=log10(t_st)+(i-1)*delta_t;
    t(i)=10^tlog(i);
end

obs_original = load([path_cdi,'\log1.dat']);
rho = load([path_cdi,'\rho1.dat']);

is_plot = 0;
for uu=1:ns
%     uu = 17 ;%
    a = signal_sample(uu,:);
    
    if uu==6
        1;
%         is_plot = 1;
    end
    obs = obs_original;
    
    
    if(is_plot==1)
        figure%(Position=[200 200 200 200])
        for i=1:1:nt
            loglog(rho,obs(:,i),'LineWidth',0.5)
            hold on
        end
        xlabel('\rho \Omega/m')
        ylabel('voltage (V)')
        grid on
        title(['测点',num2str(uu)])
        set(gca,'FontSize',12,'FontWeight','bold')
    end
    
    %
    delta = zeros(1,nt);
    cdi_rho = zeros(1,nt);
    mu0 = 4*pi*10^(-7);
    
%     ntstop = [19 16 20 25 39 40 43 34 27 22 24 20 23 25 35 28 30];
    ntstop = [27 25 34 29 37 42 44 32 31 22 34 27 34 34 33 34 34];
             % 1 /2/ 3/4 /5 /6/ 7 /8 /9 /10/11/12/13/14/15/16/17
    for i=1:ntstop(uu)%nt
        
        minerr = 1000;
        [~,midIndex] = max(obs(:,i)); % debug钰舒bug，查表查到左边去了
        for j=1:length(rho)
            err = abs(obs(j,i)-a(i));
            if(err<minerr && j > midIndex)
                minerr=err;
                index = j;
            end
        end
        
        if(rho(index)>100)
            cdi_rho(i) = 100;
        else
            cdi_rho(i) = rho(index);
        end
        
        delta(i) = (2*(t(i)-t0)*cdi_rho(i)/mu0)^0.5;  % 这里减去的是关断的时间
        
        if(is_plot==1)
            scatter(cdi_rho(i),a(i),'^')
        end
    end
    

    c = 0.7;
    b2 = 0.4;

    delta1 = delta;
    delta1(1) = c*delta(1);
    
    for i=2:length(delta1)
        if delta(i) ~= 0  % 如果只选择部分抽道时刻，那么就不计算超过此时刻的插值
            tmp = delta(i)-delta(i-1);
        else
            tmp = 0;
        end
        if(tmp<0)
            b = -b2;
        else
            b = b2;
        end
        delta1(i) = delta1(i-1)+abs(b*tmp);
    end
    
    
    depth_ = [depth_; delta1];
    cdi_rho_ = [cdi_rho_; cdi_rho];
end


%%
nolayer = nt;         % 层数

scale_factor = 100;


%%
mat = 100.*ones(ns, total_depth*scale_factor);

for x = 1:ns
    for y = 1:total_depth*scale_factor
        y1 = y/scale_factor;
        for nn=1:nolayer
            if(nn==1)
                if(y1<=depth_(x, nn))
                    mat(x,y)=cdi_rho_(x,nn);
                end
            else
                if(y1<=depth_(x, nn) && y1>depth_(x, nn-1))
                    mat(x,y)=cdi_rho_(x,nn);
                end
            end
            
        end
    end
end


dy = 1/scale_factor;
y = 0:dy:total_depth-dy;
% y = total_depth-dy:-dy:0;
pset = 1:ns+1; mat = [mat;zeros(1,total_depth*scale_factor)];
% figure(Position=[1221	188.333333333333	1244	850.666666666667]) 
figure(Position=[1221	313.666666666667	784	625.333333333333]) 

%  get(gcf,'Position')

% imagesc(delta_pset*(pset-min(pset)),y,log10(mat'))

pcolor(delta_pset*(pset),y,log10(mat'))
shading flat

% pcolor(delta_pset*(pset-min(pset)),y,log10(mat'))
% shading interp
% xlim([0 8])

colormap jet
xlabel('Measurement Line / m','FontSize',15,'FontWeight','bold')
ylabel('Depth / m','FontSize',15,'FontWeight','bold')

h=colorbar;
set(get(h,'title'),'string','log10(\rho)');
caxis([-4,2])
% caxis([-3,2])
% title('Apparent resistivity imaging')
set(gca,'FontSize',18,'FontWeight','bold')
% set(gca,'yticklabel',{'10','8','6','4','2'});
set(gca,'ydir','reverse')
hold on
for i = 1:ns
%     scatter(i-0.25,1,'^')
        text(i*0.5+0.25, 3, num2str(i), ...
        'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'bottom', 'FontSize', 12);
end
%{
    figure
    ns1 = mat(1,:);
    loglog(ns1,y,'LineWidth',1.2)
    set(gca,'ydir','reverse')
    title('1')
    grid on
    xlabel('rho(ohm*m)')
    ylabel('Depth(m)')
    figure
    ns2 = mat(2,:);
    loglog(ns2,y,'LineWidth',1.2)
    set(gca,'ydir','reverse')
    title('2')
    grid on
    xlabel('rho(ohm*m)')
    ylabel('Depth(m)')
%}

% {
save('depth_.mat','depth_')
save('cdi_rho_.mat','cdi_rho_')
mat_cdi = mat;
save('mat_cdi_time.mat','mat_cdi')
%}

