clc
clear
global density_0 density_m v_free t

h = 0.5 ; tStart =  0 ; tEnd = 800;

method = 'LFHE'; panic_config = '_panic'; plot_config = 1;
m_ave = 60; density_0 = 6; density_m = 10; tao = 0; tao_n = 3; tao_p = 0.2;
c0 = 0.6; beta = 0; density_a = 8; c1 = c0.*density_0.^beta./(density_0-density_a); b_c = 120; v_free = 1.034;
x_max = 105; y_max = 50;
area = [0,x_max;-2,y_max]; CFL = 0.2; tstep = 0.01; d_min = 0.01;
t_bound = [300 1200]; t_panic = [480, 600, 1200];
x_bound = [40,65]; y_bound = [y_max,y_max]; d_out = 6;
sigma = 10^-9; epsilon = 10^-6; gfuns = functions_given; pfuns = functions_plot;

n_OD = 7;

FD_area = [30 55;20 40];
% Computational domain
x = (area(1,1)+0.5*h):h:(area(1,2)-0.5*h); nx = length(x);
y = (area(2,2)-0.5*h):-h:(area(2,1)+0.5*h);  ny = length(y);
[X,Y] = meshgrid(x,y);

n_area(1,1) = (FD_area(1,1) - area(1,1))./h + 1;
n_area(1,2) = (FD_area(1,2) - area(1,1))./h;
n_area(2,1) = (FD_area(2,1) - area(2,1))./h + 1;
n_area(2,2) = (FD_area(2,2) - area(2,1))./h;

dir_fig = [ 'Multi-' method '_co' num2str(c0) '_dout' num2str(d_out) '_h' num2str(h) '_cd' num2str(density_0) '/'];
% dir_fig = ['C:\Users\HLiang\Desktop\Case study LP\' 'Multi-' method '_co' num2str(c0) '_te' num2str(effi_t) '_dout' num2str(d_out) '_h' num2str(h) '_cd' num2str(density_0) '/'];
% dir_fig = ['C:\Users\HOWIE-PC\Desktop\Case study LP\' 'Multi-' method '_co' num2str(c0) '_te' num2str(effi_t) '_dout' num2str(d_out) '_h' num2str(h) '_cd' num2str(density_0) '/'];

density = cell(1,n_OD);
f_x = cell(1,n_OD);
f_y = cell(1,n_OD);
for i = 1:n_OD
    density{i} = zeros(tEnd-tStart+1,(n_area(2,2)-n_area(2,1)+1),(n_area(1,2)-n_area(1,1)+1));
    f_x{i} = zeros(tEnd-tStart+1,(n_area(2,2)-n_area(2,1)+1),(n_area(1,2)-n_area(1,1)+1));
    f_y{i} = zeros(tEnd-tStart+1,(n_area(2,2)-n_area(2,1)+1),(n_area(1,2)-n_area(1,1)+1));
end
density_com = zeros(tEnd-tStart+1,(n_area(2,2)-n_area(2,1)+1),(n_area(1,2)-n_area(1,1)+1));
f_x_com = zeros(tEnd-tStart+1,(n_area(2,2)-n_area(2,1)+1),(n_area(1,2)-n_area(1,1)+1));
f_y_com = zeros(tEnd-tStart+1,(n_area(2,2)-n_area(2,1)+1),(n_area(1,2)-n_area(1,1)+1));

for i = 1:(tEnd-tStart+1)
    Q_name = [dir_fig num2str(i+tStart-1) '.mat'];
    load([dir_fig num2str(i+tStart-1) '.mat']);
    for j = 1:n_OD
        density{j}(i,:,:) = Q_cell{j}{1}((ny-n_area(2,2)+1):(ny-n_area(2,1)+1),n_area(1,1):n_area(1,2));
        f_x{j}(i,:,:) = Q_cell{j}{2}((ny-n_area(2,2)+1):(ny-n_area(2,1)+1),n_area(1,1):n_area(1,2));
        f_y{j}(i,:,:) = Q_cell{j}{3}((ny-n_area(2,2)+1):(ny-n_area(2,1)+1),n_area(1,1):n_area(1,2));
        density_com(i,:,:) = density_com(i,:,:) + density{j}(i,:,:);
        f_x_com(i,:,:) = f_x_com(i,:,:) + f_x{j}(i,:,:);
        f_y_com(i,:,:) = f_y_com(i,:,:) + f_y{j}(i,:,:);
    end
end
% density_com = mean(mean(density_com,2),3);
% f_x_com = mean(mean(f_x_com,2),3);
% f_y_com = mean(mean(f_y_com,2),3);

density_com = reshape(density_com,[1,(tEnd-tStart+1)*(n_area(1,2)-n_area(1,1)+1)*(n_area(2,2)-n_area(2,1)+1)]);
f_x_com = reshape(f_x_com,[1,(tEnd-tStart+1)*(n_area(1,2)-n_area(1,1)+1)*(n_area(2,2)-n_area(2,1)+1)]);
f_y_com = reshape(f_y_com,[1,(tEnd-tStart+1)*(n_area(1,2)-n_area(1,1)+1)*(n_area(2,2)-n_area(2,1)+1)]);
v_x_com = f_x_com./density_com;
v_y_com = f_y_com./density_com;
v_x_com(isnan(v_x_com)) = 0; v_y_com(isnan(v_y_com)) = 0;
v_com = sqrt(v_x_com.^2+v_y_com.^2);

x = 0:0.01:density_m;
t = 0;
fe = x .* gfuns.Velosity(x,0);
t = t_panic(3);
fe2 = x .* gfuns.Velosity(x,1);
% fig1 = figure(1);
% axes('Units','normalized','Position',[0.08 0.16 0.88 0.8]);
% plot(x,fe); hold on;
% plot(density,flow,'.');
% set(fig1,'Units','centimeters','Position',[10 10 12 6]);

fig2 = figure(2);
axes('Units','normalized','Position',[0.09 0.18 0.88 0.8]);
dx = 0.1;
x_ave= 0:dx:10;density_ave= 1/2*dx:dx:(10-1/2*dx);
flow_ave= zeros(1,length(x_ave)-1);
count= zeros(1,length(x_ave)-1);
for i = 2:length(x_ave)
    for j = 1:length(density_com)
        if density_com(j)>=x_ave(i-1) && density_com(j)<x_ave(i)
            flow_ave(i-1) = v_com(j).*density_com(j)+flow_ave(i-1);
            count(i-1) = 1+count(i-1);
        end
    end
end
flow_ave=flow_ave./count;

plot(x,fe,'k'); hold on;
plot(x,fe2,'r'); hold on;
plot(density_ave,flow_ave,'k*','MarkerSize',2);
xlim([0 10]);ylim([0 3.0]);xlabel('Density (ped/m^2)');ylabel('Flow (ped/(m¡¤s))')
legend('FD1','FD2','Actual flow');
set(fig2,'Units','centimeters','Position',[10 10 12 6]);

fig3 = figure(3);
axes('Units','normalized','Position',[0.09 0.18 0.88 0.8]);
plot(x,fe,'k'); hold on;
plot(x,fe2,'r'); hold on;
plot(density_com,density_com.*v_com,'.');
xlim([0 10]); ylim([0 3.0]); xlabel('density(ped/m^2)');ylabel('Flow (ped/(m¡¤s))');
legend('FD1','FD2','Simulation')
set(fig3,'Units','centimeters','Position',[10 10 12 6]);


% fig4 = figure(4);
% axes('Units','normalized','Position',[0.08 0.15 0.88 0.8]);
% ve = gfuns.Velosity(x);
% plot(x,ve); hold on;
% plot(density,v,'.');
% set(fig4,'Units','centimeters','Position',[10 10 15 8]);