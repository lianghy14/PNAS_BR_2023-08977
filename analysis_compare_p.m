clc
clear
density_0 = 6;density_m = 10; area = [0,100;0,50];
CFL = 0.2; c0 = 1.2;  method = 'CWWENO3'; panic_config = '_panic'; plot_config = 1;

section_y = 35; section_x = 55; section_t = 155;
h_seq = 0.5; size = {'-','--'};
n_y = (50 - section_y)./h_seq +1;
n_x = section_x./h_seq +1;
dir_fig = ['C:\Users\HLiang\Desktop\Case study\' method '_CFL' num2str(CFL) '_h' num2str(h_seq) '_co' num2str(c0) '_panic/'];
dir_fig_np = ['C:\Users\HLiang\Desktop\Case study\' method '_CFL' num2str(CFL) '_h' num2str(h_seq) '_co' num2str(c0) '_no_panic/'];
Q_name = [dir_fig num2str(section_t) '.mat'];
Q_name_np = [dir_fig_np num2str(section_t) '.mat'];
load([dir_fig num2str(section_t) '.mat']);
section_x_density = Q_cell{1}(n_y,:);
section_y_density = Q_cell{1}(:,n_x);
load([dir_fig_np num2str(section_t) '.mat']);
section_np_x_density = Q_cell{1}(n_y,:);
section_np_y_density = Q_cell{1}(:,n_x);
fig1 = figure(1);
axes('Units','normalized','Position',[0.1 0.17 0.88 0.8]);
x = (area(1,1)):h_seq:(area(1,2));
plot(x,section_x_density,'k','LineStyle',size{1}); hold on
plot(x,section_np_x_density,'k','LineStyle',size{2});
legend('Panic Case','Normal case');xlabel('X (m)');ylabel('Density (ped/m^2)');
xlim([0 100]);
ylim([0 10]);
set(fig1,'Units','centimeters','Position',[10 10 12 6]);
hold off
fig2 = figure(2);
axes('Units','normalized','Position',[0.1 0.17 0.88 0.8]);
y = (area(2,2)):-h_seq:(area(2,1));
plot(y,section_y_density,'k','LineStyle',size{1}); hold on
plot(y,section_np_y_density,'k','LineStyle',size{2});
legend('Panic Case','Normal case');xlabel('Y (m)');ylabel('Density (ped/m^2)');
xlim([0 50]);
ylim([0 10]);
set(fig2,'Units','centimeters','Position',[10 10 12 6]);
hold off