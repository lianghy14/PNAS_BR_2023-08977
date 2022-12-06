clc
clear
m_ave = 60; density_0 = 3.5; density_m = 7; tao = 0.5; c0 = 1.2; v_free = 1.34; fric_a = 90; fric_b = 6;
x_max = 50; y_max = 40;
area = [0,x_max;0,y_max]; CFL = 0.2; method = 'CWWENO3'; panic_config = '_panic'; plot_config = 1;

section_y = 22; section_x = 28; section_t = 180;
h_seq = [0.5]; h_size = {'-'};h_legend = cell(1,1);
section_x_density = cell(1,length(h_seq));section_y_density = cell(1,length(h_seq));
for i = 1:length(h_seq)
    dir_fig = ['C:\Users\HOWIE-PC\Desktop\Case study\' 'Multi-' method '_co' num2str(c0) panic_config '_CFL' num2str(CFL) '_h' num2str(h_seq(i)) '_fric_a' num2str(fric_a) '_b' num2str(fric_b) '/'];
    Q_name = [dir_fig num2str(section_t) '.mat'];
    load([dir_fig num2str(section_t) '.mat']);
    n_y = (50 - section_y)./h_seq(i) +1;
    n_x = section_x./h_seq(i) +1;
    section_x_density{i} = Q_cell_com{1}(n_y,:);
    section_y_density{i} = Q_cell_com{1}(:,n_x);
    h_legend{i} = ['h = ',num2str(h_seq(i))];
end
fig1 = figure(1);
% axes('Units','normalized','Position',[0.09 0.16 0.88 0.8]);
for i = 1:length(h_seq)
    x = (area(1,1)):h_seq(i):(area(1,2));
    plot(x,section_x_density{i},'k','LineStyle',h_size{i});
    hold on
end
legend(h_legend{1});
xlabel('X (m)');ylabel('Density (ped/m^2)');
xlim([area(1,1) area(1,2)]);ylim([0 density_m]);
set(fig1,'Units','centimeters','Position',[10 10 12 6]);
hold off

fig2 = figure(2);
% axes('Units','normalized','Position',[0.09 0.16 0.88 0.8]);
for i = 1:length(h_seq)
    y = (area(2,2)):-h_seq(i):(area(2,1));
    plot(y,section_y_density{i},'k','LineStyle',h_size{i});
    hold on
end
legend(h_legend{1});
xlabel('Y (m)');ylabel('Density (ped/m^2)');
xlim([area(2,1) area(2,2)]);ylim([0 density_m]);
set(fig2,'Units','centimeters','Position',[10 10 12 6]);
hold off
