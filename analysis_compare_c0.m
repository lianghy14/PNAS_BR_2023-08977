clc
clear
density_0 = 4;density_m = 10; area = [0,100;0,50];
CFL = 0.2; c0_seq = [1.0 4.0]; method = 'WENO3'; panic_config = '_panic'; plot_config = 1;

section_y = 40; section_x = 55; section_t = 150;
h = 1; c0_size = ['*','o'];c0_legend = cell(1,2);
section_x_density = cell(1,length(c0_seq));section_y_density = cell(1,length(c0_seq));
for i = 1:length(c0_seq)
    dir_fig = ['C:\Users\HLiang\Desktop\Case study\' method '_CFL' num2str(CFL) '_h' num2str(h) '_co' num2str(c0_seq(i)) '_panic/'];
    Q_name = [dir_fig num2str(section_t) '.mat'];
    load([dir_fig num2str(section_t) '.mat']);
    n_y = (50 - section_y)./h +1;
    n_x = section_x./h +1;
    section_x_density{i} = Q_cell{1}(n_y,:);
    section_y_density{i} = Q_cell{1}(:,n_x);
    c0_legend{i} = ['c_0 = ',num2str(c0_seq(i))];
end
figure
for i = 1:length(c0_seq)
    x = (area(1,1)):h:(area(1,2));
    plot(x,section_x_density{i},'-','Marker',c0_size(i));
    hold on
end
legend(c0_legend{1},c0_legend{2});
xlabel('x(m)');ylabel('density(ped/m^2)');
xlim([0 100]);
hold off

figure
for i = 1:length(c0_seq)
    y = (area(2,2)):-h:(area(2,1));
    plot(y,section_y_density{i},'-','Marker',c0_size(i));
    hold on
end
legend(c0_legend{1},c0_legend{2});
xlabel('y(m)');ylabel('density(ped/m^2)');
xlim([0 50]);
hold off