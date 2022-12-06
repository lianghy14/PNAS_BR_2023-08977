clc
clear all
global area x y h nx ny
global n_OD boundary_O boundary_D boundary_D_2h boundary_D_outside
global boundary_H_1 boundary_H_2 boundary_H boundary_Pole boundary_Container
global sigma epsilon density_m density_0 v_free m_ave c0 c1 beta density_a b_c method panic_config d_min tao_n tao_p
global gfuns pfuns
global tstep t t_bound x_bound y_bound d_out t_panic pressure D_in
global panic_xy panic_ex x_ex y_ex nx_ex ny_ex x_g k_g
%% Paramters
h = 0.5 ; t = 750; tEnd = 900;

method = 'LFHE'; panic_config = '_panic'; plot_config = 1;
m_ave = 60; density_0 = 6; density_m = 10; tao = 0; tao_n = 5; tao_p = 0.3;
c0 = 0.6; beta = 0; density_a = 8; c1 = c0./2; b_c = 120; v_free = 1.034;
x_max = 105; y_max = 50;
area = [0,x_max;-2,y_max]; CFL = 0.2; tstep = 0.01; d_min = 0.01;
t_bound = [300 1200]; t_panic = [480, 600, 1200];
x_bound = [40,65]; y_bound = [y_max,y_max]; d_out = 6;
sigma = 10^-9; epsilon = 10^-6; gfuns = functions_given; pfuns = functions_plot;

% Computational domain
x = (area(1,1)+0.5*h):h:(area(1,2)-0.5*h); nx = length(x);
y = (area(2,2)-0.5*h):-h:(area(2,1)+0.5*h);  ny = length(y);
[X,Y] = meshgrid(x,y);

n_OD = 7;
boundary_O = cell(1,n_OD);
boundary_D = cell(1,n_OD);
boundary_D_2h = cell(1,n_OD);
boundary_D_outside = cell(1,n_OD);

% Origin boundaries
boundary_O{1} = {'Rectangle',1,[-5 0;0 20]};
boundary_O{2} = {'Rectangle',1,[x_max 0;x_max+5 20]};
boundary_O{3} = {'Rectangle',1,[40 y_max;65 y_max+5]};
boundary_O{4} = {'Rectangle',1,[40 y_max;65 y_max+5]};
boundary_O{5} = {'Rectangle',0};
boundary_O{6} = {'Rectangle',0};
boundary_O{7} = {'Rectangle',0};

% Destination boundaries
h_pole = 1;
boundary_D{1} = {'Rectangle',1,[40 y_max;65 y_max+5]};
boundary_D_2h{1} = {'Rectangle',1,[40 y_max-h;65 y_max+5]};
boundary_D_outside{1} = {'Rectangle',1,[40 y_max+h;65 y_max+5]};
boundary_D{2} = {'Rectangle',1,[40 y_max;65 y_max+5]};
boundary_D_2h{2} = {'Rectangle',1,[40 y_max-h;65 y_max+5]};
boundary_D_outside{2} = {'Rectangle',1,[40 y_max+h;65 y_max+5]};
boundary_D{3} = {'Rectangle',1,[-5 0;0 20]};
boundary_D_2h{3} = {'Rectangle',1,[-5 0;h 20]};
boundary_D_outside{3} = {'Rectangle',1,[-5 0;-h 20]};
boundary_D{4} = {'Rectangle',1,[x_max 0;x_max+5 20]};
boundary_D_2h{4} = {'Rectangle',1,[x_max-h 0;x_max+5 20]};
boundary_D_outside{4} = {'Rectangle',1,[x_max+h 0;x_max+5 20]};
% boundary_D{5} = {'Rectangle',1,[64 29;65 30]};
% boundary_D_2h{5} = {'Rectangle',1,[64 29;65 30]};
% boundary_D_outside{5} = {'Rectangle',1,[64 29;65 30]};
boundary_D{5} = {'Rectangle',1,[65-h_pole 30-h_pole;65 30]};
boundary_D_2h{5} = {'Rectangle',1,[65-h_pole 30-h_pole;65 30]};
boundary_D_outside{5} = {'Rectangle',1,[65-h_pole 30-h_pole;65 30]};
boundary_D{6} = {'Rectangle',1,[50 0;55 2]};
boundary_D_2h{6} = {'Rectangle',1,[50 0;55 2]};
boundary_D_outside{6} = {'Rectangle',1,[50 0;55 2]};
boundary_D{7} = {'Rectangle',1,[40 y_max;65 y_max+5]};
boundary_D_2h{7} = {'Rectangle',1,[40 y_max-h;65 y_max+5]};
boundary_D_outside{7} = {'Rectangle',1,[40 y_max+h;65 y_max+5]};

% Physical boundary
boundary_H = {'Rectangle',5,[-5 20;40 y_max+5],[65 20;x_max+5 y_max+5],[-5 -5;x_max+5 0],[65-h_pole 30-h_pole;65 30],[50 0;55 2]};
boundary_H_panic = {'Rectangle',3,[-5 20;40 y_max+5],[65 20;x_max+5 y_max+5],[-5 -5;x_max+5 0]};
boundary_Pole = {'Rectangle',1,[65-h_pole 30-h_pole;65 30]};
boundary_Container = {'Rectangle',1,[50 0;55 2]};

% Panic area
x_ex = (x(1)-2*h):h:(x(end)+2*h); nx_ex = length(x_ex);
y_ex = (y(1)+2*h):-h:(y(end)-2*h); ny_ex = length(y_ex);
panic_ex = zeros(ny_ex,nx_ex,n_OD);
panic_xy = zeros(ny,nx,n_OD);
panic_xy(:,:,5:6) = 1;
panic_ex(:,:,5:6) = 1;
judge_Pole = zeros(ny,nx);
judge_Container = zeros(ny,nx);
for i = 1:nx_ex
    for j = 1:ny_ex 
        dis_Pole = sqrt((x_ex(i) - 65 + h_pole./2).^2+(y_ex(j) - 30 + h_pole./2).^2)./(14.5+h_pole./2);
        dis_Container = sqrt(((x_ex(i) - 52.5)./17.5).^2+((y_ex(j) - 1)./16).^2);
%         for i_OD = 1:n_OD
%             panic_ex(j,i,i_OD) = max([0,min(1,(1-dis_Pole)./1),min(1,(1-dis_Container)./1)]);
%         end
%         panic_ex(j,i,5) = max([0,min(1,(1-dis_Pole)./1)]);
%         panic_ex(j,i,6) = max([0,min(1,(1-dis_Container)./1)]);
    end
end
for i = 1:nx
    for j = 1:ny
        dis_Pole = ((x(i) - 65 + h_pole./2)./(14.5+h_pole./2)).^2+((y(j) - 30 + h_pole./2)./(14.5+h_pole./2)).^2;
        dis_Container = ((x(i) - 52.5)./17.5).^2+((y(j) - 1)./16).^2;
        judge_Pole(j,i) = max([0,min(1,(1-dis_Pole)./1)]);
        judge_Container(j,i) = max([0,min(1,(1-dis_Container)./1)]);
%         for i_OD = 1:n_OD
%             panic_xy(j,i,i_OD) = max([0,min(1,(1-dis_Pole)./1),min(1,(1-dis_Container)./1)]);
%         end
%         panic_xy(j,i,5) = max([0,min(1,(1-dis_Pole)./1)]);
%         panic_xy(j,i,6) = max([0,min(1,(1-dis_Container)./1)]);
    end
end
x_n1 = sqrt(245-14*sqrt(70))/21;
x_n2 = sqrt(245+14*sqrt(70))/21;
k_n1 = (322+13.*sqrt(70))./900;
k_n2 = (322-13.*sqrt(70))./900;
x_g = [-x_n2; -x_n1; 0; x_n1; x_n2]./2 + 0.5;
k_g = [k_n2; k_n1; 128/225; k_n1; k_n2]./2;


%% IC for density potential
% dir_fig = ['C:\Users\HLiang\Desktop\Case study LP\' 'Multi-' method '_co' num2str(c0) '_dout' num2str(d_out) '_h' num2str(h) '_cd' num2str(density_0) '_pole' num2str(h_pole) '/'];
% dir_fig = [ 'Multi-' method '_co' num2str(c0) '_dout' num2str(d_out) '_h' num2str(h) '_cd' num2str(density_0) '_pole' num2str(h_pole) '/'];
dir_fig = ['C:\Users\HOWIE-PC\Desktop\Case study LP\' 'Multi-' method '_co' num2str(c0) '_dout' num2str(d_out) '_h' num2str(h) '_cd' num2str(density_0) '_pole' num2str(h_pole) '/'];
tstart = tic;

mkdir(dir_fig);
fprintf('Start Computing---------------------------------------------\n');
% Initial Condition
Q_cell = cell(1,n_OD);
if t == 0
    Q_cell_com = {zeros(ny,nx),zeros(ny,nx),zeros(ny,nx)};
    Q_cell_com{1} = gfuns.Boundary_value(x,y,Q_cell_com{1},boundary_H,0);
    for i = 1:n_OD
        Q_cell{i} = {zeros(ny,nx),zeros(ny,nx),zeros(ny,nx)};
        Q_cell{i}{1} = gfuns.Boundary_value(x,y,Q_cell{i}{1},boundary_H,0);
        for j = 1:3
            Q_cell_com{j} = Q_cell_com{j} + Q_cell{i}{j};
        end
    end
else
    Q_name = [dir_fig num2str(t) '.mat']; 
    load([dir_fig num2str(t) '.mat']);
end

Q_1 = cell(1,n_OD);
Q_2 = cell(1,n_OD);
Q_3 = cell(1,n_OD);
for i = 1:n_OD
    Q_1{i} = {zeros(ny,nx),zeros(ny,nx),zeros(ny,nx)};
    Q_2{i} = {zeros(ny,nx),zeros(ny,nx),zeros(ny,nx)};
    Q_3{i} = {zeros(ny,nx),zeros(ny,nx),zeros(ny,nx)};
end

alpha = 5;
%% LF scheme
while (t<=tEnd)

%% Add Panic flow    
    if t == t_panic(1)
        Q_old = Q_cell;
        for i = 1:3
            Q_cell{1}{i} = Q_old{1}{i}.*0;
            Q_cell{2}{i} = Q_old{2}{i}.*0;
            Q_cell{3}{i} = Q_old{3}{i}.*0;
            Q_cell{4}{i} = Q_old{4}{i}.*0;
            Q_cell{5}{i} = (Q_old{1}{i} + Q_old{2}{i} + Q_old{3}{i} + Q_old{4}{i}).*judge_Pole;
            Q_cell{6}{i} = (Q_old{1}{i} + Q_old{2}{i} + Q_old{3}{i} + Q_old{4}{i}).*judge_Container;
            Q_cell{7}{i} = (Q_old{1}{i} + Q_old{2}{i} + Q_old{3}{i} + Q_old{4}{i}) .* (1 - judge_Pole - judge_Container);
%             for j = 1:4
%                 Q_cell{j}{i} =Q_cell{j}{i}.*(1 - judge_Pole - judge_Container);
%             end
        end
    end

%% Plot&Save
    if(mod(t,1)==0)
        if (mod(t,1)==0) && (plot_config~=0)
            density_f = zeros(ny,nx,8);
            fig = figure(1);
            plotrange=[0,x_max,0,y_max,0,10];
            for i = 1:n_OD
                density_f(:,:,i) = gfuns.Boundary_value(x,y,Q_cell{i}{1},boundary_H,nan);
                subplot(5,2,i); mesh(X,Y,density_f(:,:,i)); colormap Jet; axis(plotrange);
                title([method ', dx = ',num2str(h),', dy = ',num2str(h),', time: ',num2str(t)]); xlabel('x(m)'); ylabel('y(m)'); zlabel('density(ped/m^2)');
            end
            density_com = gfuns.Boundary_value(x,y,Q_cell_com{1},boundary_H,nan);
            
            subplot(5,2,9); mesh(X,Y, density_com); colormap Jet; axis(plotrange);
            title([method ', dx = ',num2str(h),', dy = ',num2str(h),', time: ',num2str(t)]);
            xlabel('x(m)'); ylabel('y(m)'); zlabel('density(ped/m^2)');
            subplot(5,2,10); imagesc(x,y, density_com,'alphadata',~isnan(density_com)); colormap Jet; set(gca,'YDir','normal','color',0*[1 1 1]);caxis([0 10]); colorbar;
            axis([0,x_max,0,y_max]); xlabel('x(m)'); ylabel('y(m)');
            
            set(fig,'unit','centimeters','position',[10 5 20 18]);
            saveas(fig,[dir_fig num2str(t) '.png']);
        end
        tend = toc(tstart);
        fprintf([method '.t = %d. '],t)
        fprintf('Computing time is %f, Alpha is %f\n',tend,alpha);
        Q_name = [dir_fig num2str(t) '.mat'];
        save(Q_name,'Q_cell*');
        tstep = min([CFL*h/alpha, CFL]);
    else
        tstep = min([CFL*h/alpha, ceil(t)-t, CFL]);
    end

    [Q_cell,alpha] = TVD_RK3(Q_cell,t);
    
    Q_cell_com = {zeros(ny,nx),zeros(ny,nx),zeros(ny,nx)};
    for i = 1:3
        for j = 1:n_OD
            Q_cell{j}{i} = gfuns.Boundary_value(x,y,Q_cell{j}{i},boundary_H,0);
            Q_cell_com{i} = Q_cell_com{i} + Q_cell{j}{i};
        end
    end
    t = t+tstep;
end
