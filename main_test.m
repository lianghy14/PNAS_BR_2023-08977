clc
clear all
global area x y h nx ny
global n_OD boundary_O boundary_D boundary_D_2h boundary_D_outside
global boundary_H_1 boundary_H_2 boundary_H boundary_Pole boundary_Container 
global sigma epsilon density_m density_0 v_free m_ave c0 c1 beta density_a b_c method panic_config d_min tao_n tao_p effi_t thita
global gfuns pfuns
global tstep t t_bound x_bound y_bound d_out t_panic
global panic_xy panic_ex x_ex y_ex nx_ex ny_ex x_g k_g
%% Paramters
h = 2 ; t =  0; tEnd = 300;

method = 'LFHE'; panic_config = '_panic'; plot_config = 1;
m_ave = 60; density_0 = 6; density_m = 10; tao = 0; tao_n = 1; tao_p = 0.2; effi_t = 0.2;
c0 = 0.1; beta = 1; density_a = 7.5; c1 = c0.*density_0.^beta./(density_0-density_a); b_c = 120; v_free = 1.034; thita = 0;
x_max = 120; y_max = 30;
area = [0,x_max;-2,y_max+2]; CFL = 0.2; tstep = 0.01; d_min = 0.01;
t_bound = [tEnd tEnd]; t_panic = [900, 900, 900];
x_bound = [30,55]; y_bound = [y_max,y_max]; d_out = 0;
sigma = 10^-9; epsilon = 10^-6; gfuns = functions_given; pfuns = functions_plot;

% Computational domain
x = (area(1,1)+0.5*h):h:(area(1,2)-0.5*h); nx = length(x);
y = (area(2,2)-0.5*h):-h:(area(2,1)+0.5*h);  ny = length(y);
[X,Y] = meshgrid(x,y);

n_OD = 2;
boundary_O = cell(1,n_OD);
boundary_D = cell(1,n_OD);
boundary_D_2h = cell(1,n_OD);
boundary_D_outside = cell(1,n_OD);

% Origin boundaries
boundary_O{1} = {'Rectangle',1,[-5 0;0 y_max]};
boundary_O{2} = {'Rectangle',1,[x_max 0;x_max+5 y_max]};

% Destination boundaries
boundary_D{1} = {'Rectangle',1,[x_max 0;x_max+5 y_max]};
boundary_D_2h{1} = {'Rectangle',1,[x_max-h 0;x_max+5 y_max]};
boundary_D_outside{1} = {'Rectangle',1,[x_max+h 0;x_max+5 y_max]};
boundary_D{2} = {'Rectangle',1,[-5 0;0 y_max]};
boundary_D_2h{2} = {'Rectangle',1,[-5 0;h y_max]};
boundary_D_outside{2} = {'Rectangle',1,[-5 0;-h y_max]};


% Physical boundary
boundary_H = {'Rectangle',2,[-5 -5;x_max+5 0],[-5 y_max;x_max+5 y_max+2]};

% Panic area
x_ex = (x(1)-2*h):h:(x(end)+2*h); nx_ex = length(x_ex);
y_ex = (y(1)+2*h):-h:(y(end)-2*h); ny_ex = length(y_ex);
panic_xy = zeros(ny,nx);
panic_ex = zeros(ny_ex,nx_ex);

x_n1 = sqrt(245-14*sqrt(70))/21;
x_n2 = sqrt(245+14*sqrt(70))/21;
k_n1 = (322+13.*sqrt(70))./900;
k_n2 = (322-13.*sqrt(70))./900;
x_g = [-x_n2; -x_n1; 0; x_n1; x_n2]./2 + 0.5;
k_g = [k_n2; k_n1; 128/225; k_n1; k_n2]./2;


%% IC for density potential
dir_fig = ['C:\Users\HLiang\Desktop\Case study Multi\' 'Multi-' method '_co' num2str(c0) '_te' num2str(effi_t) '_dout' num2str(d_out) '_h' num2str(h) '_cd' num2str(density_0) '/'];
% dir_fig = ['C:\Users\HOWIE-PC\Desktop\Case study Multi\' 'Multi-' method '_co' num2str(c0) '_te' num2str(effi_t) '_dout' num2str(d_out) '_h' num2str(h) '_cd' num2str(density_0) '/'];
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

%% Plot&Save
    if(mod(t,1)==0)
        if (mod(t,10)==0) && (plot_config~=0)
            density_f = zeros(ny,nx,2);
            fig = figure(1);
            plotrange=[0,x_max,0,y_max,0,10];
            for i = 1:n_OD
                density_f(:,:,i) = gfuns.Boundary_value(x,y,Q_cell{i}{1},boundary_H,nan);
                subplot(2,2,i); mesh(X,Y,density_f(:,:,i)); colormap Jet; axis(plotrange);
                title([method ', dx = ',num2str(h),', dy = ',num2str(h),', time: ',num2str(t)]); xlabel('x(m)'); ylabel('y(m)'); zlabel('density(ped/m^2)');
            end
            density_com = gfuns.Boundary_value(x,y,Q_cell_com{1},boundary_H,nan);
            
            subplot(2,2,3); mesh(X,Y, density_com); colormap Jet; axis(plotrange);
            title([method ', dx = ',num2str(h),', dy = ',num2str(h),', time: ',num2str(t)]);
            xlabel('x(m)'); ylabel('y(m)'); zlabel('density(ped/m^2)');
            subplot(2,2,4); imagesc(x,y, density_com,'alphadata',~isnan(density_com)); colormap Jet; set(gca,'YDir','normal','color',0*[1 1 1]);caxis([0 10]); colorbar;
            axis([0,x_max,0,y_max]); xlabel('x(m)'); ylabel('y(m)');
            
            set(fig,'unit','centimeters','position',[10 5 20 10]);
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
