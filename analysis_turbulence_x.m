global density_m density_0 m_ave tao c0 beta method panic_config
tstart = tic;


h = 0.5; tStart = 1; tEnd = 900;

method = 'LFHE'; panic_config = '_panic'; plot_config = 1;
m_ave = 60; density_0 = 6; density_m = 10; tao = 0; tao_n = 1; tao_p = 1; effi_t = 0.3;
c0 = 0.3; beta = 1; density_a = 7.5; c1 = c0.*density_0.^beta./(density_0-density_a); b_c = 120; v_free = 1.034; thita = 0;
x_max = 100; y_max = 50;
area = [0,x_max;0,y_max+2]; CFL = 0.2; tstep = 0.01; d_min = 0.01;
t_bound = [200 tEnd]; t_panic = [200, 400, 900, 900];
x_bound = [x_max-2,x_max]; y_bound = [40,50]; d_out = 0;
sigma = 10^-9; epsilon = 10^-6; gfuns = functions_given; pfuns = functions_plot;

FD_area = [50 35;60 40];
n_area(1,1) = (FD_area(1,1) - area(1,1))./h + 1;
n_area(1,2) = (FD_area(1,2) - area(2,1))./h + 1;
n_area(2,1) = (FD_area(2,1) - area(1,1))./h;
n_area(2,2) = (FD_area(2,2) - area(2,1))./h;

t_p = t_bound(1):t_bound(2);
% Computational domain
x = (area(1,1)+0.5*h):h:(area(1,2)-0.5*h); nx = length(x);
y = (area(2,2)-0.5*h):-h:(area(2,1)+0.5*h);  ny = length(y);
[X,Y] = meshgrid(x,y);
pe = 3;

% Origin boundaries
n_OD = 2;
boundary_O = cell(1,n_OD);
boundary_D = cell(1,n_OD);
boundary_D_2h = cell(1,n_OD);
boundary_D_outside = cell(1,n_OD);

% Origin boundaries
boundary_O{1} = {'Rectangle',1,[-5 40;0 50]};
boundary_O{2} = {'Rectangle',1,[50 -5;60 0]};

% Destination boundaries
boundary_D{1} = {'Rectangle',1,[x_max 40;x_max+5 50]};
boundary_D_2h{1} = {'Rectangle',1,[x_max-h 40;x_max+5 50]};
boundary_D_outside{1} = {'Rectangle',1,[x_max+h 40;x_max+5 50]};
boundary_D{2} = {'Rectangle',1,[x_max 40;x_max+5 50]};
boundary_D_2h{2} = {'Rectangle',1,[x_max-h 40;x_max+5 50]};
boundary_D_outside{2} = {'Rectangle',1,[x_max+h 40;x_max+5 50]};


% Physical boundary
boundary_H = {'Rectangle',3,[-5 -5;50 40],[-5 y_max;x_max+5 y_max+5],[60 -5;x_max+5 40]};
boundary_H_bound = {'Rectangle',4,[x_max-5 43;x_max-2 47],[-5 -5;50 40],[-5 y_max;x_max+5 y_max+5],[60 -5;x_max+5 40]};

% Panic area

x_ex = (x(1)-2*h):h:(x(end)+2*h); nx_ex = length(x_ex);
y_ex = (y(1)+2*h):-h:(y(end)-2*h); ny_ex = length(y_ex);
panic_ex = zeros(ny_ex,nx_ex);
for i = 1:nx_ex
    for j = 1:ny_ex
        panic_ex(j,i) = max(0,min(1,(90-x_ex(i))./20));
    end
end
for i = 1:nx
    for j = 1:ny
        panic_xy(j,i) = max(0,min(1,(90-x(i))./20));
    end
end


turbu = zeros(tEnd-tStart+1);
density = zeros((n_area(2,2)-n_area(1,2)+1),(n_area(2,1)-n_area(1,1)+1), tEnd-tStart+1);
v = zeros((n_area(2,2)-n_area(1,2)+1),(n_area(2,1)-n_area(1,1)+1), tEnd-tStart+1);
f_x = zeros((n_area(2,2)-n_area(1,2)+1),(n_area(2,1)-n_area(1,1)+1), tEnd-tStart+1);
f_y = zeros((n_area(2,2)-n_area(1,2)+1),(n_area(2,1)-n_area(1,1)+1), tEnd-tStart+1);


dir_fig = ['C:\Users\HLiang\Desktop\Case study HJ\' 'Multi-' method '_co' num2str(c0) '_te' num2str(effi_t) '_dout' num2str(d_out) '_h' num2str(h) '_cd' num2str(density_0) '/'];
% dir_fig = ['C:\Users\HOWIE-PC\Desktop\Case study HJ\' 'Multi-' method '_co' num2str(c0) '_te' num2str(effi_t) '_dout' num2str(d_out) '_h' num2str(h) '_cd' num2str(density_0) '/'];

for i = (tStart+1):(tEnd+1)
    load([dir_fig num2str(i-1) '.mat']);
    Q_cell_com = {zeros(ny,nx),zeros(ny,nx),zeros(ny,nx)};
    Q_cell_com{1} = gfuns.Boundary_value(x,y,Q_cell_com{1},boundary_H,0);
    for j = 1:n_OD
        for k = 1:3
            Q_cell_com{k} = Q_cell_com{k} + Q_cell{j}{k};
        end
    end
    density(:,:,i) = Q_cell_com{1}((ny-n_area(2,2)+1):(ny-n_area(1,2)+1),n_area(1,1):n_area(2,1));
    f_x(:,:,i) = Q_cell_com{2}((ny-n_area(2,2)+1):(ny-n_area(1,2)+1),n_area(1,1):n_area(2,1));
    f_y(:,:,i) = Q_cell_com{3}((ny-n_area(2,2)+1):(ny-n_area(1,2)+1),n_area(1,1):n_area(2,1));
    v(:,:,i) = sqrt(f_x(:,:,i).^2+f_y(:,:,i).^2)./max(10^-4,abs(density(:,:,i)));
end
for i = (tStart+1):(tEnd+1)
    v_var = var(v(:,:,i),0,'all');
    density_ave = mean(mean(density(:,:,i)));
    turbu(i-1) = density_ave.*v_var;
end
fig = figure();
plot(t_p,turbu(t_p));
xlabel('Time (s)');ylabel('Pressure(var of t (/s^2)')
set(fig,'Units','centimeters','Position',[10 10 15 8]);

% t_plot = 350;
% turbu_plot = gfuns.Boundary_value(x,y,turbu(:,:,t_plot),boundary_H,nan);
% fig = figure();
% imagesc(x,y,turbu_plot,'alphadata',~isnan(turbu_plot)); colormap Jet; set(gca,'YDir','normal','color',0*[1 1 1]);colorbar;
% xlabel('x(m)'); ylabel('y(m)');
% title([method ', dx = ',num2str(h),', dy = ',num2str(h),', time: ',num2str(t_plot)])
% set(fig,'unit','centimeters','position',[10 5 15 10]);

TP_max = max(turbu(t_p));
tend = toc(tstart);
fprintf('Computing time is %f, D_out is %f, TE is %f\n',tend, d_out, effi_t);





