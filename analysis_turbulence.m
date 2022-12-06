global area x y h nx ny
global n_OD boundary_O boundary_D boundary_D_2h boundary_D_outside
global boundary_H_1 boundary_H_2 boundary_H boundary_Pole boundary_Container
global sigma epsilon density_m density_0 v_free m_ave c0 c1 beta density_a b_c method panic_config d_min tao_n tao_p
global gfuns pfuns
global tstep t t_bound x_bound y_bound d_out t_panic pressure D_in i_OD
global panic_xy panic_ex x_ex y_ex nx_ex ny_ex x_g k_g

tstart = tic;

h = 0.5 ; tStart = 0; tEnd = 820;

method = 'LFHE'; panic_config = '_panic'; plot_config = 1;
m_ave = 60; density_0 = 6; density_m = 10; tao = 0; tao_n = 5; tao_p = 0.3;
c0 = 0.6; beta = 0; density_a = 8; c1 = c0./2; b_c = 120; v_free = 1.034;
x_max = 105; y_max = 50;
area = [0,x_max;-2,y_max]; CFL = 0.2; tstep = 0.01; d_min = 0.01;
t_bound = [300 1200]; t_panic = [480, 600, 1200];
x_bound = [40,65]; y_bound = [y_max,y_max]; d_out = 6;
sigma = 10^-9; epsilon = 10^-6; gfuns = functions_given; pfuns = functions_plot;

Turbu_area = [40,65;15,45];
n_area(1,1) = (Turbu_area(1,1) - area(1,1))./h + 1;
n_area(1,2) = (Turbu_area(1,2) - area(1,1))./h;
n_area(2,1) = (area(2,2) - Turbu_area(2,2))./h + 1;
n_area(2,2) = (area(2,2) - Turbu_area(2,1))./h;
radius = 1.5; nr = (radius-h)/h/2;
pe = 2;

% Computational domain
x = (area(1,1)+0.5*h):h:(area(1,2)-0.5*h); nx = length(x);
y = (area(2,2)-0.5*h):-h:(area(2,1)+0.5*h);  ny = length(y);
[X,Y] = meshgrid(x,y);
h_i = 0.1;
x_i = (area(1,1)+0.5*h_i):h_i:(area(1,2)-0.5*h_i);
y_i = (area(2,2)-0.5*h_i):-h_i:(area(2,1)+0.5*h_i);
[X_i,Y_i] = meshgrid(x_i,y_i);

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




turbu = zeros(ny,nx,tEnd-tStart+1);
density_com = zeros(ny, nx, tEnd-tStart+1);
v = zeros(ny, nx, tEnd-tStart+1);

% dir_fig = ['Multi-' method '_co' num2str(c0) '_dout' num2str(d_out) '_h' num2str(h) '_cd' num2str(density_0) '/'];
dir_fig = ['C:\Users\HLiang\Desktop\Case study LP\' 'Multi-' method '_co' num2str(c0) '_dout' num2str(d_out) '_h' num2str(h) '_cd' num2str(density_0) '_pole' num2str(h_pole) '/'];
% dir_fig = ['C:\Users\HOWIE-PC\Desktop\Case study LP\' 'Multi-' method '_co' num2str(c0) '_dout' num2str(d_out) '_h' num2str(h) '_cd' num2str(density_0) '_pole' num2str(h_pole) '/' '/'];

boundary_judge_H = gfuns.Boundary_value(x,y,ones(ny,nx),boundary_H,0);
turbu_nan = nan(ny,nx,tEnd-tStart+1);
for i = (tStart+1):(tEnd+1)
    load([dir_fig num2str(i-1) '.mat']);
    Q_cell_com = {zeros(ny,nx),zeros(ny,nx),zeros(ny,nx)};
    Q_cell_com{1} = gfuns.Boundary_value(x,y,Q_cell_com{1},boundary_H,0);
    for j = 1:n_OD
        for k = 1:3
            Q_cell_com{k} = Q_cell_com{k} + Q_cell{j}{k};
        end
    end
    for iy = (nr+n_area(2,1)):(n_area(2,2)-nr)
        for ix = (nr+n_area(1,1)):(n_area(1,2)-nr)
            if all(boundary_judge_H((iy-nr):(iy+nr),(ix-nr):(ix+nr)) == 1)
                density_com(iy,ix,i) = mean(mean(Q_cell_com{1}((iy-nr):(iy+nr),(ix-nr):(ix+nr))));
                fx_com = mean(mean(Q_cell_com{2}((iy-nr):(iy+nr),(ix-nr):(ix+nr))));
                fy_com = mean(mean(Q_cell_com{3}((iy-nr):(iy+nr),(ix-nr):(ix+nr))));
                v(iy,ix,i) = sqrt(fx_com^2+fy_com^2)/max(10^-4,abs(density_com(iy,ix,i)));
                turbu_nan(iy,ix,i) = 1;
            end
        end
    end
    if mod(i,10)==0
        tend = toc(tstart);
        fprintf('Computing time is %f, t is %d\n',tend, i);
    end
end
for i = (pe+tStart+1):(tEnd-pe+1)
    v_ave = mean(v(:,:,(i-pe):(i+pe)),3);
    turbu(:,:,i) = turbu_nan(:,:,i).*density_com(:,:,i).*(v(:,:,i) - v_ave).^2;
%     v_ave = mean(v(n_area(2,1):n_area(2,2),n_area(1,1):n_area(1,2),(i-pe):(i+pe)),3);
%     turbu(n_area(2,1):n_area(2,2),n_area(1,1):n_area(1,2),i) = turbu_nan(n_area(2,1):n_area(2,2),n_area(1,1):n_area(1,2),i).*density_com(n_area(2,1):n_area(2,2),n_area(1,1):n_area(1,2),i).*(v(n_area(2,1):n_area(2,2),n_area(1,1):n_area(1,2),i) - v_ave).^2;
end
turbu_t = reshape(max(max(turbu,[],1),[],2),[1,tEnd-tStart+1]);
turbu_t_ave = turbu_t;


writerObj = VideoWriter([dir_fig 'video_P.avi']);
writerObj.Quality = 100;
writerObj.FrameRate = 5;
open(writerObj);
f_num = 1;
x1_i = 10/h_i+1;
x2_i = 95/h_i;
x_plot = 1/2*h_i:h_i:(85-1/2*h_i);
pressure_m_t = zeros(1,tEnd+1);
for t = 500:800 
    load([dir_fig num2str(t) '.mat']);
    D_in = cell(1,n_OD);
    for i_OD = 1:n_OD
        D_in{i_OD} = gfuns.D_in(t,i_OD);
    end
    tao_ex = tao_n.*(1-gfuns.Panic(t).*panic_ex)+tao_p.*gfuns.Panic(t).*panic_ex;
    [vep_x,vep_y] = RHS_FSM(Q_cell, tao_ex);
    pressure(pressure>=10^6) = 0;
    pressure_m_t(t+1) = max(max(pressure));
    fprintf([method '.t = %d. \n'],t)
    pressure = gfuns.Boundary_value(x,y,pressure,boundary_H,nan);
    fig = figure(1);
    pressure_interp = interp2(X,Y,pressure,X_i,Y_i);
    imagesc(x_plot,y_i, pressure_interp(:,x1_i:x2_i),'alphadata',~isnan(pressure_interp(:,x1_i:x2_i))); colormap Jet; set(gca,'YDir','normal','color',0*[1 1 1]); caxis([0 1000]);
    colorbar;
    axis([40,65,15,40]); xlabel('x (m)'); ylabel('y (m)');
    set(fig,'unit','centimeters','position',[10 10 10 8]);
    title(['Time(s): ',num2str(t)])
    F=getframe(gcf);
    writeVideo(writerObj,F);
end
close(writerObj);

fig = figure();
t_p = 400:800;
axes('Units','normalized','Position',[0.1 0.15 0.78 0.8]);
yyaxis left;
plot(t_p,turbu_t_ave(t_p+1),'b');
xlabel('Time (s)');ylabel('TP');ylim([0 0.06]);set(gca,'YColor','k');
hold on;
yyaxis right;
plot(480:800,pressure_m_t(481:801),'r','LineStyle','--');
xlim([480 800]);ylabel('Max. Pressure (N/m)');ylim([0 1400]);set(gca,'YColor','k');
hold off;
legend('TP','P_2')
set(fig,'Units','centimeters','Position',[10 10 15 8]);



dx = 0.1;
x_ave= 0:dx:10;density_ave= 1/2*dx:dx:(10-1/2*dx);
turbu_ave= zeros(1,length(x_ave)-1);
count= zeros(1,length(x_ave)-1);
density_s = reshape(density_com,[1,(tEnd-tStart+1)*ny*nx]);
turbu_s = reshape(turbu,[1,(tEnd-tStart+1)*ny*nx]);
density_s(isnan(turbu_s)) = [];
turbu_s(isnan(turbu_s)) = [];
for i = 2:length(x_ave)
    for j = 1:length(density_s)
        if density_s(j)>=x_ave(i-1) && density_s(j)<x_ave(i)
            turbu_ave(i-1) = turbu_s(j)+turbu_ave(i-1);
            count(i-1) = 1+count(i-1);
        end
    end
end
turbu_ave=turbu_ave./count;
fig2 = figure();
axes('Units','normalized','Position',[0.09 0.18 0.88 0.8]);

% plot(density_s,turbu_s,'k.','MarkerSize',1); hold on;
plot(density_ave,turbu_ave,'r-','MarkerSize',2);

xlim([0 10]);ylim([0 0.06]);xlabel('Density (ped/m^2)');ylabel('Turbulence Pressre')
set(fig2,'Units','centimeters','Position',[10 10 12 6]);


t_plot = 800;
turbu_plot = gfuns.Boundary_value(x,y,turbu(:,:,t_plot+1),boundary_H,nan);
fig = figure();
imagesc(x,y,turbu_plot,'alphadata',~isnan(turbu_plot)); colormap Jet; set(gca,'YDir','normal','color',0*[1 1 1]);caxis([0 0.06]); colorbar;
xlabel('x (m)'); ylabel('y (m)');
xlim([30 75]);ylim([10 50]);
set(fig,'unit','centimeters','position',[10 5 18 10]);

TP_max = max(turbu_t(t_p));
tend = toc(tstart);
fprintf('Computing time is %f, D_out is %f\n',tend, d_out);








