global density_m density_0 m_ave tao c0 beta method panic_config
tstart = tic;

h = 0.5; tStart = 0; tEnd = 810;

method = 'LFHE'; panic_config = '_panic'; plot_config = 1;
m_ave = 60; density_0 = 6; density_m = 10; tao = 0; tao_n = 3; tao_p = 0.2;
c0 = 0.6; beta = 0; density_a = 8; c1 = c0.*density_0.^beta./(density_0-density_a); b_c = 120; v_free = 1.034;
x_max = 105; y_max = 50;
area = [0,x_max;-2,y_max]; CFL = 0.2; tstep = 0.01; d_min = 0.01;
t_bound = [300 1200]; t_panic = [480, 600, 1200];
x_bound = [40,65]; y_bound = [y_max,y_max]; d_out = 6;
sigma = 10^-9; epsilon = 10^-6; gfuns = functions_given; pfuns = functions_plot;

Turbu_area = [61,64;28,31];
area_x = 16:21; area_y = 11:16; ne = 0;
speed_edges = 0:0.01:0.1;
angle_edges = -180:10:180;
n_area(1,1) = (Turbu_area(1,1) - area(1,1))./h + 1;
n_area(1,2) = (Turbu_area(1,2) - area(1,1))./h;
n_area(2,1) = (area(2,2) - Turbu_area(2,2))./h + 1;
n_area(2,2) = (area(2,2) - Turbu_area(2,1))./h;
pe = 0;
t_p = 0:800;

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
boundary_D{5} = {'Rectangle',1,[64 29;65 30]};
boundary_D_2h{5} = {'Rectangle',1,[64 29;65 30]};
boundary_D_outside{5} = {'Rectangle',1,[64 29;65 30]};
boundary_D{6} = {'Rectangle',1,[50 0;55 2]};
boundary_D_2h{6} = {'Rectangle',1,[50 0;55 2]};
boundary_D_outside{6} = {'Rectangle',1,[50 0;55 2]};
boundary_D{7} = {'Rectangle',1,[40 y_max;65 y_max+5]};
boundary_D_2h{7} = {'Rectangle',1,[40 y_max-h;65 y_max+5]};
boundary_D_outside{7} = {'Rectangle',1,[40 y_max+h;65 y_max+5]};

% Physical boundary
boundary_H = {'Rectangle',5,[-5 20;40 y_max+5],[65 20;x_max+5 y_max+5],[-5 -5;x_max+5 0],[64 29;65 30],[50 0;55 2]};
boundary_H_panic = {'Rectangle',3,[-5 20;40 y_max+5],[65 20;x_max+5 y_max+5],[-5 -5;x_max+5 0]};
boundary_Pole = {'Rectangle',1,[64 29;65 30]};
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
        dis_Pole = sqrt((x_ex(i) - 64.5).^2+(y_ex(j) - 29.5).^2)./15;
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
        dis_Pole = ((x(i) - 64.5)./15).^2+((y(j) - 29.5)./15).^2;
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

% density_com = zeros(n_area(2,2)-n_area(2,1)+1, n_area(1,2)-n_area(1,1)+1, tEnd-tStart+1);
% u_simu = zeros(n_area(2,2)-n_area(2,1)+1, n_area(1,2)-n_area(1,1)+1, tEnd-tStart+1);
% v_simu = zeros(n_area(2,2)-n_area(2,1)+1, n_area(1,2)-n_area(1,1)+1, tEnd-tStart+1);
% entropy_m_simu = zeros(1, tEnd-tStart+1);
% entropy_a_simu = zeros(1, tEnd-tStart+1);

% dir_fig = [ 'Multi-' method '_co' num2str(c0) '_dout' num2str(d_out) '_h' num2str(h) '_cd' num2str(density_0) '/'];
% dir_fig = ['C:\Users\HLiang\Desktop\Case study LP\' 'Multi-' method '_co' num2str(c0) '_dout' num2str(d_out) '_h' num2str(h) '_cd' num2str(density_0) '/'];
% dir_fig = ['C:\Users\HOWIE-PC\Desktop\Case study LP\' 'Multi-' method '_co' num2str(c0) '_dout' num2str(d_out) '_h' num2str(h) '_cd' num2str(density_0) '/'];
% 
% boundary_judge_H = gfuns.Boundary_value(x,y,ones(ny,nx),boundary_H,0);
% for i = (tStart+1):(tEnd+1)
%     load([dir_fig num2str(i-1) '.mat']);
%     Q_cell_com = {zeros(ny,nx),zeros(ny,nx),zeros(ny,nx)};
%     Q_cell_com{1} = gfuns.Boundary_value(x,y,Q_cell_com{1},boundary_H,0);
%     for j = 1:n_OD
%         for k = 1:3
%             Q_cell_com{k} = Q_cell_com{k} + Q_cell{j}{k};
%         end
%     end
%     density_com(:,:,i) = Q_cell_com{1}(n_area(2,1):n_area(2,2),n_area(1,1):n_area(1,2));
%     fx_com = Q_cell_com{2}(n_area(2,1):n_area(2,2),n_area(1,1):n_area(1,2));
%     fy_com = Q_cell_com{3}(n_area(2,1):n_area(2,2),n_area(1,1):n_area(1,2));
%     u_simu(:,:,i) = fx_com./max(10^-4,abs(density_com(:,:,i)));
%     v_simu(:,:,i) = fy_com./max(10^-4,abs(density_com(:,:,i)));
%     if mod(i,100)==0
%         tend = toc(tstart);
%         fprintf('Computing time is %f, t is %d\n',tend, i);
%     end
%     % entropy
%     speed = reshape(sqrt(u_simu(:,:,i).^2+v_simu(:,:,i).^2),[1,(n_area(2,2)-n_area(2,1)+1).*(n_area(1,2)-n_area(1,1)+1)]);
%     angle = atan2d(reshape(v_simu(:,:,i),[1,(n_area(2,2)-n_area(2,1)+1).*(n_area(1,2)-n_area(1,1)+1)]),reshape(u_simu(:,:,i),[1,(n_area(2,2)-n_area(2,1)+1).*(n_area(1,2)-n_area(1,1)+1)]));
%     p_speed_dis = histcounts(speed,speed_edges)./(length(area_y).*length(area_x));
%     p_angle_dis = histcounts(angle,angle_edges)./(length(area_y).*length(area_x));
%     entropy_m_simu(i) = -sum(p_speed_dis(p_speed_dis>0).*log2(p_speed_dis(p_speed_dis>0)));
%     entropy_a_simu(i) = -sum(p_angle_dis(p_angle_dis>0).*log2(p_angle_dis(p_angle_dis>0)));
%     
% end
% speed = sqrt(u_simu.^2+v_simu.^2);
% speed_ave = mean(mean(speed,1),2);
% u_simu_ave = mean(mean(u_simu,1),2);
% v_simu_ave = mean(mean(v_simu,1),2);
% speed_ave = reshape(speed_ave,[1,tEnd-tStart+1]);
% u_simu_ave = reshape(u_simu_ave,[1,tEnd-tStart+1]);
% v_simu_ave = reshape(v_simu_ave,[1,tEnd-tStart+1]);
%% Groups
% 35s   16:30:00-16:30:35 

%% Filter

% filename = 'Videos\results_frames.mat';
filename = 'Videos\results_frames_35test.mat';
load(filename);

%% calculate t sequence
t_seq_point = [0 35;0 0];
tstep = 0.2;
t_seq = 0.1:0.2:(600-0.1);
u_video = zeros(600/tstep,1);
v_video = zeros(600/tstep,1);
entropy_m_video = zeros(600/tstep,1);
entropy_a_video = zeros(600/tstep,1);
n = 24;
for i = 1:(600/tstep)
    t = tstep*i-0.1;
    for j = 1:1
        if t>=(t_seq_point(j,1)+5) && t<=(t_seq_point(j,2)-5)
            k = j-1+sum(t_seq_point(1:(j-1),2)-t_seq_point(1:(j-1),1))*5+round((t-t_seq_point(j,1)+0.1)/0.2);
            u_filtered_sum = zeros(length(area_y),length(area_x));
            v_filtered_sum = zeros(length(area_y),length(area_x));
            for l = k-n:k+n
                u_filtered_sum = u_filtered_sum+u_filtered{l}(area_y,area_x);
                v_filtered_sum = v_filtered_sum+v_filtered{l}(area_y,area_x);
            end
            u_filtered_1 = u_filtered_sum ./ (2*n+1);
            v_filtered_1 = v_filtered_sum ./ (2*n+1);
            u_video(i) = mean(mean(u_filtered_1));
            v_video(i) = mean(mean(v_filtered_1));
            
            % entropy
            speed = reshape(sqrt(u_filtered_1.^2+v_filtered_1.^2),[1,length(area_y).*length(area_x)]);
            angle = atan2d(reshape(v_filtered_1,[1,length(area_y).*length(area_x)]),reshape(u_filtered_1,[1,length(area_y).*length(area_x)]));
            p_speed_dis = histcounts(speed,speed_edges)./(length(area_y).*length(area_x));
            p_angle_dis = histcounts(angle,angle_edges)./(length(area_y).*length(area_x));
            entropy_m_video(i) = -sum(p_speed_dis(p_speed_dis>0).*log2(p_speed_dis(p_speed_dis>0)));
            entropy_a_video(i) = -sum(p_angle_dis(p_angle_dis>0).*log2(p_angle_dis(p_angle_dis>0)));
            break;
        else
            u_video(i) = nan;
            v_video(i) = nan;
            entropy_m_video(i) = nan;
            entropy_a_video(i) = nan;
        end
    end
end
%% Derive speed
% entropy
fig6 = figure(6);
% plot(t_seq+350,entropy_a_video.*entropy_m_video,'b'); hold on;
% plot(t_p,entropy_a_simu(t_p+1).*entropy_m_simu(t_p+1),'r--'); hold off;
plot(t_seq+350,entropy_a_video,'b'); hold on;
plot(t_p,entropy_a_simu(t_p+1),'r--'); hold off;
% plot(t_seq+350,entropy_m_video,'b'); hold on;
% plot(t_p,entropy_m_simu(t_p+1),'r--'); hold off;
xlim([0 t_p(end)]);
xlabel('Time (s)');ylabel('Entropy');legend('Video','Simulation');
set(fig6,'Units','centimeters','Position',[10 10 15 8]);

