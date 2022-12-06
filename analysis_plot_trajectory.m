% clc
% clear all
h = 0.5 ; tStart =  480 ; tEnd = 900;

method = 'LFHE'; panic_config = '_panic'; plot_config = 1;
m_ave = 60; density_0 = 6; density_m = 10; tao = 0; tao_n = 3; tao_p = 0.2;
c0 = 0.6; beta = 0; density_a = 8; c1 = c0./2; b_c = 120; v_free = 1.034;
x_max = 105; y_max = 50;
area = [0,x_max;-2,y_max]; CFL = 0.2; tstep = 0.01; d_min = 0.01;
t_bound = [180 1200]; t_panic = [480, 600, 1200];
x_bound = [40,65]; y_bound = [y_max,y_max]; d_out = 6;
sigma = 10^-9; epsilon = 10^-6; gfuns = functions_given; pfuns = functions_plot;

track_loc = zeros(2,tEnd-tStart+2);
track_v = zeros(2,tEnd-tStart+1);
track_loc(1,1) = 63; track_loc(2,1) = 30;

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
x_n1 = sqrt(245-14*sqrt(70))/21;
x_n2 = sqrt(245+14*sqrt(70))/21;
k_n1 = (322+13.*sqrt(70))./900;
k_n2 = (322-13.*sqrt(70))./900;
x_g = [-x_n2; -x_n1; 0; x_n1; x_n2]./2 + 0.5;
k_g = [k_n2; k_n1; 128/225; k_n1; k_n2]./2;


% dir_fig = ['C:\Users\HLiang\Desktop\Case study LP\' 'Multi-' method '_co' num2str(c0) '_dout' num2str(d_out) '_h' num2str(h) '_cd' num2str(density_0) '_pole' num2str(h_pole) '/'];
% dir_fig = ['C:\Users\HLiang\Desktop\Case study LP\' 'Multi-' method '_co' num2str(c0) '_dout' num2str(d_out) '_h' num2str(h) '_cd' num2str(density_0) '/'];
dir_fig = ['C:\Users\HOWIE-PC\Desktop\Case study LP\' 'Multi-' method '_co' num2str(c0) '_dout' num2str(d_out) '_h' num2str(h) '_cd' num2str(density_0) '_pole' num2str(h_pole) '/'];

%% IC for density potential
fprintf('Start Ploting---------------------------------------------\n');
x1 = 10/h+1;
x2 = 95/h;
x_plot = 1/2*h:h:(85-1/2*h);

for t = tStart:tEnd
        Q_name = [dir_fig num2str(t) '.mat']; 
        load([dir_fig num2str(t) '.mat']);
        vx = Q_cell{5}{2}./Q_cell{5}{1};
        vy = Q_cell{5}{3}./Q_cell{5}{1};
        track_v(1,t-tStart+1) = interp2(x,y,vx,track_loc(1,t-tStart+1),track_loc(2,t-tStart+1));
        track_v(2,t-tStart+1) = interp2(x,y,vy,track_loc(1,t-tStart+1),track_loc(2,t-tStart+1));
        track_loc(1,t-tStart+2) = track_v(1,t-tStart+1).*1+track_loc(1,t-tStart+1);
        track_loc(2,t-tStart+2) = track_v(2,t-tStart+1).*1+track_loc(2,t-tStart+1);
end 
det_t = 2;

% R:t = 500s; G:t = 800s
fig = figure(1);
r = 255:-round(255/(tEnd-tStart+1).*det_t):0;
g = 0:round(255/(tEnd-tStart+1).*det_t):255;
b = 0;
for i = 1:det_t:(tEnd-tStart+1-det_t)
    plot(track_v(1,i:det_t:(i+det_t)),track_v(2,i:det_t:(i+det_t)),'-o','color',[r((i-1)/det_t+1) g((i-1)/det_t+1) 0]./255,'MarkerSize',1); hold on;
end
axis([-0.05,0.05,-0.05,0.05]); xlabel('u (m/s)'); ylabel('v (m/s)');
set(fig,'unit','centimeters','position',[10 10 10 9]);

fig = figure(2);
for i = 1:det_t:(tEnd-tStart+1-det_t)
    plot(track_loc(1,i:det_t:(i+det_t)),track_loc(2,i:det_t:(i+det_t)),'-o','color',[r((i-1)/det_t+1) g((i-1)/det_t+1) 0]./255,'MarkerSize',1); hold on;
end
rectangle('Position',[64,29,h_pole,h_pole],'EdgeColor','none','FaceColor','k'); 

axis([60,65,28,33]); xlabel('x (m)'); ylabel('y (m/s)');
set(fig,'unit','centimeters','position',[10 10 10 9]);

