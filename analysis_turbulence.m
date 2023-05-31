clc
clear all
close all
%% Run fortran codes
if exist('F90_FDM_God.mexw64','file')==0
    error('F90_FDM.mexw64 does not exist.')
end
% mex F90_FDM_God.f90;

%% Parameters and functions
gfuns   =   functions_given;
pfuns   =   functions_plot;
% pfuns.potential(P2_state,x_ex3,y_ex3,boundary_H{1})
%% Parameters
global area_cal h t_panic n_OD boundary_O boundary_D boundary_D_2h boundary_H
global dir_data Record_dt tEnd
global Q_state wm1 wm2 Record_OD xdet ydet

gfuns.Para('LTY20230413');

%% Layout settings
% 0 is bound_H, 1 is bound_D, 2 is bound_O, 3 is 2h near bound_D
x               =   (area_cal(1,1)+0.5*h):h:(area_cal(1,2)-0.5*h);      nx = length(x);
y               =   (area_cal(2,2)-0.5*h):-h:(area_cal(2,1)+0.5*h);     ny = length(y);
x_ex3           =   (x(1)-3*h):h:(x(end)+3*h);              nx_ex3 = length(x_ex3);
y_ex3           =   (y(1)+3*h):-h:(y(end)-3*h);             ny_ex3 = length(y_ex3);
[bj_ex3ODH]     =   gfuns.Layout(n_OD,boundary_O,boundary_D,boundary_D_2h,boundary_H,x_ex3,y_ex3);
[bj_calODH]     =   gfuns.Layout(n_OD,boundary_O,boundary_D,boundary_D_2h,boundary_H,x,y);
[X,Y]           =   meshgrid(x,y);

xdet = zeros(6,1); ydet = zeros(6,1);
xdet(1) =  (20 - x_ex3(1))/h;                   ydet(1) =  (y_ex3(1) - wm1/2)/h;
xdet(2) =  (area_cal(1,2) - 20 - x_ex3(1))/h;   ydet(2) =  (y_ex3(1) - wm1/2)/h;
xdet(3) =  (area_cal(1,2) - 30 - x_ex3(1))/h;   ydet(3) =  (y_ex3(1) - wm1/2)/h;
xdet(4) =  (30 - x_ex3(1))/h;                   ydet(4) =  (y_ex3(1) - area_cal(2,2) + wm2/2)/h;
xdet(5) =  (20 - x_ex3(1))/h;                   ydet(5) =  (y_ex3(1) - area_cal(2,2) + wm2/2)/h;
xdet(6) =  (area_cal(1,2) - 20 - x_ex3(1))/h;   ydet(6) =  (y_ex3(1) - area_cal(2,2) + wm2/2)/h;
xdet = round(xdet);
ydet = round(ydet);

%% Simulation period
tStart      =   0; % 60,120,180,...
tEnd        =   15000;
gfuns       =   functions_given;
pfuns       =   functions_plot;

%% Record data
Record_Q    =   zeros(ny,nx,3,n_OD,Record_dt);
Record_Ve   =   zeros(ny,nx,2,n_OD,Record_dt);
Record_P2   =   zeros(ny,nx,Record_dt);
Record_F2   =   zeros(ny,nx,2,Record_dt);
Record_OD   =   zeros(tEnd/10,n_OD);

%% Initialization
tstep       =   0.01;
CFL         =   0.2;
cpu_tstart  =   tic;
t           =   tStart;
Q_state     =   zeros(ny_ex3,nx_ex3,3,n_OD);
Ve_state    =   zeros(ny_ex3,nx_ex3,2,n_OD);
panic       =   zeros(ny_ex3,nx_ex3,n_OD);
bj_out      =   zeros(2,1);
P2_state    =   zeros(ny_ex3,nx_ex3);
F2_state    =   zeros(ny_ex3,nx_ex3,2);
alpha_LF    =   5;
rn          =   60;
mkdir(dir_data);
if t > 0
    load ([dir_data num2str(Record_dt) '_' num2str(t/Record_dt) '.mat']);
    Q_state(4:end-3,4:end-3,:,:) = Record_Q(:,:,:,:,Record_dt);
end


t_p = 0:10:240;
Turbu_area = [40,44;wm1,area_cal(2,2)-wm2]; pe = 2;
speed_edges = 0:0.01:0.1;
angle_edges = -180:10:180;


%% LTY Slope
global wm1 wm2
slope   =   zeros(ny_ex3,nx_ex3);
judge   =   zeros(ny_ex3,nx_ex3,n_OD);
for i = 1:nx_ex3
    for j = 1:ny_ex3
        if x_ex3(i)>=40 && x_ex3(i)<=44 && y_ex3(j)>=(wm1) && y_ex3(j)<=(area_cal(2,2)-wm2)
            slope(j,i) = 0.15;
        end
%         if bj_ex3ODH(j,i,1)~=0
%             dis = ((y_ex3(j)-1/2*(area_cal(2,2)+wm1/2-wm2/2))^2+(x_ex3(i)-42)^2)/((1/2*(area_cal(2,2)-wm1-wm2-4))^2);
%             judge(j,i,:) = ones(1,6).* max( 0,1 - dis);
%         end
        judge(j,i,3:4) = ones(1,2);
    end
end


%% Initialization
tstart = tic;
n_area(1,1) = (Turbu_area(1,1) - area_cal(1,1))./h + 1;
n_area(1,2) = (Turbu_area(1,2) - area_cal(1,1))./h;
n_area(2,1) = (area_cal(2,2) - Turbu_area(2,2))./h + 1;
n_area(2,2) = (area_cal(2,2) - Turbu_area(2,1))./h;
if exist([dir_data 'turbu.mat'])==0
    Turbuloc_F_fric = zeros(ny, nx, (tEnd-tStart)/Record_dt+1);
    Turbuloc_F_fgrad = zeros(ny, nx, (tEnd-tStart)/Record_dt+1);
    Turbuloc_F_fp2 = zeros(ny, nx, (tEnd-tStart)/Record_dt+1);
    Turbuloc_P = zeros(ny, nx, (tEnd-tStart)/Record_dt+1);
    Turbuloc_Den = zeros(ny, nx, (tEnd-tStart)/Record_dt+1);
    Turbuloc_Vx = zeros(ny, nx, (tEnd-tStart)/Record_dt+1);
    Turbuloc_Vy = zeros(ny, nx, (tEnd-tStart)/Record_dt+1);
    t_min = 0;
    for t = (tStart):Record_dt:(tEnd)
        t_min = t_min + 1;
        load ([dir_data num2str(Record_dt) '_' num2str(t/Record_dt) '.mat']);
        Turbuloc_P(:,:,t_min) = Record_P2(:,:,Record_dt);
        Turbuloc_Den(:,:,t_min) = sum(Record_Q(:,:,1,:,Record_dt),4);
        Turbuloc_Vx(:,:,t_min) = sum(Record_Q(:,:,2,:,Record_dt),4)./sum(Record_Q(:,:,1,:,Record_dt),4);
        Turbuloc_Vy(:,:,t_min) = sum(Record_Q(:,:,3,:,Record_dt),4)./sum(Record_Q(:,:,1,:,Record_dt),4);

        panic           =   judge(4:end-3,4:end-3,:) .* gfuns.Panic(t);
        tao             =   5.*(1-panic) + 0.5.*panic;
        for i_OD = 1:n_OD
            f_x = Record_Q(:,:,1,i_OD,Record_dt).*Record_Ve(:,:,1,i_OD,Record_dt)-Record_Q(:,:,2,i_OD,Record_dt);
            f_y = Record_Q(:,:,1,i_OD,Record_dt).*Record_Ve(:,:,2,i_OD,Record_dt)-Record_Q(:,:,3,i_OD,Record_dt);
            Turbuloc_F_fric(:,:,t_min) = Turbuloc_F_fric(:,:,t_min) + 65.*sqrt(f_x.^2+f_y.^2)./tao(:,:,i_OD);
        end

        alpha       =   min(1,max(0, (Turbuloc_Den(:,:,t_min) - 6) ./ (10 - 6)));
        Turbuloc_F_fgrad(:,:,t_min) = 65.*Turbuloc_Den(:,:,t_min).*slope(4:end-3,4:end-3).*alpha;


        Turbuloc_F_fp2(:,:,t_min) = sqrt(Record_F2(:,:,1,Record_dt).^2+Record_F2(:,:,2,Record_dt).^2);

        Turbuloc_F_fric(:,:,t_min) = Turbuloc_F_fric(:,:,t_min)./(max(10^-6,Turbuloc_Den(:,:,t_min)).*(h.^2));
        Turbuloc_F_fgrad(:,:,t_min) = Turbuloc_F_fgrad(:,:,t_min)./(max(10^-6,Turbuloc_Den(:,:,t_min)).*(h.^2));
        Turbuloc_F_fp2(:,:,t_min) = Turbuloc_F_fp2(:,:,t_min)./(max(10^-6,Turbuloc_Den(:,:,t_min)).*(h.^2));

        if mod(t_min,10)==0
            tend = toc(tstart);
            fprintf('Computing time is %f, t is %d\n',tend, t);
        end
    end
    save([dir_data 'turbu.mat'],'Turbuloc*');
else
    load([dir_data 'turbu.mat']);
end
%% Animation
t_min = tEnd/Record_dt+1;
pressure_t = zeros(3,t_min);
density_t = zeros(3,t_min);
F_fric_t = zeros(3,t_min);
F_fgrad_t = zeros(3,t_min);
F_fp2_t = zeros(3,t_min);
VE_t = zeros(1,t_min);

for t = 1:t_min
    pressure_t(1,t) = mean(mean(Turbuloc_P(n_area(2,1):n_area(2,2),n_area(1,1):n_area(1,2),t)));
    pressure_t(2,t) = max(max(Turbuloc_P(n_area(2,1):n_area(2,2),n_area(1,1):n_area(1,2),t)));
    pressure_t(3,t) = min(min(Turbuloc_P(n_area(2,1):n_area(2,2),n_area(1,1):n_area(1,2),t)));
    
    density_t(1,t) = mean(mean(Turbuloc_Den(n_area(2,1):n_area(2,2),n_area(1,1):n_area(1,2),t)));
    density_t(2,t) = max(max(Turbuloc_Den(n_area(2,1):n_area(2,2),n_area(1,1):n_area(1,2),t)));
    density_t(3,t) = min(min(Turbuloc_Den(n_area(2,1):n_area(2,2),n_area(1,1):n_area(1,2),t)));
    
    F_fric_t(1,t) = mean(mean(Turbuloc_F_fric(n_area(2,1):n_area(2,2),n_area(1,1):n_area(1,2),t)));
    F_fric_t(2,t) = max(max(Turbuloc_F_fric(n_area(2,1):n_area(2,2),n_area(1,1):n_area(1,2),t)));
    F_fric_t(3,t) = min(min(Turbuloc_F_fric(n_area(2,1):n_area(2,2),n_area(1,1):n_area(1,2),t)));
    
    F_fgrad_t(1,t) = mean(mean(Turbuloc_F_fgrad(n_area(2,1):n_area(2,2),n_area(1,1):n_area(1,2),t)));
    F_fgrad_t(2,t) = max(max(Turbuloc_F_fgrad(n_area(2,1):n_area(2,2),n_area(1,1):n_area(1,2),t)));
    F_fgrad_t(3,t) = min(min(Turbuloc_F_fgrad(n_area(2,1):n_area(2,2),n_area(1,1):n_area(1,2),t)));
    
    F_fp2_t(1,t) = mean(mean(Turbuloc_F_fp2(n_area(2,1):n_area(2,2),n_area(1,1):n_area(1,2),t)));
    F_fp2_t(2,t) = max(max(Turbuloc_F_fp2(n_area(2,1):n_area(2,2),n_area(1,1):n_area(1,2),t)));
    F_fp2_t(3,t) = min(min(Turbuloc_F_fp2(n_area(2,1):n_area(2,2),n_area(1,1):n_area(1,2),t)));

    speed = reshape(sqrt(Turbuloc_Vx(n_area(2,1):n_area(2,2),n_area(1,1):n_area(1,2),t).^2+Turbuloc_Vy(n_area(2,1):n_area(2,2),n_area(1,1):n_area(1,2),t).^2),[1,(n_area(2,2)-n_area(2,1)+1).*(n_area(1,2)-n_area(1,1)+1)]);
    angle = atan2d(reshape(Turbuloc_Vy(n_area(2,1):n_area(2,2),n_area(1,1):n_area(1,2),t),[1,(n_area(2,2)-n_area(2,1)+1).*(n_area(1,2)-n_area(1,1)+1)]),reshape(Turbuloc_Vx(n_area(2,1):n_area(2,2),n_area(1,1):n_area(1,2),t),[1,(n_area(2,2)-n_area(2,1)+1).*(n_area(1,2)-n_area(1,1)+1)]));
    p_speed_dis = histcounts(speed,speed_edges)./((n_area(1,2)-n_area(1,1)+1).*(n_area(2,2)-n_area(2,1)+1));
    p_angle_dis = histcounts(angle,angle_edges)./((n_area(1,2)-n_area(1,1)+1).*(n_area(2,2)-n_area(2,1)+1));
    entropy_m_simu = -sum(p_speed_dis(p_speed_dis>0).*log2(p_speed_dis(p_speed_dis>0)));
    entropy_a_simu = -sum(p_angle_dis(p_angle_dis>0).*log2(p_angle_dis(p_angle_dis>0)));
    VE_t(t) = entropy_m_simu.*entropy_a_simu;
end

%% Evolution

fig = figure();

a_fric = area([t_p;t_p]',[F_fric_t(3,t_p+1);F_fric_t(2,t_p+1)-F_fric_t(3,t_p+1)]','FaceColor','blue','EdgeColor','none'); hold on;
set(a_fric(1),'FaceColor','none');
set(a_fric(2),'FaceAlpha',0.5);
plot(t_p,F_fric_t(1,t_p+1),'blue-'); hold on;

a_fgrad = area([t_p;t_p]',[F_fgrad_t(3,t_p+1);F_fgrad_t(2,t_p+1)-F_fgrad_t(3,t_p+1)]','FaceColor','yellow','EdgeColor','none'); hold on;
set(a_fgrad(1),'FaceColor','none');
set(a_fgrad(2),'FaceAlpha',0.5);
plot(t_p,F_fgrad_t(1,t_p+1),'yellow--'); hold on;

ypos = F_fp2_t(2,t_p+1) - F_fp2_t(1,t_p+1);
yneg = F_fp2_t(1,t_p+1) - F_fp2_t(3,t_p+1);

l_fp2 = errorbar(t_p,F_fp2_t(1,t_p+1),yneg,ypos,[],[],'-.or'); hold off;

legend([a_fric(2),a_fgrad(2),l_fp2],'$F_1$ (Friction force)','$F_2$ (Leaning force)','$F_3$ (Pushing force)','Location','northwest','Interpreter','latex');
ylabel('Force (N)');
xlabel('Time (min)');
ylim([0 250]); xlim([min(t_p)-10 max(t_p)+10]);
set(fig,'Units','centimeters','Position',[10 10 15 7]);
set(gca, 'Position', [0.08,0.15,0.90,0.83]);
print(fig, '-depsc','C:\Users\Administrator\OneDrive - The University of Hong Kong\03_PhD_1\05 Paper04\manuscript\figures\4forces');

%% Evolution

plot_data = {density_t,pressure_t,VE_t};
save_name = {'1density','2pressure','3VE'};
y_label = {'Density (ped/m^2)','Pressure (N/m)','VE'};
y_lim = {[0 10],[0 2500],[0 15]};


fig = figure();
j = 1;
yyaxis left;
ypos = plot_data{j}(2,t_p+1) - plot_data{j}(1,t_p+1);
yneg = plot_data{j}(1,t_p+1) - plot_data{j}(3,t_p+1);
errorbar(t_p,plot_data{j}(1,t_p+1),yneg,ypos,[],[],'bo-'); hold on;
ylabel(y_label{j});
set(gca,'ycolor','b');
xlabel('Time (min)');
ylim(y_lim{j});xlim([min(t_p)-10 max(t_p)+10]);

j = 2;
yyaxis right;
ypos = plot_data{j}(2,t_p+1) - plot_data{j}(1,t_p+1);
yneg = plot_data{j}(1,t_p+1) - plot_data{j}(3,t_p+1);
errorbar(t_p,plot_data{j}(1,t_p+1),yneg,ypos,[],[],'ro-'); hold on;
set(fig,'Units','centimeters','Position',[10 10 15 7]);
set(gca, 'Position', [0.08,0.15,0.82,0.83]);
ylabel(y_label{j});
set(gca,'ycolor','r');
ylim(y_lim{j});
print(fig, '-depsc','C:\Users\Administrator\OneDrive - The University of Hong Kong\03_PhD_1\05 Paper04\manuscript\figures\densitypressure');

fig = figure();
j = 3;
plot(t_p,plot_data{j}(t_p+1),'o-'); hold on;
set(fig,'Units','centimeters','Position',[10 10 15 7]);
set(gca, 'Position', [0.08,0.15,0.90,0.83]);
ylabel(y_label{j});
xlabel('Time (min)');
ylim(y_lim{j});xlim([min(t_p)-10 max(t_p)+10]);
print(fig, '-depsc',['C:\Users\Administrator\OneDrive - The University of Hong Kong\03_PhD_1\05 Paper04\manuscript\figures\',save_name{j}]);








