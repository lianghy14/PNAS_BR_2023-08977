clc
clear

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
tStart      =   13800; % 60,120,180,...
tEnd        =   13800;
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

%% LTY Slope
global wm1 wm2
slope   =   zeros(ny_ex3,nx_ex3);
judge   =   zeros(ny_ex3,nx_ex3,n_OD);
for i = 1:nx_ex3
    for j = 1:ny_ex3
        if x_ex3(i)>=40 && x_ex3(i)<=44 && y_ex3(j)>=wm1/1.5 && y_ex3(j)<=(area_cal(2,2)-wm2/1.5)
            slope(j,i) = 0.10;
        end
%         if bj_ex3ODH(j,i,1)~=0
%             dis = ((y_ex3(j)-1/2*(area_cal(2,2)+wm1/2-wm2/2))^2+(x_ex3(i)-42)^2)/((1/2*(area_cal(2,2)-wm1-wm2-4))^2);
%             judge(j,i,:) = ones(1,6).* max( 0,1 - dis);
%         end
        judge(j,i,3:4) = ones(1,2);
    end
end


%% IC for density potential

fprintf('Start Ploting---------------------------------------------\n');
color_map_1=[linspace(0,1,100)',linspace(1,0,100)',linspace(0,0,100)'];
writerObj = VideoWriter([dir_data '00simulation.avi']);
writerObj.Quality = 100;
writerObj.FrameRate = 5;
open(writerObj);
f_num = 1;
x1 = 0/h+1;
x2 = area_cal(1,2)/h;
x_plot = 1/2*h:h:(area_cal(1,2)-1/2*h);
for t = tStart:300:tEnd
    load ([dir_data num2str(Record_dt) '_' num2str(t/Record_dt) '.mat']);
    density_com_1 = zeros(ny,nx);
    density_com_2 = zeros(ny,nx);
    for i = 1:n_OD
        switch i
            case {4,5,6}
                density_com_1 = density_com_1 + Record_Q(:,:,1,i,Record_dt);
            case {1,2,3}
                density_com_2 = density_com_2 + Record_Q(:,:,1,i,Record_dt);
        end
    end
    density_com_1 = gfuns.Boundary_value(x,y,density_com_1,boundary_H{1},nan);
    density_com_2 = gfuns.Boundary_value(x,y,density_com_2,boundary_H{1},nan);
    density_bi = (density_com_2 - density_com_1)./(density_com_1 + density_com_2);


    fig = figure(); imagesc(x_plot,y,density_bi(:,x1:x2),'alphadata',~isnan(density_bi(:,x1:x2)));
    colormap(color_map_1); set(gca,'YDir','normal','color',0*[1 1 1]);caxis([-1 1]);
    fig_colorbar = colorbar;
    axis([30,55,0,area_cal(2,2)]);
%         title(['Time(s): ',num2str(t)])
    xlabel('x (m)');
    ylabel('y (m)');
    set(fig,'unit','centimeters','position',[10 10 4 5]);
    set(gca, 'LooseInset', [0,0,0,0]);   
    print(fig, '-depsc',['C:\Users\Administrator\OneDrive - The University of Hong Kong\03_PhD_1\05 Paper04\manuscript\figures\bidirec' num2str(t)]);
    
    F=getframe(gcf);
    [I,map]=rgb2ind(frame2im(getframe(gcf)),256);
    if f_num == 1
        imwrite(I,map,[dir_data 'GIF.gif'],'gif', 'Loopcount',inf,'DelayTime',0.1);
    else
        imwrite(I,map,[dir_data 'GIF.gif'],'gif','WriteMode','append','DelayTime',0.1);
    end
    f_num=f_num+1;
    
    writeVideo(writerObj,F);
end
close(writerObj);
