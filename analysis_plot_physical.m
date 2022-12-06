clc
clear all
global x y h nx ny
global boundary_O boundary_D boundary_D_2h boundary_O_2h boundary_D_outside Point_panic
global sigma epsilon density_m density_0 m_ave tao c0 beta method panic_config
global gfuns
global tstep tStart

h = 0.5; tStart = 0; tEnd = 0; method = 'CWWENO3'; panic_config = '_no_panic'; plot_config = 1;

m_ave = 60; density_0 = 6; density_m = 10; tao = 0.61; c0 = 1.2;
area = [0,100;0,50]; CFL = 0.2; tstep = 0.01; d_min = -10;
sigma = 10^-9; epsilon = 10^-6; gfuns = functions_given;

tstart = tic;
% Computational domain
x = area(1,1):h:area(1,2); nx = length(x);
y = area(2,2):-h:area(2,1);  ny = length(y);
[X,Y] = meshgrid(x,y);

% boundaries
boundary_O = {'Rectangle',1,[-5 0+0.5*h;-0.5*h 50-0.5*h]};
boundary_O_2h = {'Rectangle',1,[-5 0+0.5*h;0+2*h 50-0.5*h]};

boundary_D = {'Rectangle',1,[100 0.5*h;105 50-0.5*h]};
boundary_D_2h = {'Rectangle',1,[100-h -0.5*h;105 50+0.5*h]};
boundary_D_outside = {'Rectangle',1,[100+h 0.5*h;105 50-0.5*h]};
boundary_H_release = {'Rectangle',5,[-5 -5;105 0.5*h],[-5 50-0.5*h;105 55],[60-0.5*h 0.5*h;65+0.5*h 10+0.5*h],[60-0.5*h 40-0.5*h;65+0.5*h 50-0.5*h],[60-0.5*h 20-0.5*h;65+0.5*h 30+0.5*h]};
boundary_H_bound = {'Rectangle',5,[-5 -5;105 0.5*h],[-5 50-0.5*h;105 55],[60-0.5*h 0.5*h;65+0.5*h 10+0.5*h],[60-0.5*h 40-0.5*h;65+0.5*h 50-0.5*h],[60-0.5*h 20-0.5*h;65+0.5*h 33+0.5*h]};
Point_panic = [60 31.5];

    %% Plot&Save
fig1 = figure(1);
axes('Units','normalized','Position',[0 0.16 1.0 0.8]);
rectangle('Position',[area(1,1),area(2,1),area(1,2)-area(1,1),area(2,2)-area(2,1)],'LineWidth',1,'LineStyle','-'); hold on;


r = 20; theta = 0:pi/100:2*pi;
x = Point_panic(1) + r*cos(theta); y = Point_panic(2) + r*sin(theta);
hold on; axis equal;
% panic_area = fill(x,y,'r');
% set(panic_area,'edgealpha',0,'facealpha',0.3);

n = 1;
boundary = boundary_H_release;
while(n<=length(boundary))
    for k = (n+2):(n+1+boundary{n+1})
        rectangle('Position',[boundary{k}(1,1),boundary{k}(1,2),boundary{k}(2,1)-boundary{k}(1,1),boundary{k}(2,2)-boundary{k}(1,2)],'LineWidth',1,'LineStyle','-'); hold on;
        x = [boundary{k}(1,1) boundary{k}(2,1) boundary{k}(2,1) boundary{k}(1,1)];
        y = [boundary{k}(1,2) boundary{k}(1,2) boundary{k}(2,2) boundary{k}(2,2)];
        fill(x,y,'k');hold on;
    end
    n = n + boundary{n+1}+2;
end
% k = 7;
% x = [boundary{k}(1,1) boundary{k}(2,1) boundary{k}(2,1) boundary{k}(1,1)];
% y = [boundary{k}(2,2) boundary{k}(2,2) boundary{k}(2,2)+3 boundary{k}(2,2)+3];
% fill(x,y,'r'); hold on;

xlabel('X (m)'); ylabel('Y (m)'); xlim([area(1,1) area(1,2)]); ylim([area(2,1) area(2,2)]);
set(fig1,'unit','centimeters','position',[10 10 12 6]);
hold off