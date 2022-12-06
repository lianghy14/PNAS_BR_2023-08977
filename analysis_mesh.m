h = 0.5; 
x_max = 85; y_max = 50;
area = [0,x_max;-2,y_max+2];gfuns = functions_given; pfuns = functions_plot;
boundary_H = {'Rectangle',5,[-5 20;40 y_max+5],[65 20;x_max+5 y_max+5],[-5 -5;x_max+5 0],[64 29;65 30],[50 0;55 2]};
boundary_Pole = {'Rectangle',1,[64 29;65 30]};
boundary_Container = {'Rectangle',1,[50 0;55 2]};

% Computational domain
x = (area(1,1)):h:(area(1,2)); nx = length(x);
y = (area(2,2)):-h:(area(2,1));  ny = length(y);
[X,Y] = meshgrid(x,y);

fig = figure();
for i = 1:length(x)
    plot([x(i),x(i)],[area(2,1),area(2,2)],'k'); hold on;
end
for i = 1:length(y)
    plot([area(1,1),area(1,2)],[y(i),y(i)],'k'); hold on;
end
% rectangle('Position',[50 10 20 20],'FaceColor','w');
rectangle('Position',[0 20 30 30],'FaceColor','w');
rectangle('Position',[55 20 40 30],'FaceColor','w');
rectangle('Position',[54 29 1 1],'FaceColor','w');
rectangle('Position',[40 0 5 2],'FaceColor','w');
hold off;
axis([0,x_max,0,y_max]); xlabel('x (m)'); ylabel('y (m)');
set(fig,'Units','centimeters','Position',[10 5 15 8]);
set (gca,'position',[0.08 0.15 0.90 0.8] );





