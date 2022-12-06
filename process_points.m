function dot_crop = process_points()

dot_crop = cell(1,5);
%左上，左下，右下,右上
dot_crop{1} = [197,434; 268,617; 1151,617;820,434];
dot_crop{2} = [124,505; 174,683; 1100,683;773,505];
dot_crop{3} = [61,448; 98,555; 696,555;484,448];
dot_crop{4} = [63,445; 95,552; 695,552;486,445];
dot_crop{5} = [155,443; 180,552; 788,552;580,443];

%% plots
filename = 'Videos\Kamera13_1620_1640.mp4';
v = VideoReader(filename);
n_Fstart = [15001 17214 21114 24314 26114];
ng = 1;
i_frame = n_Fstart(ng);
f = read(v,i_frame);
fig1 = figure(1);
imshow(f); hold on;
y=[dot_crop{ng}(1,1),dot_crop{ng}(2,1),dot_crop{ng}(3,1),dot_crop{ng}(4,1)];
x=[dot_crop{ng}(1,2),dot_crop{ng}(2,2),dot_crop{ng}(3,2),dot_crop{ng}(4,2)];
plot_x=[x,x(1)];
plot_y=[y,y(1)];
for i = 1:11
    line_y = [y(1)+i/12*(y(4)-y(1)),y(2)+i/12*(y(3)-y(2))];
    line_x = [x(1)+i/12*(x(4)-x(1)),x(2)+i/12*(x(3)-x(2))];
    plot(line_y,line_x,'r','LineWidth',1);hold on;
end
for i = 1:11
    line_y = [y(1)+i/12*(y(2)-y(1)),y(4)+i/12*(y(3)-y(4))];
    line_x = [x(1)+i/12*(x(2)-x(1)),x(4)+i/12*(x(3)-x(4))];
    plot(line_y,line_x,'r','LineWidth',1);hold on;
end
plot(plot_y,plot_x,'r','LineWidth',1);hold off;
end