%% Groups
% 35s   16:30:00-16:30:35 
% 61s
% 44s   16:31:36-16:32:20 
% 112s
% 98s   16:34:12-16:35:50 
% 30s
% 58s   16:36:20-16:37:18 
% 14s
% 53s   16:37:32-16:38:25 
% 95s
%%
clc
clear
global v
filename = 'Videos\Kamera13_1620_1640.mp4';
% dir_fig = 'C:\Users\HLiang\Desktop\Case study LP\VideoFrames\';
dir_fig = 'C:\Users\HOWIE-PC\\Desktop\Case study LP\VideoFrames\';

mkdir(dir_fig);

v = VideoReader(filename);
FR = v.FrameRate;
t_seq = [0 96 252 380 452];
n_Fstart = [15001 17214 21114 24314 26114];
n_Fend = [15876 18314 23564 25764 27439];

%G0: 0,0; 600,0; 600,0; 600,600
%G1: 446,434; 820,434; 621,617; 1151;617
%G2: 

dot_crop = process_points();




%% Group 1
for ng = 1:5
    for i = n_Fstart(ng):FR/5:n_Fend(ng)
        t_frame = t_seq(ng)+(i-n_Fstart(ng))/FR/1.0;
        process_frame(dot_crop{ng},t_frame,i);
        if mod(t_frame,10)==0
            fprintf('Process frame: %d\n',i);
        end
    end
end


%% Select frame
% f_start_b = read(v,n_Fstart(5)-1);
% f_end_b = read(v,n_Fend(5)-1);
% f_start = read(v,n_Fstart(5));
% f_end = read(v,n_Fend(5));
% fig1 = figure();
% subplot(2,2,1);imshow(imcrop(f_start_b,[1100,660,180,60]));
% subplot(2,2,2);imshow(imcrop(f_start,[1100,660,180,60]));
% subplot(2,2,3);imshow(imcrop(f_end_b,[1100,660,180,60]));
% subplot(2,2,4);imshow(imcrop(f_end,[1100,660,180,60]));

%% Solve transformation
%G0: 0,0; 600,0; 600,0; 600,600
%G1: 446,434; 820,434; 621,617; 1151;617
% syms t11 t12 t13 t21 t22 t23 t31 t32 t33
% x1 = 446; y1 = 434;
% x2 = 820; y2 = 434;
% x3 = 621; y3 = 617;
% x4 = 1151; y4 = 617;
% eq1 = t11*x1+t12*y1+t13*1==0;
% eq2 = t21*x1+t22*y1+t23*1==0;
% eq3 = t31*x1+t32*y1+t33*1==600;
% eq4 = t11*x2+t12*y2+t13*1==0;
% eq5 = t21*x2+t22*y2+t23*1==0;
% eq6 = t31*x2+t32*y2+t33*1==1;
% eq7 = t11*x3+t12*y3+t13*1==0;
% eq8 = t21*x3+t22*y3+t23*1==-1;
% eq9 = t31*x3+t32*y3+t33*1==1;
% sol = solve(eq1,eq2,eq3,eq4,eq5,eq6,eq7,eq8,eq9);












