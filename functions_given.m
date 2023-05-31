function gfuns = functions_given
gfuns.Para              = @Para_cal;
gfuns.Layout            = @Layout_cal;
gfuns.Q_in              = @Q_in_cal;
gfuns.Panic             = @Panic_cal;
gfuns.Velosity          = @Velosity_cal;
gfuns.Boundary_value    = @Boundary_value_cal;
end
%% Parameters
function Para_cal(day)
global area_cal h t_panic t_bound OD_bound OD_dir n_OD boundary_O boundary_D boundary_D_2h boundary_H
global dir_data Record_dt wm1 wm2 Q_max Q_lim

switch day
    case 'LTY20230413'
        h           =   0.5;
        Record_dt   =   60;
        x_max       =   84; 
        y_max       =   61;
        area_cal    =   [0,x_max;0,y_max]; 
        t_panic     =   3600.*[1, 2, 3, 4];
        t_bound     =   3600.*[0, 0.5, 1, 2, 3, 5];
            %    left  right  up   down  left  right
        Q_max = [0.40, 0.40, 0.40, 0.40, 0.40, 0.40];
        Q_lim = [  3,   1.5,   2,   2,  1.5,   3];
        wm1 = 15;   wm2 = 6;
        % f1:A to B
        % f2:B to A
        % f3:B to C
        % f4:C to B
        % f5:C to D
        % f6:D to C      left  right  up   down  left  right
        OD_dir      = [    1,    2,    2,    1,    1,    2];
        OD_bound    = [    0,    0,    0,    0,    0,    0;%0
                        0.30, 0.30, 0.30, 0.30, 0.20, 0.40;%0.5
                        0.05, 0.05, 0.05, 0.10, 0.12, 0.18;%1

                        0.05, 0.05, 0.05, 0.10, 0.12, 0.18;%1.2
                        0.05, 0.05, 0.05, 0.10, 0.12, 0.18;%2
                        0.05, 0.05, 0.05, 0.10, 0.12, 0.18;%3
                        0.05, 0.05, 0.05, 0.10, 0.12, 0.18];%4
        n_OD = 6;
        boundary_O = cell(1,n_OD);
        boundary_D = cell(1,n_OD);
        boundary_H = cell(1,n_OD);
        boundary_D_2h = cell(1,n_OD);
        boundary_O{1} = {'Rectangle',1,[-5 0;0 wm1]}; % f1:A to B
        boundary_O{2} = {'Rectangle',1,[x_max 0;x_max+5 wm1]};% f2:B to A
        boundary_O{3} = {'Rectangle',1,[x_max 0;x_max+5 wm1]};% f3:B to C
        boundary_O{4} = {'Rectangle',1,[-5 y_max-wm2;0 y_max]}; % f4:C to B
        boundary_O{5} = {'Rectangle',1,[-5 y_max-wm2;0 y_max]}; % f5:C to D
        boundary_O{6} = {'Rectangle',1,[x_max y_max-wm2;x_max+5 y_max]};% f6:D to C
        
        % Destination boundaries
        boundary_D{1} = {'Rectangle',1,[x_max 0;x_max+5 wm1]};
        boundary_D_2h{1} = {'Rectangle',1,[x_max-h 0;x_max+5 wm1]};
        
        boundary_D{2} = {'Rectangle',1,[-5 0;0 wm1]};
        boundary_D_2h{2} = {'Rectangle',1,[-5 0;-h wm1]};
        
        boundary_D{3} = {'Rectangle',1,[-5 y_max-wm2;0 y_max]};
        boundary_D_2h{3} = {'Rectangle',1,[-5 y_max-wm2;h y_max]};
        
        boundary_D{4} = {'Rectangle',1,[x_max 0;x_max+5 wm1]};
        boundary_D_2h{4} = {'Rectangle',1,[x_max-h 0;x_max+5 wm1]};
        
        boundary_D{5} = {'Rectangle',1,[x_max y_max-wm2;x_max+5 y_max]};
        boundary_D_2h{5} = {'Rectangle',1,[x_max-h y_max-wm2;x_max+5 y_max]};
        
        boundary_D{6} = {'Rectangle',1,[-5 y_max-wm2;0 y_max]};
        boundary_D_2h{6} = {'Rectangle',1,[-5 y_max-wm2;h y_max]};
        
        % Physical boundary
        for i = 1:6
            boundary_H{i} = {'Rectangle',4,[-5 -5;x_max+5 0],[-5 y_max;x_max+5 y_max+5],[-5 wm1;40 y_max-wm2],[44 wm1;x_max+5 y_max-wm2]};
        end
        dir_data = [getenv('UserProfile') '\Desktop\Case study LTY\LTY20230413\'];
end
end
%% Layout
function [bj_cal] = Layout_cal(n_OD,boundary_O,boundary_D,boundary_D_2h,boundary_H,x,y)
    bj_cal = zeros(length(y),length(x),n_OD);
    for i = 1:n_OD
        bj_cal(:,:,i) = ones(length(y),length(x)).*9;
        bj_cal(:,:,i) = Boundary_value_cal(x,y,bj_cal(:,:,i),boundary_D_2h{i},3);
        bj_cal(:,:,i) = Boundary_value_cal(x,y,bj_cal(:,:,i),boundary_O{i},2);
        bj_cal(:,:,i) = Boundary_value_cal(x,y,bj_cal(:,:,i),boundary_D{i},1);
        bj_cal(:,:,i) = Boundary_value_cal(x,y,bj_cal(:,:,i),boundary_H{i},0);
    end
end
%% Panic
function Panic = Panic_cal(t)
global t_panic

if t<=t_panic(1)
    Panic = 0;
else
    if t<=t_panic(2)
        Panic = 0.7.*(t-t_panic(1)) ./ (t_panic(2)-t_panic(1));
    else
        if t<=t_panic(3)
            Panic = 0.7;
        else
            if t<=t_panic(4)
                Panic = 0.7 + 0.3.*(t - t_panic(3)) ./ (t_panic(4)-t_panic(3));
            else
                Panic = 1;
            end
        end
    end
end
end
%% Q_in
function Q_in = Q_in_cal(t)
global t_panic n_OD t_bound OD_bound OD_dir Q_state tEnd
global Record_OD xdet ydet Q_max Q_lim
D_in = zeros(3,1);
Q_in = zeros(3,n_OD);

if mod(t,10) == 0
    for i_OD = 1:n_OD
        if Q_state(ydet(i_OD),xdet(i_OD),1,i_OD)<Q_lim(i_OD)
            OD_judge = 1;
        else
            OD_judge = -1;
        end
        Record_OD(t/10+2,i_OD)   =   min(Q_max(i_OD),max(0,Record_OD(t/10+1,i_OD) + OD_judge*0.05));
    end
end
n1 = floor(t/10)+1;
n2 = ceil(t/10)+1;

for i_OD = 1:n_OD
      D_in(1) = (Record_OD(n2,i_OD) - Record_OD(n1,i_OD)) ./ 10 .* mod(t,10) + Record_OD(n1,i_OD);

%     for i = 1:(length(t_bound)-1)
%         if (t>=t_bound(i)) && (t<=t_bound(i+1))
%             D_in(1) = (OD_bound(i+1,i_OD) - OD_bound(i,i_OD)) ./ (t_bound(i+1) - t_bound(i)) .* (t - t_bound(i)) + OD_bound(i,i_OD);
%         end
%     end
    
    switch OD_dir(i_OD)
        case 1
            D_in(2) = Velosity_cal(D_in(1),Panic_cal(t)) .* D_in(1);
            D_in(3) = 0;
        case 2
            D_in(2) = -Velosity_cal(D_in(1),Panic_cal(t)) .* D_in(1);
            D_in(3) = 0;
        case 3
            D_in(2) = 0;
            D_in(3) = Velosity_cal(D_in(1),Panic_cal(t)) .* D_in(1);
        case 4
            D_in(2) = 0;
            D_in(3) = -Velosity_cal(D_in(1),Panic_cal(t)) .* D_in(1);
    end
    Q_in(:,i_OD) = D_in;
end

end
%% Fundmental diagram
function Velosity = Velosity_cal(density,panic)
gama_1      = -0.075.*(1-panic) - 0.05.*panic;
alpha       = 2;
Velosity    = 1.034.*exp(gama_1.*max(0,density).^alpha);
end
%% Boundary value
function cell = Boundary_value_cal(x,y,cell,boundary,value)
n = 1;
while(n<=length(boundary))
    switch boundary{n}
        case 'Rectangle'
            for k = (n+2):(n+1+boundary{n+1})
                for i = length(cell(1,:)):-1:1
                    for j = 1:length(cell(:,1))
                        if (x(i)>=boundary{k}(1,1))&&(x(i)<=boundary{k}(2,1))&&(y(j)>=boundary{k}(1,2))&&(y(j)<=boundary{k}(2,2))
                                cell(j,i) = value;
                        end
                    end
                end
            end
            n = n + boundary{n+1}+2;
        case 'Circle'
            for k = (n+2):(n+1+boundary{n+1})
                for i = length(cell(1,:)):-1:1
                    for j = 1:length(cell(:,1))
                        if ((x(i)-boundary{k}(1))^2+(y(j)-boundary{k}(2))^2) <= boundary{k}(3)^2
                                cell(j,i) = value;
                        end
                    end
                end
            end
            n = n + boundary{n+1}+2;
    end
end
end
