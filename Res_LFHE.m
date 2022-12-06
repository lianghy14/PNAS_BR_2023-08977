function [alpha,div] = Res_LFHE(Q_cell,density_com,boundary_O,boundary_D,D_in)
global area x y h nx ny gfuns boundary_H boundary_Pole boundary_Container
global t t_bound t_panic x_bound y_bound d_out n_OD i_OD
% tstart=tic;
div = {zeros(ny,nx),zeros(ny,nx),zeros(ny,nx)};
% Boundary expand
% Q_vector
boundary_judge_x = ones(ny,nx);
boundary_judge_y = ones(ny,nx);
Q_cell_xex = cell(1,3);Q_cell_xex_D = cell(1,3);
Q_cell_yex = cell(1,3);Q_cell_yex_D = cell(1,3);
n = 2;
boundary_judge_H = gfuns.Boundary_value(x,y,boundary_judge_x,boundary_H,0);
Q_cell_xex{1} = gfuns.Boundary_ex(Q_cell{1},boundary_judge_H,'x',n,'symmetric+');
Q_cell_yex{1} = gfuns.Boundary_ex(Q_cell{1},boundary_judge_H,'y',n,'symmetric+');
density_com_xex = gfuns.Boundary_ex(density_com,boundary_judge_H,'x',n,'symmetric+');
density_com_yex = gfuns.Boundary_ex(density_com,boundary_judge_H,'y',n,'symmetric+');
Q_cell_xex{2} = gfuns.Boundary_ex(Q_cell{2},boundary_judge_H,'x',n,'symmetric-');
Q_cell_yex{2} = gfuns.Boundary_ex(Q_cell{2},boundary_judge_H,'y',n,'symmetric+');
Q_cell_xex{3} = gfuns.Boundary_ex(Q_cell{3},boundary_judge_H,'x',n,'symmetric+');
Q_cell_yex{3} = gfuns.Boundary_ex(Q_cell{3},boundary_judge_H,'y',n,'symmetric-');
for i = 1:3
    Q_cell_xex_D{i} = gfuns.Boundary_ex(Q_cell{i},boundary_judge_x,'x',n,'nearest');
    Q_cell_yex_D{i} = gfuns.Boundary_ex(Q_cell{i},boundary_judge_y,'y',n,'nearest');
end
density_com_xex_D = gfuns.Boundary_ex(density_com,boundary_judge_x,'x',n,'nearest');
density_com_yex_D = gfuns.Boundary_ex(density_com,boundary_judge_y,'y',n,'nearest');
x_ex = (x(1)-n*h):h:(x(end)+n*h);nx_ex = length(x_ex);
y_ex = (y(1)+n*h):-h:(y(end)-n*h);ny_ex = length(y_ex);
boundary_judge_x = ones(ny,nx_ex);
boundary_judge_y = ones(ny_ex,nx);
boundary_judge_D_x = gfuns.Boundary_value(x_ex,y,boundary_judge_x,boundary_D,0);
boundary_judge_D_y = gfuns.Boundary_value(x,y_ex,boundary_judge_y,boundary_D,0);
boundary_judge_OD_x = gfuns.Boundary_value(x_ex,y,boundary_judge_D_x,boundary_O,0);
boundary_judge_OD_y = gfuns.Boundary_value(x,y_ex,boundary_judge_D_y,boundary_O,0);
for i = 1:3
    Q_cell_xex{i} = Q_cell_xex{i}.*boundary_judge_OD_x+Q_cell_xex_D{i}.*(1-boundary_judge_OD_x);
    Q_cell_yex{i} = Q_cell_yex{i}.*boundary_judge_OD_y+Q_cell_yex_D{i}.*(1-boundary_judge_OD_y);
end
density_com_xex = density_com_xex.*boundary_judge_OD_x+density_com_xex_D.*(1-boundary_judge_OD_x);
density_com_yex = density_com_yex.*boundary_judge_OD_y+density_com_yex_D.*(1-boundary_judge_OD_y);

temp = Q_cell_yex{2};
Q_cell_yex{2} = Q_cell_yex{3};
Q_cell_yex{3} = temp;

% F_vector
F_cell_ex = cell(1,3);G_cell_ex = cell(1,3);
F_cell_ex{1} = Q_cell_xex{2};
F_cell_ex{2} = Q_cell_xex{2}.^2 ./ Q_cell_xex{1};
F_cell_ex{3} = Q_cell_xex{2}.*Q_cell_xex{3}./ Q_cell_xex{1};
F_cell_ex{2}(isnan(F_cell_ex{2})) = 0;F_cell_ex{3}(isnan(F_cell_ex{3})) = 0;
F_cell_ex{2} = F_cell_ex{2} + gfuns.Force_ic(Q_cell_xex{1},density_com_xex);
% G_vector
G_cell_ex{1} = Q_cell_yex{2};
G_cell_ex{2} = Q_cell_yex{2}.^2 ./ Q_cell_yex{1};
G_cell_ex{3} = Q_cell_yex{2}.*Q_cell_yex{3}./ Q_cell_yex{1};
G_cell_ex{2}(isnan(G_cell_ex{2})) = 0;G_cell_ex{3}(isnan(G_cell_ex{3})) = 0;
G_cell_ex{2} = G_cell_ex{2} + gfuns.Force_ic(Q_cell_yex{1},density_com_yex);
% Obtain deviation of density
Q_xex = zeros(ny,nx_ex,3);F_vector = zeros(ny,nx_ex,3);
Q_yex = zeros(ny_ex,nx,3);G_vector = zeros(ny_ex,nx,3);
F_cell_mid = cell(1,3);G_cell_mid = cell(1,3);
for i = 1:3
    Q_xex(:,:,i) = Q_cell_xex{i};
    Q_yex(:,:,i) = Q_cell_yex{i};
    F_vector(:,:,i) = F_cell_ex{i};
    G_vector(:,:,i) = G_cell_ex{i};
end

[F_right_mid,alpha_F] = Reconstruction_Flow(Q_xex,density_com_xex,F_vector,[0 -1 0]);
[G_right_mid,alpha_G] = Reconstruction_Flow(Q_yex,density_com_yex,G_vector,[1  0 0]);
alpha = max(alpha_F,alpha_G);
for i = 1:3
    F_cell_mid{i} = F_right_mid(:,2:end-2,i);
    G_cell_mid{i} = G_right_mid(3:end-1,:,i);
end
% BC on physical boundary
F_cell_mid{1}(:,2:end) = gfuns.Boundary_value(x,y,F_cell_mid{1}(:,2:end),boundary_H,0);
F_cell_mid{1}(:,1:end-1) = gfuns.Boundary_value(x,y,F_cell_mid{1}(:,1:end-1),boundary_H,0);
G_cell_mid{1}(1:end-1,:) = gfuns.Boundary_value(x,y,G_cell_mid{1}(1:end-1,:),boundary_H,0);
G_cell_mid{1}(2:end,:) = gfuns.Boundary_value(x,y,G_cell_mid{1}(2:end,:),boundary_H,0);
% BC on origin boundary

if D_in(1) ~= 0
    % boundary_judge_x = ones(ny,nx+1);
    x_ex = (x(1)-0.5*h):h:(x(end)+0.5*h);
    y_ex = (y(1)+0.5*h):-h:(y(end)-0.5*h);
    % boundary_judge_Ox = gfuns.Boundary_value(x,y,boundary_judge_x,boundary_H,0);
    F_cell_mid{1}  = gfuns.Boundary_value(x_ex,y,F_cell_mid{1},boundary_O,D_in(2));
    G_cell_mid{1}  = gfuns.Boundary_value(x,y_ex,G_cell_mid{1},boundary_O,D_in(3));
    F_cell_mid{2} = gfuns.Boundary_value(x_ex,y,F_cell_mid{2},boundary_O,D_in(2).^2 ./ D_in(1)+gfuns.Force_ic(D_in(1),D_in(1)));
    F_cell_mid{3} = gfuns.Boundary_value(x_ex,y,F_cell_mid{3},boundary_O,D_in(2).*D_in(3)./D_in(1));
    G_cell_mid{2} = gfuns.Boundary_value(x,y_ex,G_cell_mid{2},boundary_O,D_in(3).^2 ./ D_in(1)+gfuns.Force_ic(D_in(1),D_in(1)));
    G_cell_mid{3} = gfuns.Boundary_value(x,y_ex,G_cell_mid{3},boundary_O,D_in(2).*D_in(3)./D_in(1));
end
% Boundary conditions
if (t>=t_bound(1))&&(t<=t_bound(2))
    if any(i_OD == [1,2,7])
        y_ob1 = (area(2,2) - y_bound(2))./h +1;
        x_ob1 = x_bound(1)./h + 1;
        x_ob2 = x_bound(2)./h + 0;
        q1_out = Q_cell{1}(y_ob1,x_ob1:x_ob2)./density_com(y_ob1,x_ob1:x_ob2).*d_out;
        v_out = Q_cell{3}(y_ob1,x_ob1:x_ob2)./q1_out; v_out(isnan(v_out)) = 0;
%         v_out = 0./d_out;
        G_cell_mid{1}(y_ob1,x_ob1:x_ob2) = v_out.*q1_out;
        G_cell_mid{2}(y_ob1,x_ob1:x_ob2) = (v_out.^2).*q1_out+gfuns.Force_ic(q1_out,d_out);
        G_cell_mid{3}(y_ob1,x_ob1:x_ob2) = 0;
    end
end

if (t>=t_panic(1))&&(t<=t_panic(end))
    x_ex = (x(1)-0.5*h):h:(x(end)+0.5*h);
    y_ex = (y(1)+0.5*h):-h:(y(end)-0.5*h);
    switch i_OD
        case 5
            boundary_judge_PC_mid_x = zeros(ny,nx+1);
            boundary_judge_PC_mid_y = zeros(ny+1,nx);
            boundary_judge_PC_mid_x = gfuns.Boundary_value(x_ex,y,boundary_judge_PC_mid_x,boundary_Pole,1);
            boundary_judge_PC_mid_y = gfuns.Boundary_value(x,y_ex,boundary_judge_PC_mid_y,boundary_Pole,1);
            F_cell_mid{1}(boundary_judge_PC_mid_x==1) = max(-1/10,min(1/10,F_cell_mid{1}(boundary_judge_PC_mid_x==1)));
            G_cell_mid{1}(boundary_judge_PC_mid_y==1) = max(-1/10,min(1/10,G_cell_mid{1}(boundary_judge_PC_mid_y==1)));
        case 6
            boundary_judge_PC_mid_x = zeros(ny,nx+1);
            boundary_judge_PC_mid_y = zeros(ny+1,nx);
            boundary_judge_PC_mid_x = gfuns.Boundary_value(x_ex,y,boundary_judge_PC_mid_x,boundary_Container,1);
            boundary_judge_PC_mid_y = gfuns.Boundary_value(x,y_ex,boundary_judge_PC_mid_y,boundary_Container,1);
            F_cell_mid{1}(boundary_judge_PC_mid_x==1) = max(-1/10,min(1/10,F_cell_mid{1}(boundary_judge_PC_mid_x==1)));
            G_cell_mid{1}(boundary_judge_PC_mid_y==1) = max(-1/10,min(1/10,G_cell_mid{1}(boundary_judge_PC_mid_y==1)));
    end
%     F_cell_mid{2}(boundary_judge_PC_mid_x==1) = (v_out.^2).*q_out+gfuns.Force_ic(Q_cell{1}(y_ob1,x_ob1:x_ob2),density_com(y_ob1,x_ob1:x_ob2));
%     F_cell_mid{3}(boundary_judge_PC_mid_x==1) = 0;
end



div{1} = - 1 ./ h .* (F_cell_mid{1}(:,2:end) - F_cell_mid{1}(:,1:end-1) + G_cell_mid{1}(1:end-1,:) - G_cell_mid{1}(2:end,:));
div{2} = - 1 ./ h .* (F_cell_mid{2}(:,2:end) - F_cell_mid{2}(:,1:end-1) + G_cell_mid{3}(1:end-1,:) - G_cell_mid{3}(2:end,:));
div{3} = - 1 ./ h .* (F_cell_mid{3}(:,2:end) - F_cell_mid{3}(:,1:end-1) + G_cell_mid{2}(1:end-1,:) - G_cell_mid{2}(2:end,:));



for i = 1:3
    div{i} = gfuns.Boundary_value(x,y,div{i},boundary_H,0);
    div{i} = gfuns.Boundary_value(x,y,div{i},boundary_D,0);
    
end

% tend = toc(tstart);
% fprintf('Reconstruction. Computing time is %f\n',tend);
end

function [F_p,alpha_m] = Reconstruction_Flow(Q,density_com,F,ind)
global density_m b_c density_a v_free gfuns
% Compute nvmerical fluxes at cell 'i' interfaces.
% Set Variables
ny = length(F(:,1,1)); nx = length(F(1,:,1));
c_cell = gfuns.Force_ig(Q(:,:,1),density_com);

if ind(1)==0
    F_p1 =circshift(F,ind);
    Q_p1 =circshift(Q,ind);
    density_com_p1 = circshift(density_com,ind);
    density_p = (Q(:,:,1)+Q_p1(:,:,1))/2;
    density_com_p = (density_com_p1+density_com)/2;
    vx = Q(:,:,2)./Q(:,:,1); vx(isnan(vx)) = 0; 
    vx_p =( Q(:,:,2)+Q_p1(:,:,2) ) ./density_p ./2;  vx_p(density_p==0) = 0;
    vy_p =( Q(:,:,3)+Q_p1(:,:,3) ) ./density_p ./2;  vy_p(density_p==0) = 0;
    
    c_p = gfuns.Force_ig(density_p,density_com_p);
    % LLF
    alpha_H = abs(c_cell).*(density_com<=density_a) + abs(vx).*(density_com<=density_a);
    alpha_H = max(alpha_H,circshift(alpha_H,ind));
%     alpha_H = (max(abs(c_cell).*(density_com<=density_a),[],2) + max(abs(vx).*(density_com<=density_a),[],2)) * ones(1,nx);
    
    M_E = max(max(0, (-2.*vx+2.*sqrt(max(0,vx.^2-c_cell))).*(density_com>density_a) ),[],2) * ones(1,nx);
    alpha_E = max(abs(sqrt(max(0,M_E.^2-4.*M_E.*vx+4.*c_cell))+2.*vx-M_E)./2.*(density_com>density_a),[],2) * ones(1,nx);
    
else
    
    temp = Q(:,:,2); Q(:,:,2) = Q(:,:,3); Q(:,:,3) = temp;
    temp = F(:,:,2); F(:,:,2) = F(:,:,3); F(:,:,3) = temp;
    density_com_p1 = circshift(density_com,ind);
    F_p1 =circshift(F,ind);
    Q_p1 =circshift(Q,ind);
    density_p = (Q(:,:,1)+Q_p1(:,:,1))/2;
    density_com_p = (density_com_p1+density_com)/2;
    vx = Q(:,:,2)./Q(:,:,1); vx(isnan(vx)) = 0; 
    vx_p =( Q(:,:,2)+Q_p1(:,:,2) ) ./density_p ./2;  vx_p(density_p==0) = 0;
    vy_p =( Q(:,:,3)+Q_p1(:,:,3) ) ./density_p ./2;  vy_p(density_p==0) = 0;
    
    c_p = gfuns.Force_ig(density_p,density_com_p);
    % LLF
    alpha_H = abs(c_cell).*(density_com<=density_a) + abs(vx).*(density_com<=density_a);
    alpha_H = max(alpha_H,circshift(alpha_H,ind));
%     alpha_H = ones(ny,1) * (max(abs(c_cell).*(density_com<=density_a),[],1) + max(abs(vx).*(density_com<=density_a),[],1));
    
    M_E = ones(ny,1) * max(max(0, (-2.*vx+2.*sqrt(max(0,vx.^2-c_cell))).*(density_com>density_a) ),[],1);
    alpha_E = ones(ny,1) * max(abs(sqrt(max(0,M_E.^2-4.*M_E.*vx+4.*c_cell))+2.*vx-M_E)./2.*(density_com>density_a),[],1);
    
end

F_p = zeros(size(F)); 
for i = 1:ny
    for j = 1:nx
        % LF splitting
        if (density_com_p(i,j)~=0)
            if (density_com_p(i,j)<=density_a)
    %             F_p1_s_m = 1/2* ( F_p1(i,j,:) - alpha_H(i,j).*Q_p1(i,j,:) );
    %             F_s_p    = 1/2* ( F(i,j,:)    + alpha_H(i,j).*Q(i,j,:)    );
    %             F_p(i,j,:) = F_s_p + F_p1_s_m;

                c_p(i,j) = max(10^-6,c_p(i,j));
                L_p1 = [(vx_p(i,j)+c_p(i,j))./2./c_p(i,j) -1./2./c_p(i,j) 0; -vy_p(i,j) 0 1; -(vx_p(i,j)-c_p(i,j))./2./c_p(i,j) 1./2./c_p(i,j) 0];
                F_s    = L_p1*reshape(F(i,j,:),[3,1]);    Q_s    = L_p1*reshape(Q(i,j,:),[3,1]); 
                F_p1_s = L_p1*reshape(F_p1(i,j,:),[3,1]); Q_p1_s = L_p1*reshape(Q_p1(i,j,:),[3,1]); 

                F_p1_s_m = 1/2* ( F_p1_s - alpha_H(i,j).*Q_p1_s );
                F_s_p    = 1/2* ( F_s    + alpha_H(i,j).*Q_s    );
                K_p = F_s_p + F_p1_s_m;

                R_p1 = [1 0 1; vx_p(i,j)-c_p(i,j) 0 vx_p(i,j)+c_p(i,j); vy_p(i,j) 1 vy_p(i,j)];
                F_p(i,j,:) = reshape(R_p1*K_p,[1,1,3]);

            else

                alpha_temp = reshape([alpha_E(i,j)+M_E(i,j); alpha_E(i,j) ; 0],[1,1,3]);
                F_p1_s_m = 1/2* ( F_p1(i,j,:) - alpha_temp.*Q_p1(i,j,:) );
                F_s_p    = 1/2* ( F(i,j,:)    + alpha_temp.*Q(i,j,:)    );
                F_p(i,j,:) = F_s_p + F_p1_s_m;

            end
        end
    end
end
alpha_m = max(max(max(alpha_H)),max(max(alpha_E)));

if ind(2)==0
    temp = F_p(:,:,2); F_p(:,:,2) = F_p(:,:,3); F_p(:,:,3) = temp;
end


end
