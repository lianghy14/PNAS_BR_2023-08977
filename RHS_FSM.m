function [vep_x,vep_y] = RHS_FSM(Q_cell, tao_ex)
global x y h nx ny
global n_OD boundary_D boundary_H boundary_O boundary_D_outside
global density_0 density_m m_ave sigma v_free thita
global gfuns pfuns
global t D_in pressure
global panic_xy panic_ex x_ex y_ex nx_ex ny_ex

%% Cost potential: FSM
% Obtain magnitude of velosity
method = 'zero';% 'linear'/'nearest'/'pchip'/'spline'
density_com_ex = zeros(ny_ex,nx_ex);
density_weight_ex = zeros(ny_ex,nx_ex);
density_ex = cell(1,n_OD);
velosity_ex = cell(1,n_OD);
potential = cell(1,n_OD);
for i = 1:n_OD
    density_ex{i} = gfuns.Boundary_ex(Q_cell{i}{1},ones(ny,nx),'x',2,method);
    density_ex{i} = gfuns.Boundary_ex(density_ex{i},ones(ny,nx+4),'y',2,method);
    density_ex{i} = gfuns.Boundary_value(x_ex,y_ex,density_ex{i},boundary_O{i},D_in{i}(1));
    density_com_ex = density_com_ex + density_ex{i};
    density_weight_ex = density_weight_ex + density_ex{i}.*panic_ex(:,:,i);
end

for i = 1:n_OD
    velosity_ex{i} = gfuns.Velosity(density_com_ex,panic_ex(:,:,i).*gfuns.Panic(t));
    for j = 1:n_OD
        cosangle = (Q_cell{i}{2}.*Q_cell{j}{2} + Q_cell{i}{3}.*Q_cell{j}{3}) ./ max(10^-12,(sqrt(Q_cell{j}{2}.^2+Q_cell{j}{3}.^2) .* sqrt(Q_cell{i}{2}.^2+Q_cell{i}{3}.^2)));
        cosangle(isnan(cosangle)) = 1;
        cosangle_ex = gfuns.Boundary_ex(cosangle,ones(ny,nx),'x',2,'linear');
        cosangle_ex = gfuns.Boundary_ex(cosangle_ex,ones(ny,nx+4),'y',2,'linear');
        velosity_ex{i} = velosity_ex{i}.*gfuns.Velosity_c(density_ex{j},cosangle_ex);
    end
    cost_1 = 1./velosity_ex{i};
    cost_1(isnan(cost_1)) = 10.^4;
%     cost_2 = 0.2 .* max(0,density_com_ex);
    cost_2 = 0.02 .* (max(0,density_com_ex)).^2;
    cost =  cost_1 + cost_2;
    cost = gfuns.Boundary_value(x_ex,y_ex,cost,boundary_H,10.^4);
    potential{i} = fast_Godunov(cost,boundary_D{i},x_ex,y_ex);
end

% number = pfuns.potential(potential{1}(3:end-2,3:end-2),'potential');
%% Pressure distribution
boundary_judge = ones(ny_ex,nx_ex);
method = 'linear';% 'linear'/'nearest'/'pchip'/'spline'
vel_x_com  = zeros(ny_ex,nx_ex);
vel_y_com  = zeros(ny_ex,nx_ex);
ven_x_com  = zeros(ny_ex,nx_ex);
ven_y_com  = zeros(ny_ex,nx_ex);
ve_x = cell(1,n_OD);
ve_y = cell(1,n_OD);
boundary_judge_H = gfuns.Boundary_value(x_ex,y_ex,boundary_judge,boundary_H,0);
for i = 1:n_OD
    boundary_judge_DH = gfuns.Boundary_value(x_ex,y_ex,boundary_judge_H,boundary_D_outside{i},0);
    % Derive speed direction from potential
    potential_x_ex = gfuns.Boundary_ex(potential{i},boundary_judge_DH,'x',1,method);
    potential_y_ex = gfuns.Boundary_ex(potential{i},boundary_judge_DH,'y',1,method);
    phi_x = -(potential_x_ex(:,3:end)-potential_x_ex(:,1:end-2))./2./h;
    phi_y = -(potential_y_ex(1:end-2,:)-potential_y_ex(3:end,:))./2./h;
    phy_grad = sqrt(phi_x.^2+phi_y.^2);
    ve_x{i} = phi_x./phy_grad;
    ve_y{i} = phi_y./phy_grad;
    ve_x{i}(isnan(ve_x{i})) = 0;
    ve_y{i}(isnan(ve_y{i})) = 0;
    ven_x_com = ven_x_com + density_ex{i}.*ve_x{i}./ density_com_ex;
    ven_y_com = ven_y_com + density_ex{i}.*ve_y{i}./ density_com_ex;
    vel_x_com = vel_x_com + density_ex{i}.*ve_x{i}.*velosity_ex{i}./ density_com_ex;
    vel_y_com = vel_y_com + density_ex{i}.*ve_y{i}.*velosity_ex{i}./ density_com_ex;
end
vel_x_com(isnan(vel_x_com)) = 0;
vel_y_com(isnan(vel_y_com)) = 0;
ve_x_com = vel_x_com ./ sqrt(vel_x_com.^2+vel_y_com.^2);
ve_y_com = vel_y_com ./ sqrt(vel_x_com.^2+vel_y_com.^2);

for i = 1:n_OD
    velosity_ex{i} = gfuns.Velosity(density_com_ex,panic_ex(:,:,i).*gfuns.Panic(t));
    for j = 1:n_OD
        cosangle = ve_x{i}.*ve_x{j} + ve_y{i}.*ve_y{j};
        velosity_ex{i} = velosity_ex{i}.*gfuns.Velosity_c(density_ex{j},cosangle);
    end
end


% Derive density gradient
density_com_x_ex = gfuns.Boundary_ex(density_com_ex,boundary_judge_H,'x',1,method);
density_com_y_ex = gfuns.Boundary_ex(density_com_ex,boundary_judge_H,'y',1,method);
rho_x = (density_com_x_ex(:,3:end)-density_com_x_ex(:,1:end-2))./2./h;
rho_y = (density_com_y_ex(1:end-2,:)-density_com_y_ex(3:end,:))./2./h;
rho_grad = sqrt(rho_x.^2+rho_y.^2);
rho_x = rho_x./rho_grad;
rho_y = rho_y./rho_grad;

% Derive pushing force direction
Angle = acos(ve_x_com.*rho_x+ve_y_com.*rho_y); Angle(isnan(Angle)) = 0;
Angle = gfuns.Boundary_value(x_ex,y_ex,Angle,boundary_H,0);
Re = zeros(ny_ex,nx_ex);
Re(Angle<=pi/2) = 1; Re(Angle>pi/2) = 2;

% Derive pushing force magnitude
alpha = max(0, (density_com_ex - density_0) / (density_m - density_0));
% alpha = ones(size(density_com_ex));
% alpha = 1+Angle./pi.*(beta-1);
comflict = sqrt(ven_x_com.^2+ven_y_com.^2);

% comflict = 1;
force_p = gfuns.Panic(t).*gfuns.Force_p(density_com_ex).*comflict.*max(panic_ex(:,:,5),panic_ex(:,:,6));
force_pa = force_p./alpha; force_pa(isnan(force_pa)) = 0;
pressure_old = zeros(ny_ex,nx_ex);
pressure = 10.^12 .* ones(ny_ex,nx_ex);
pressure_a = 10.^12 .* ones(ny_ex,nx_ex);
pressure = gfuns.Boundary_value(x_ex,y_ex,pressure,boundary_O{1},0);
pressure_a = gfuns.Boundary_value(x_ex,y_ex,pressure_a,boundary_O{1},0);
pressure = gfuns.Boundary_value(x_ex,y_ex,pressure,boundary_O{2},0);
pressure_a = gfuns.Boundary_value(x_ex,y_ex,pressure_a,boundary_O{2},0);
boundary_judge_HO = gfuns.Boundary_value(x_ex,y_ex,boundary_judge_H,boundary_O{1},0);
boundary_judge_HO = gfuns.Boundary_value(x_ex,y_ex,boundary_judge_HO,boundary_O{2},0);
iteration = 0;
GS_it = {[1 nx_ex;1 ny_ex],[nx_ex 1;1 ny_ex],[nx_ex 1;ny_ex 1],[1 nx_ex;ny_ex 1]};
GS_i = [1 1;-1 1;-1 -1;1 -1];
while (norm(pressure - pressure_old)>=sigma)
    iteration = iteration + 1;
    pressure_old = pressure;
    for it = 1:4
        for i = GS_it{it}(2,1):GS_i(it,2):GS_it{it}(2,2)
            for j = GS_it{it}(1,1):GS_i(it,1):GS_it{it}(1,2)
                if boundary_judge_HO(i,j) ~= 0
                    switch Re(i,j)
                        case 1
                            px = [pressure(i,max(1,j-1)) pressure(i,j) pressure(i,min(nx_ex,j+1))];
                            py = [pressure(max(1,i-1),j) pressure(i,j) pressure(min(ny,i+1),j)];
                            px_min = min(px(1),px(3));
                            py_min = min(py(1),py(3));
                            judge = (abs(px_min-py_min)>=(force_p(i,j).*h));
                            pressure_new_1 = min(px_min,py_min) + force_p(i,j).*h;
                            pressure_new_2 = (px_min+py_min+(2.*(force_p(i,j).^2).*(h.^2)-(px_min-py_min).^2).^0.5)./2;
                            pressure(i,j) = judge .* pressure_new_1 + (1-judge) .* pressure_new_2;
                            if alpha(i,j) == 0
                                pressure_a(i,j) = 0;
                                pressure(i,j) = 0;
                            else
                                pressure_a(i,j) = pressure(i,j)/alpha(i,j);
                            end
                        case 2
                            px = [pressure_a(i,max(1,j-1)) pressure_a(i,j) pressure_a(i,min(nx_ex,j+1))];
                            py = [pressure_a(max(1,i-1),j) pressure_a(i,j) pressure_a(min(ny,i+1),j)];
                            px_min = min(px(1),px(3));
                            py_min = min(py(1),py(3));
                            judge = (abs(px_min-py_min)>=(force_pa(i,j).*h));
                            pressure_new_1 = min(px_min,py_min) + force_pa(i,j).*h;
                            pressure_new_2 = (px_min+py_min+(2.*(force_pa(i,j).^2).*(h.^2)-(px_min-py_min).^2).^0.5)./2;
                            pressure_a(i,j) = judge .* pressure_new_1 + (1-judge) .* pressure_new_2;
                            pressure(i,j) = pressure_a(i,j)*alpha(i,j);
                    end
                end
            end
        end
    end
%     if mod(iteration,10) == 0
%         error = norm(pressure - pressure_old);
%         fprintf('Iteration number is %d,Improved error is %f\n',iteration,error);
%     end
    if iteration >= 10000
        pressure = null;
        return;
    end
end
% tend = toc(tstart);
% fprintf('IC for potential. Computing time is %f\n-----------------------------------\n',tend);
% fprintf('IC for pressure. Iteration number is %d\n-----------------------------------\n',iteration);
pressure_x_expand = gfuns.Boundary_ex(pressure,boundary_judge_H,'x',1,'linear');
pressure_y_expand = gfuns.Boundary_ex(pressure,boundary_judge_H,'y',1,'linear');
[ap_x,ap_y] = gfuns.Velosity_ap(density_com_ex,(pressure_x_expand(:,3:end)-pressure_x_expand(:,1:end-2))./2./h,(pressure_y_expand(1:end-2,:)-pressure_y_expand(3:end,:))./2./h);
%% Final speed direction on RHS
vep_x = cell(1,n_OD); 
vep_y = cell(1,n_OD);
for i = 1:n_OD
    ve_x{i} = ve_x{i}.*velosity_ex{i};
    ve_y{i} = ve_y{i}.*velosity_ex{i};
    vep_x{i} = ve_x{i}(3:end-2,3:end-2) - ap_x(3:end-2,3:end-2).*tao_ex(3:end-2,3:end-2,i);
    vep_y{i} = ve_y{i}(3:end-2,3:end-2) - ap_y(3:end-2,3:end-2).*tao_ex(3:end-2,3:end-2,i);
    vep_x{i} = gfuns.Boundary_value(x,y,vep_x{i},boundary_H,0);
    vep_y{i} = gfuns.Boundary_value(x,y,vep_y{i},boundary_H,0);
end
pressure = pressure(3:end-2,3:end-2);
% number = pfuns.potential(pressure(3:end-2,3:end-2),'pressure(N/m)');
end