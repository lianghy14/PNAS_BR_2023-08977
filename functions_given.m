function gfuns = functions_given
gfuns.D_in = @D_in_cal;
gfuns.Panic = @Panic_cal;
gfuns.Velosity = @Velosity_cal;
gfuns.Velosity_c = @Velosity_c_cal;
gfuns.Force_ig = @Force_ig_cal;
gfuns.Force_ic = @Force_ic_cal;
gfuns.Force_p = @Force_p_cal;
gfuns.Velosity_ap = @Velosity_ap_cal;
gfuns.Panic = @Panic_cal;
gfuns.Boundary_ex = @Boundary_ex_cal;
gfuns.Boundary_ex2 = @Boundary_ex2_cal;
gfuns.Boundary_value = @Boundary_value_cal;
end
%% Panic
function Panic = Panic_cal(t)
global t_panic

if t<=t_panic(1)
    Panic = 0;
else
    if t<=t_panic(2)
        Panic = 1.*(t-t_panic(1)) ./ (t_panic(2)-t_panic(1));
    else
        if t<=t_panic(3)
            Panic = 1;
        else
            if t<=t_panic(4)
                Panic = 1.*(t_panic(4)-t) ./ (t_panic(4)-t_panic(3));
            else
                Panic = 0;
            end
        end
    end
end
end

%% Fundmental diagram
function Velosity = Velosity_cal(density,panic)
global v_free
gama_1 = -0.09.*(1-panic) - 0.04.*panic;
alpha = 2;
Velosity = v_free.*exp(gama_1.*max(0,density).^alpha);
% Velosity = min( v_free, max( 0, v_free.*(1-density./density_m) ) );
end
function Velosity = Velosity_c_cal(density_c,angle)
gama_2 = -0.025;
Velosity = exp(gama_2.*(1-angle).*max(0,density_c).^2);
% Velosity = min( v_free, max( 0, v_free.*(1-density./density_m) ) );
end
%% Flow in
function D_in = D_in_cal(t,i_OD)
global t_bound t_panic
D_in = zeros(1,3);
switch i_OD
    case 1 % A to C
        switch true
            case t<=60
                D_in(1) = 1.0*(t / 60);
                D_in(2) = Velosity_cal(D_in(1),0) .* D_in(1);
                D_in(3) = 0;
            case t>60 && t<=300
                D_in(1) = 1.0;
                D_in(2) = Velosity_cal(D_in(1),0) .* D_in(1);
                D_in(3) = 0;
            case t>300 && t<t_panic(1)
                D_in(1) = max(0.5,1.0*(6-t / 60));
                D_in(2) = Velosity_cal(D_in(1),0) .* D_in(1);
                D_in(3) = 0;
            case t>=t_panic(1)
                D_in(1) = 0;
                D_in(2) = 0;
                D_in(3) = 0;
        end
    case 2 % B to C
        switch true
            case t<=60
                D_in(1) = 1.1*(t / 60);
                D_in(2) = - Velosity_cal(D_in(1),0) .* D_in(1);
                D_in(3) = 0;
            case t>60 && t<=300
                D_in(1) = 1.1;
                D_in(2) = - Velosity_cal(D_in(1),0) .* D_in(1);
                D_in(3) = 0;
            case t>300 && t<t_panic(1)
                D_in(1) = max(0.5,1.1*(6-t / 60));
                D_in(2) = - Velosity_cal(D_in(1),0) .* D_in(1);
                D_in(3) = 0;
            case t>=t_panic(1)
                D_in(1) = 0;
                D_in(2) = 0;
                D_in(3) = 0;
        end
    case 3 % C to A
        switch true
            case t<=60
                D_in(1) = 0.4*(t / 60);
                D_in(2) = 0;
                D_in(3) = - Velosity_cal(D_in(1),0) .* D_in(1);
            case (t>60)&&(t<=t_bound(1)-30)
                D_in(1) = 0.4;
                D_in(2) = 0;
                D_in(3) = - Velosity_cal(D_in(1),0) .* D_in(1);
            case t>t_bound(1)-30
                D_in(1) =  max(0,0.4.*(t_bound(1)./30 - t ./ 30));
                D_in(2) = 0;
                D_in(3) = - Velosity_cal(D_in(1),0) .* D_in(1);
        end
    case 4 % C to B
        switch true
            case t<=60
                D_in(1) = 0.2*(t / 60);
                D_in(2) = 0;
                D_in(3) = - Velosity_cal(D_in(1),0) .* D_in(1);
            case (t>60)&&(t<=t_bound(1)-30)
                D_in(1) = 0.2;
                D_in(2) = 0;
                D_in(3) = - Velosity_cal(D_in(1),0) .* D_in(1);
            case t>t_bound(1)-30
                D_in(1) = max(0,0.2.*(t_bound(1)./30 - t ./ 30));
                D_in(2) = 0;
                D_in(3) = - Velosity_cal(D_in(1),0) .* D_in(1);
        end
end

end

% sqrt(h'(rho))
function Force_ig = Force_ig_cal(density,density_com)
global c0 c1 beta density_m density_0 density_a i_OD
switch i_OD
    case {1,2,3,4,7}
%         Force_ig = c0.*max(0,density_com).^beta;
        Force_ig = c0.*(density_com<=density_0) ...
            + c1.*(density_com>density_0).*(density_com<=density_a) ...
            + 0;% HE model
%         Force_ig = (c0.*max(0,density_com).^beta).*sqrt(max(0,density./density_com)).*(density_com<=density_0) ...
%             + (c1.*(density_com-density_a)).*sqrt(max(0,density./density_com)).*(density_com>density_0).*(density_com<=density_a) ...
%             + 0.*(density_com>density_a);% HE model
        Force_ig(isnan(Force_ig)) = 0;
    case {5,6}
%         Force_ig = c0.*max(0,density_com).^beta;
%         Force_ig = (c0.*max(0,density_com).^beta).*sqrt(max(0,density./density_com));
        Force_ig = c0.*(density_com<=density_0) ...
            + c1.*(density_com>density_0).*(density_com<=density_a) ...
            + 0;% HE model
%         Force_ig = (c0.*max(0,density_com).^beta).*sqrt(max(0,density./density_com)).*(density_com<=density_0) ...
%             + (c1.*(density_com-density_a)).*sqrt(max(0,density./density_com)).*(density_com>density_0).*(density_com<=density_a) ...
%             + 0.*(density_com>density_a);% HE model
%         Force_ig(isnan(Force_ig)) = 0;
end
%         Force_ig = (c0.*max(0,density).^beta).*(density<=density_0) ...
%             + (c1.*(density-density_a)).*(density>density_0).*(density<=density_a) ...
%             + 0.*(density>density_a);% New model
%         Force_ig(isnan(Force_ig)) = 0;
end

function Force_ic = Force_ic_cal(density,density_com)
global k_g x_g beta density_0 density_a c0 c1 i_OD

switch i_OD
    case {1,2,3,4,7}
%         Force_ic = (c0.^2).*max(0,density);
        Force_ic = (c0.^2).*max(0,density).*(density<=density_0) ...
            +((c1.^2).*density + (c0.^2).*density_0 - (c1.^2).*density_0).*(density>density_0).*(density<=density_a) ...
            +((c1.^2).*density_a + (c0.^2).*density_0 - (c1.^2).*density_0).*(density>density_a);% New model
    case {5,6}
        Force_ic = (c0.^2).*max(0,density).*(density<=density_0) ...
            +((c1.^2).*density + (c0.^2).*density_0 - (c1.^2).*density_0).*(density>density_0).*(density<=density_a) ...
            +((c1.^2).*density_a + (c0.^2).*density_0 - (c1.^2).*density_0).*(density>density_a);% New model
end
% Force_ic = (k_g(1).*Force_ig_cal(x_g(1).*density,density_com-(1-x_g(1)).*density).^2 +...
%     k_g(2).*Force_ig_cal(x_g(2).*density,density_com-(1-x_g(2)).*density).^2+...
%     k_g(3).*Force_ig_cal(x_g(3).*density,density_com-(1-x_g(3)).*density).^2+...
%     k_g(4).*Force_ig_cal(x_g(4).*density,density_com-(1-x_g(4)).*density).^2+...
%     k_g(5).*Force_ig_cal(x_g(5).*density,density_com-(1-x_g(5)).*density).^2).*density;% New model

% Force_ic = (c0.^2).*(density + (density_com-density).*log((density_com-density)./density_com));
% Force_ic(isnan(Force_ic)) = 0;
% Force_ic = 1./2.*(c0.^2).*(density.^2);
% Force_ic = (c0.^2).*(density);
% Force_ic = (c0.^2).*max(0,density).*(density<=density_0) ...
%             +((c1.^2).*density + (c0.^2).*density_0 - (c1.^2).*density_0).*(density>density_0).*(density<=density_a) ...
%             +((c1.^2).*density_a + (c0.^2).*density_0 - (c1.^2).*density_0).*(density>density_a);% New model
% Force_ic = (1./(2.*beta+1).*c0.^2.*max(0,density).^(2.*beta+1)).*(density<=density_0) ...
%             +(1./(2.*beta+1).*c1.^2.*(density-density_a).^(2.*beta+1) + (1./(2.*beta+1).*c0.^2.*density_0.^(2.*beta+1)) - (1./(2.*beta+1).*c1.^2.*(density_0-density_a).^(2.*beta+1))).*(density>density_0).*(density<=density_a) ...
%             +(1./(2.*beta+1).*c0.^2.*density_0.^(2.*beta+1) - (1./(2.*beta+1).*c1.^2.*(density_0-density_a).^(2.*beta+1))).*(density>density_a);% New model
% Force_ic = 1./6.*0 + 2./3.*Force_ig_cal(1./2.*density,density_com-1./2.*density)+ 1./6.*Force_ig_cal(density,density_com);% New model
end


%% Pushing force
function Force_p = Force_p_cal(density,panic_xy)
global m_ave density_0 t;
% Force_p = 0;
% Force_p = m_ave .* 5.*(density-density_0).^0.5 .* (density>=density_0);
% Force_p = m_ave .* 2 .*sqrt(density-density_0) .* (density>=density_0);
Force_p = 120 .* Panic_cal(t).*sqrt(density-density_0) .* (density>=density_0);
end

function [ap_x,ap_y] = Velosity_ap_cal(density,pressure_x,pressure_y)
global m_ave
pressure_grad = sqrt(pressure_x.^2 + pressure_y.^2)./density./m_ave;
ap_x = pressure_x./density./m_ave;
ap_y = pressure_y./density./m_ave;
ap_x(isnan(ap_x)) = 0; ap_y(isnan(ap_y)) = 0;

end


%% Boundary expand
function cell_ex = Boundary_ex_cal(cell,cell_judge,dim,n,method)
nx = length(cell(1,:));
ny = length(cell(:,1));
switch dim
    case 'x'
        cell_ex = zeros(ny,nx+2*n);
        cell_judge = [cell_judge,zeros(ny,1)];
        for i = 1:ny
            n_loc = n+1;
            while n_loc<=(nx+n)
                if cell_judge(i,n_loc-n) == 0
                    cell_ex(i,n_loc) = cell(i,n_loc-n);
                    n_loc = n_loc+1;
                else
                    start_loc = n_loc-n;
                    end_loc = find(cell_judge(i,(n_loc-n):(nx+1))==0,1,'first')-1+n_loc-n-1;

                    switch method
                        case 'zero'
                            cell_ex(i,(start_loc):(start_loc+n-1)) = 0;
                            cell_ex(i,(start_loc+n):(end_loc+n)) =  cell(i,(start_loc):(end_loc));
                            cell_ex(i,(end_loc+n+1):(end_loc+2*n)) = 0;
                        case 'symmetric+'
                            cell_ex(i,(start_loc):(start_loc+n-1)) = cell(i,(start_loc+n-1):-1:(start_loc));
                            cell_ex(i,(start_loc+n):(end_loc+n)) =  cell(i,(start_loc):(end_loc));
                            cell_ex(i,(end_loc+n+1):(end_loc+2*n)) = cell(i,(end_loc):-1:(end_loc-n+1));
                        case 'symmetric-'
                            cell_ex(i,(start_loc):(start_loc+n-1)) = - cell(i,(start_loc+n-1):-1:(start_loc));
                            cell_ex(i,(start_loc+n):(end_loc+n)) =  cell(i,(start_loc):(end_loc));
                            cell_ex(i,(end_loc+n+1):(end_loc+2*n)) = - cell(i,(end_loc):-1:(end_loc-n+1));
                        otherwise
                            x = (n+1):(n+end_loc-start_loc+1);
                            x_expand = 1:(n+end_loc-start_loc+1+n);
                            cell_ex(i,(start_loc):(end_loc+2*n)) = interp1(x,cell(i,start_loc:end_loc),x_expand,method,'extrap');
                    end
                    n_loc = end_loc+2*n+1;
                end
            end
        end
    case 'y'
        cell_ex = zeros(ny+2*n,nx);
        cell_judge = [cell_judge;zeros(1,nx)];
        for i = 1:nx
            n_loc = n+1;
            while n_loc<=(ny+n)
                if cell_judge(n_loc-n,i) == 0
                    cell_ex(n_loc,i) = cell(n_loc-n,i);
                    n_loc = n_loc+1;
                else
                    start_loc = n_loc-n;
                    end_loc = find(cell_judge((n_loc-n):(ny+1),i)==0,1,'first')-1+n_loc-n-1;
                    switch method
                        case 'zero'
                            cell_ex((start_loc):(start_loc+n-1),i) = 0;
                            cell_ex((start_loc+n):(end_loc+n),i) =  cell((start_loc):(end_loc),i);
                            cell_ex((end_loc+n+1):(end_loc+2*n),i) = 0;
                        case 'symmetric+'
                            cell_ex((start_loc):(start_loc+n-1),i) = cell((start_loc+n-1):-1:(start_loc),i);
                            cell_ex((start_loc+n):(end_loc+n),i) =  cell((start_loc):(end_loc),i);
                            cell_ex((end_loc+n+1):(end_loc+2*n),i) = cell((end_loc):-1:(end_loc-n+1),i);
                        case 'symmetric-'
                            cell_ex((start_loc):(start_loc+n-1),i) = - cell((start_loc+n-1):-1:(start_loc),i);
                            cell_ex((start_loc+n):(end_loc+n),i) =  cell((start_loc):(end_loc),i);
                            cell_ex((end_loc+n+1):(end_loc+2*n),i) = - cell((end_loc):-1:(end_loc-n+1),i);
                        otherwise
                            y = (n+1):(n+end_loc-start_loc+1);
                            y_expand = 1:(n+end_loc-start_loc+1+n);
                            cell_ex((start_loc):(end_loc+2*n),i) = interp1(y,cell(start_loc:end_loc,i),y_expand,method,'extrap');
                    end
                    n_loc = end_loc+2*n+1;
                end
            end
        end
end

%% Simple zero
% switch dim
%     case 1
%         flow_expand = zeros(ny,nx+2*n);
%         flow_expand(:,(n+1):(n+nx)) = flow;
%         flow_expand(:,[1:n (n+nx+1):(n+nx+n)]) = 0;
%     case 2
%         flow_expand = zeros(ny+2*n,nx);
%         flow_expand((n+1):(n+ny),:) = flow;
%         flow_expand([1:n (n+ny+1):(n+ny+n)],:) = 0;
% end;

end

function cell_ex2 = Boundary_ex2_cal(cell,cell_judge,n,method)
nx = length(cell(1,:));
ny = length(cell(:,1));
cell_judge(cell_judge~=0) = 1;

cell_vector = zeros(sum(sum(cell_judge)),3);
count = 1;
for i = 1:ny
    for j = 1:nx
        if cell_judge(i,j) ~= 0
            cell_vector(count,:) = [j,i,cell(i,j)];
            count = count+1;
        end
    end
end
[X,Y] = meshgrid((1-n):(nx+n),(1-n):(ny+n));
cell_ex2 = griddata(cell_vector(:,1),cell_vector(:,2),cell_vector(:,3),X,Y,method);
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
