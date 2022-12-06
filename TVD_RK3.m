function [Q,alpha] = TVD_RK3(Q_cell,t)
global gfuns tstep i_OD n_OD h area nx ny tao_n tao_p
global boundary_D boundary_O
global D_in
global t_bound x_bound y_bound d_out
global panic_ex panic_xy


density_com = zeros(ny,nx);
for i_OD = [1,2,7]
    density_com = Q_cell{i_OD}{1} + density_com;
end
% for k = [1,2,7]
%     if (t>=t_bound(1))&&(t<=t_bound(2))
%         y_ob1 = (area(2,2) - y_bound(2))./h +1;
%         y_ob2 = (area(2,2) - y_bound(1))./h +0;
%         x_ob1 = x_bound(1)./h +1;
%         x_ob2 = x_bound(2)./h +1;
%     %     q_out = Q_cell{1}(y_ob1:y_ob2,x_ob1-1)./density_com(y_ob1:y_ob2,x_ob1-1).*d_out;
%     %     v_out = gfuns.Velosity(d_out,gfuns.Panic(t));
%         Q_cell{k}{1}(y_ob1:y_ob2,(x_ob1-1):(x_ob2-1)) = Q_cell{k}{1}(y_ob1:y_ob2,(x_ob1-1):(x_ob2-1))./density_com(y_ob1:y_ob2,(x_ob1-1):(x_ob2-1)).*d_out;
% %         Q_cell{k}{1}(y_ob1:y_ob2,(x_ob1-1):(x_ob2-1)) = max(Q_cell{k}{1}(y_ob1:y_ob2,(x_ob1-1):(x_ob2-1)),Q_cell{k}{1}(y_ob1:y_ob2,(x_ob1-1):(x_ob2-1))./density_com(y_ob1:y_ob2,(x_ob1-1):(x_ob2-1)).*d_out);
%     end
% end
% for k = [3,4,5,6]
%     if (t>=t_bound(1))&&(t<=t_bound(2))
%         y_ob1 = (area(2,2) - y_bound(2))./h +1;
%         y_ob2 = (area(2,2) - y_bound(1))./h +0;
%         x_ob1 = x_bound(1)./h +1;
%         x_ob2 = x_bound(2)./h +1;
%     %     q_out = Q_cell{1}(y_ob1:y_ob2,x_ob1-1)./density_com(y_ob1:y_ob2,x_ob1-1).*d_out;
%     %     v_out = gfuns.Velosity(d_out,gfuns.Panic(t));
% %         Q_cell{k}{1}(y_ob1:y_ob2,(x_ob1-1):(x_ob2-1)) = 0;
% %         Q_cell{k}{2}(y_ob1:y_ob2,(x_ob1-1):(x_ob2-1)) = 0;
% %         Q_cell{k}{3}(y_ob1:y_ob2,(x_ob1-1):(x_ob2-1)) = 0;
% %         Q_cell{k}{1}(y_ob1:y_ob2,(x_ob1-1):(x_ob2-1)) = max(Q_cell{k}{1}(y_ob1:y_ob2,(x_ob1-1):(x_ob2-1)),Q_cell{k}{1}(y_ob1:y_ob2,(x_ob1-1):(x_ob2-1))./density_com(y_ob1:y_ob2,(x_ob1-1):(x_ob2-1)).*d_out);
%     end
% end


tstart = tic;
D_in = cell(1,n_OD);
for i_OD = 1:n_OD
    D_in{i_OD} = gfuns.D_in(t,i_OD);
end

tao_ex = tao_n.*(1-gfuns.Panic(t).*panic_ex)+tao_p.*gfuns.Panic(t).*panic_ex;
tao = tao_n.*(1-gfuns.Panic(t).*panic_xy)+tao_p.*gfuns.Panic(t).*panic_xy;
[vep_x,vep_y] = RHS_FSM(Q_cell,tao_ex);

for k = 1:n_OD
    Q_cell{k}{2}(Q_cell{k}{1}<=0.1) = Q_cell{k}{1}(Q_cell{k}{1}<=0.1).*vep_x{k}(Q_cell{k}{1}<=0.1);
    Q_cell{k}{3}(Q_cell{k}{1}<=0.1) = Q_cell{k}{1}(Q_cell{k}{1}<=0.1).*vep_y{k}(Q_cell{k}{1}<=0.1);
end

tend_FSM = toc(tstart);
%% Numerical scheme

Q = cell(1,n_OD); alpha = zeros(1,n_OD);
density_com = zeros(ny,nx);
for i_OD = 1:n_OD
    density_com = Q_cell{i_OD}{1} + density_com;
end

for i_OD = 1:n_OD
    [alpha(i_OD),div_F] = Res_LFHE(Q_cell{i_OD},density_com,boundary_O{i_OD},boundary_D{i_OD},D_in{i_OD});
    for j = 1:3
        Q{i_OD}{j} = Q_cell{i_OD}{j}+tstep.*div_F{j};
    end
    Q{i_OD}{2} = Q{i_OD}{2}+tstep.*(Q_cell{i_OD}{1}.*vep_x{i_OD}-Q_cell{i_OD}{2})./tao(:,:,i_OD);
    Q{i_OD}{3} = Q{i_OD}{3}+tstep.*(Q_cell{i_OD}{1}.*vep_y{i_OD}-Q_cell{i_OD}{3})./tao(:,:,i_OD);
end

alpha = max(alpha);

tend_LFHE = toc(tstart);
fprintf('LF, time is %.2f, Computing time (FSM,LF) is %.2f %.2f\n',t,tend_FSM,tend_LFHE);
end

