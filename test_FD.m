clc
clear
global density_m density_0 m_ave c0 c1 beta density_a b_c v_free k_g x_g i_OD
warning off

density_0 = 6; density_m = 10; tao = 0.61; gfuns = functions_given;
v_free = 1.034; gama = -0.09;  alpha = 2; i_OD = 6;
x_n1 = sqrt(245-14*sqrt(70))/21;
x_n2 = sqrt(245+14*sqrt(70))/21;
k_n1 = (322+13.*sqrt(70))./900;
k_n2 = (322-13.*sqrt(70))./900;
x_g = [-x_n2; -x_n1; 0; x_n1; x_n2]./2 + 0.5;
k_g = [k_n2; k_n1; 128/225; k_n1; k_n2]./2;


density = 0:0.1:density_m;
velosity = v_free.*exp(gama.*max(0,density).^alpha);

force_ex = max(0,90.*(density-density_0));
flow = velosity.*density;
ve_rho = alpha.*gama.*((density.^(alpha-1))).*min( v_free, max( 0, v_free.*exp(gama.*(density).^alpha) ) );
% force_i = m_ave.*c0.^2.*(3./2.*density+2.*density_m./pi.*sin(pi.*density./density_m)+1./4.*density_m./pi.*sin(2.*pi.*density./density_m));% PW model;
% force_i = 1./(2.*beta+1).*c0.^2.*density.^(2.*beta+1)./density_m.^(2.*beta);% PW model
% c_rho = c0.*(cos(pi.*density./density_m)+1);
% PW model
% force_i = m_ave.*1./(2.*beta+1).*c0.^2.*density.^(2.*beta+1)./density_m.^(2.*beta);
% c_rho = c0.*(density./density_m).^beta;

% Eliptic model
% RT = 1.0; b = 0.1; a = 0.4;
% force_i = m_ave.*(RT.*density./(1-b.*density)-a.*density.^2);
% c_rho = b.*RT.*density./(1-b.*density).^2 + RT./(1-b.*density) -2.*a.*density;
% c0 = 0.05; a = 4.0; b = 7.0;
% force_i =  m_ave.*c0.*(1./3.*density.^2-1./2.*(a+b).*density+a.*b).*density;
% c_rho = c0.*(density.^2-(a+b).*density+a.*b);

% c0 = 1.5; c1 = 1; density_a = 3.0; b_c = 5.0;
% force_i =  m_ave.*(c0.^2.*(density).*(density<=density_a) + (-c1.^2.*(density-density_a)+c0.^2.*density_a).*(density>density_a).*(density<=b_c) + (c0.^2.*(density-b_c)-c1.^2.*(b_c-density_a)+c0.^2.*density_a).*(density>b_c));
% c_rho = (c0).*(density<=density_a) - c1.*(density>density_a).*(density<=b_c) + (c0).*(density>b_c);

% c0 = 0.15; density_a = 4.0; b_c = 5.0;
% force_i =  m_ave.*c0.*((-(density-density_a).^2+density_a^2).*(density<=density_a) + density_a^2.*(density>density_a).*(density<=b_c) + ((density-b_c).^2+density_a^2).*(density>b_c));
% c_rho = c0.*((-(density-density_a).*2).*(density<=density_a) + 0.*(density>density_a).*(density<=b_c) + ((density-b_c).*2).*(density>b_c));
% 
m_ave = 60; density_0 = 6; density_m = 10; gfuns = functions_given;
c0 = 0.6; beta = 0; density_a = 7; c1 = c0./2; 

force_i = gfuns.Force_ic(density,density).*m_ave;
c_rho = gfuns.Force_ig(density,density);

figure;
plot(density,density.*velosity);hold on;
plot(density,velosity);
xlabel('density (ped/m^2)');ylabel('FD');
figure;
plot(density,abs(density.*ve_rho));
hold on;
plot(density,c_rho);
xlabel('density (ped/m^2)');ylabel('dv/d\rho & dc/d\rho');

figure;
plot(density,force_i);
% hold on;
% plot(density,force_ex);
% hold on;
% plot(density,density.*(velosity-c_rho));
xlabel('density (ped/m^2)');ylabel('P_1');
