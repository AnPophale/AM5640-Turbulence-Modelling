clc; clear all; close all; format compact; fopen shortg; 

%Mesh Parameters
L = 1;
N = 100; 

%Node locations
ynode = zeros(N+2,1);
ynode(2) = 0.15; %First node away from wall (yplus between 30 and 100)
for i = 3:N+2
    ynode(i) = ynode(i-1) + (1-ynode(2))/(N); %All other nodes uniformly spaced
 end

%Cell Face Locations
yface = zeros(N+1,1);
yface(1) = ynode(1);
yface(end) = ynode(end);
for i = 2:N
    yface(i) = (ynode(i) + ynode(i+1))/2;
end

%Delta y, Delta yn and Delta ys
dy = zeros(N,1);
dy_n = zeros(N,1);
dy_s = zeros(N,1);
for i = 1:N
    dy(i) = yface(i+1) - yface(i);
    dy_n(i) = ynode(i+2) - ynode(i+1);
    dy_s(i) = ynode(i+1) - ynode(i);
end

%Interpolation Functions
fn = zeros(N,1); fs = zeros(N,1);
for i = 1:N
    fn(i) = 0.5*dy(i)/dy_n(i);
    fs(i) = 0.5*dy(i)/dy_s(i);
end

%Physical Parameters
nu = 1/395;
rho = 1;

%Model Parameters
E = 9;
B = 5.2;
kappa = 0.41;
c_mu = 0.09;
c1 = 1.8;
c2 = 0.6;
c1_p = 0.5;
c2_p = 0.3;
c1_eps = 1.44;
c2_eps = 1.92;
sigma_k = 1;
sigma_eps = 1.3;

%Under Relaxation Factors
urf = 0.2;

%Residual limit
R = 1;
R_max = 0.0001;

%Initial Guess and Boundary Conditions
ustar = 1;

%Loading k w model results
load("y_dns.dat");
u_kw = readmatrix("u_kw.txt");
u2_kw = readmatrix("u2_kw.txt");
v2_kw = readmatrix("v2_kw.txt");
w2_kw = readmatrix("w2_kw.txt");
uv_kw = -1*readmatrix("uv_kw.txt");
eps_kw = readmatrix("eps_kw.txt");

%Inital guess interpolated from k w model results
%To run this code, the results of the k omega model must be stored in the same folder
%u2, v2 and w2 taken as k/3 from k omega model
u = interp1(y_dns, u_kw, ynode, 'linear'); u(1) = 0; 
u2 = interp1(y_dns, u2_kw, ynode, 'linear'); u2(1) = 0; 
v2 = interp1(y_dns, v2_kw, ynode, 'linear'); v2(1) = 0;
w2 = interp1(y_dns, w2_kw, ynode, 'linear'); w2(1) = 0;
uv = -1*interp1(y_dns, uv_kw, ynode, 'linear'); uv(1) = 0; uv(end) = 0;
eps = interp1(y_dns, eps_kw, ynode, 'linear');

nu_t = zeros(N+2,1);
nu_t_n = zeros(N,1); nu_t_s = zeros(N,1);

aP_u = zeros(N,1); aN_u = zeros(N,1); aS_u = zeros(N,1); Su_u = zeros(N,1);
aP_u2 = zeros(N,1); aN_u2 = zeros(N,1); aS_u2 = zeros(N,1); Su_u2 = zeros(N,1); Sp_u2 = zeros(N,1);
aP_v2 = zeros(N,1); aN_v2 = zeros(N,1); aS_v2 = zeros(N,1); Su_v2 = zeros(N,1); Sp_v2 = zeros(N,1);
aP_w2 = zeros(N,1); aN_w2 = zeros(N,1); aS_w2 = zeros(N,1); Su_w2 = zeros(N,1); Sp_w2 = zeros(N,1);
aP_uv = zeros(N,1); aN_uv = zeros(N,1); aS_uv = zeros(N,1); Su_uv = zeros(N,1);
aP_eps = zeros(N,1); aN_eps = zeros(N,1); aS_eps = zeros(N,1); Su_eps = zeros(N,1); Sp_eps = zeros(N,1);

Pk = zeros(N,1); P_11 = zeros(N,1); P_22 = zeros(N,1); P_33 = zeros(N,1); P_12 = zeros(N,1);
Pi_11_r = zeros(N,1); Pi_11_s = zeros(N,1); Pi_11_r_p = zeros(N,1); Pi_11_s_p = zeros(N,1);
Pi_22_r = zeros(N,1); Pi_22_s = zeros(N,1); Pi_22_r_p = zeros(N,1); Pi_22_s_p = zeros(N,1);
Pi_33_r = zeros(N,1); Pi_33_s = zeros(N,1); Pi_33_r_p = zeros(N,1); Pi_33_s_p = zeros(N,1);
Pi_12_r = zeros(N,1); Pi_12_s = zeros(N,1); Pi_12_r_p = zeros(N,1); Pi_12_s_p = zeros(N,1);
eps_11 = zeros(N,1); eps_22 = zeros(N,1); eps_33 = zeros(N,1); eps_12 = zeros(N,1);
dudy = zeros(N,1); duvdy = zeros(N,1); f = zeros(N,1);

iter = 0;
%Main iteration loop
while(R > R_max)
    iter = iter + 1;
    u_old = u; u2_old = u2; v2_old = v2; w2_old = w2; uv_old = uv; eps_old = eps;
    k_old = 0.5*(u2_old + v2_old + w2_old);
    
    residual_u = 0; residual_u2 = 0; residual_v2 = 0; residual_w2 = 0; residual_uv = 0; residual_eps = 0;

    %Solving for ustar iteratively
    err = 1;
    while err > 0.001
        ustar_old = ustar;
        ustar = kappa*u_old(2)/log(E*ustar_old*ynode(2)/nu);
        err = abs(ustar - ustar_old)/ustar_old;
    end

    %Wall Functions
    tau_wall = rho*ustar^2;
    u2(2) = 3.67*ustar^2;
    v2(2) = 0.83*ustar^2;
    w2(2) = 2.17*ustar^2;
    uv(2) = -ustar^2;
    eps(2) = ustar^3/(kappa*ynode(2));

    %Calculating turbulence viscosity
    nut(1) = 0;
    nu_t(2:end) = c_mu*k_old(2:end).^2./eps_old(2:end);
    
    %Interpolating Turbulence Viscosity to Cell Faces
    for i = 1:N
        nu_t_n(i) = fn(i)*nu_t(i+2) + (1-fn(i))*nu_t(i+1);
        nu_t_s(i) = fs(i)*nu_t(i) + (1-fs(i))*nu_t(i+1);
    end
    nu_t_n(end) = nu_t(end);
    nu_t_s(1) = nu_t(1);

    %Calculating all source terms
    for i = 1:N
        dudy(i) = (u_old(i+2)-u_old(i))/(dy_n(i) + dy_s(i));
        duvdy(i) = (uv_old(i+2)-uv_old(i))/(dy_n(i) + dy_s(i));

        %P_ij and Pk
        Pk(i) = -uv_old(i+1)*dudy(i);
        P_11(i) = -2*uv_old(i+1)*dudy(i);
        P_22(i) = 0;
        P_33(i) = 0;
        P_12(i) = -v2_old(i+1)*dudy(i);

        %Epsilon ij
        eps_11(i) = 2*eps_old(i+1)/3;
        eps_22(i) = 2*eps_old(i+1)/3;
        eps_33(i) = 2*eps_old(i+1)/3;
        eps_12(i) = 0;

        %Pi_ij rapid
        Pi_11_r(i) = -c2*rho*(P_11(i) - 2*Pk(i)/3);
        Pi_22_r(i) = -c2*rho*(P_22(i) - 2*Pk(i)/3);
        Pi_33_r(i) = -c2*rho*(P_33(i) - 2*Pk(i)/3);
        Pi_12_r(i) = -c2*rho*P_12(i);

        %Pi_ij slow
        Pi_11_s(i) = -c1*rho*(eps_old(i+1)/k_old(i+1))*(u2_old(i+1) - 2*k_old(i+1)/3);
        Pi_22_s(i) = -c1*rho*(eps_old(i+1)/k_old(i+1))*(v2_old(i+1) - 2*k_old(i+1)/3);
        Pi_33_s(i) = -c1*rho*(eps_old(i+1)/k_old(i+1))*(w2_old(i+1) - 2*k_old(i+1)/3);
        Pi_12_s(i) = -c1*rho*(eps_old(i+1)/k_old(i+1))*uv_old(i+1);

        %Pi_ij rapid correction
        f(i) = (k_old(i+1)^(3/2))/(2.55*eps_old(i+1)*ynode(i+1));
            
        Pi_11_r_p(i) = c2_p*f(i)*Pi_11_r(i);
        Pi_22_r_p(i) = -2*c2_p*f(i)*Pi_22_r(i);
        Pi_33_r_p(i) = c2_p*f(i)*Pi_33_r(i);
        Pi_12_r_p(i) = -(3/2)*c2_p*f(i)*Pi_12_r(i);

        %Pi_ij slow correction
        Pi_11_s_p(i) = c1_p*rho*(eps_old(i+1)/k_old(i+1))*v2_old(i+1)*f(i);
        Pi_22_s_p(i) = -2*c1_p*rho*(eps_old(i+1)/k_old(i+1))*v2_old(i+1)*f(i);
        Pi_33_s_p(i) = c1_p*rho*(eps_old(i+1)/k_old(i+1))*v2_old(i+1)*f(i);
        Pi_12_s_p(i) = -(3/2)*c1_p*rho*(eps_old(i+1)/k_old(i+1))*uv_old(i+1)*f(i);
    end

    %Solving the u momentum equation
    for i = 1:N
        aN_u(i) = nu/dy_n(i);
        aS_u(i) = nu/dy_s(i);
        aP_u(i) = aN_u(i) + aS_u(i);
        Su_u(i) = -duvdy(i)*dy(i) + dy(i)/rho; 
    end
    Su_u(1) = Su_u(1) - tau_wall;

    %Gauss Seidel method
    for i = 1:N
        u(i+1) = (aN_u(i)*u(i+2) + aS_u(i)*u(i) + Su_u(i))/aP_u(i);
    end
    u(end) = u(end-1); %Neumann BC at channel centerline
    u = (1-urf)*u + urf*u_old; %Under Relaxation
    %Residual Calcuation
    for i = 1:N
        residual_u = residual_u + abs(aN_u(i)*u(i+2) + aS_u(i)*u(i) + Su_u(i) - aP_u(i)*u(i+1));
    end

    %Solving the u2 equation
    for i = 2:N
        aN_u2(i) = (nu + nu_t_n(i)/sigma_k)/dy_n(i);
        aS_u2(i) = (nu + nu_t_s(i)/sigma_k)/dy_s(i);
        aP_u2(i) = aN_u2(i) + aS_u2(i);
        Su_u2(i) = (P_11(i) + Pi_11_s(i) + Pi_11_r(i) + Pi_11_s_p(i) + Pi_11_r_p(i))*dy(i);
        Sp_u2(i) = eps_11(i)*dy(i)/u2_old(i+1);
    end

    %Gauss Seidel method
    for i = 2:N
        u2(i+1) = (aN_u2(i)*u2(i+2) + aS_u2(i)*u2(i) + Su_u2(i))/(aP_u2(i) + Sp_u2(i));
    end
    u2(end) = u2(end-1); %Neumann BC at channel centerline
    u2 = (1-urf)*u2 + urf*u2_old; %Under Relaxation
    %Residual Calcuation
    for i = 2:N
        residual_u2 = residual_u2 + abs(aN_u2(i)*u2(i+2) + aS_u2(i)*u2(i) + Su_u2(i) - (aP_u2(i) + Sp_u2(i))*u2(i+1));
    end
    
    %Solving the v2 equation
    for i = 2:N
        aN_v2(i) = (nu + nu_t_n(i)/sigma_k)/dy_n(i);
        aS_v2(i) = (nu + nu_t_s(i)/sigma_k)/dy_s(i);
        aP_v2(i) = aN_v2(i) + aS_v2(i);
        Su_v2(i) = (P_22(i) + Pi_22_s(i) + Pi_22_r(i) + Pi_22_s_p(i) + Pi_22_r_p(i))*dy(i);
        Sp_v2(i) = eps_22(i)*dy(i)/v2_old(i+1);
    end

    %Gauss Seidel method
    for i = 2:N
        v2(i+1) = (aN_v2(i)*v2(i+2) + aS_v2(i)*v2(i) + Su_v2(i))/(aP_v2(i) + Sp_v2(i));
    end
    v2(end) = v2(end-1); %Neumann BC at channel centerline
    v2 = (1-urf)*v2 + urf*v2_old; %Under Relaxation
    %Residual Calcuation
    for i = 2:N
        residual_v2 = residual_v2 + abs(aN_v2(i)*v2(i+2) + aS_v2(i)*v2(i) + Su_v2(i) - (aP_v2(i) + Sp_v2(i))*v2(i+1));
    end

    %Solving the w2 equation
    for i = 2:N
        aN_w2(i) = (nu + nu_t_n(i)/sigma_k)/dy_n(i);
        aS_w2(i) = (nu + nu_t_s(i)/sigma_k)/dy_s(i);
        aP_w2(i) = aN_w2(i) + aS_w2(i);
        Su_w2(i) = (P_33(i) + Pi_33_s(i) + Pi_33_r(i) + Pi_33_s_p(i) + Pi_33_r_p(i))*dy(i);
        Sp_w2(i) = eps_33(i)*dy(i)/w2_old(i+1);
    end

    %Gauss Seidel method
    for i = 2:N
        w2(i+1) = (aN_w2(i)*w2(i+2) + aS_w2(i)*w2(i) + Su_w2(i))/(aP_w2(i) + Sp_w2(i));
    end
    w2(end) = w2(end-1); %Neumann BC at channel centerline
    w2 = (1-urf)*w2 + urf*w2_old; %Under Relaxation
    %Residual Calcuation
    for i = 2:N
        residual_w2 = residual_w2 + abs(aN_w2(i)*w2(i+2) + aS_w2(i)*w2(i) + Su_w2(i) - (aP_w2(i) + Sp_w2(i))*w2(i+1));
    end

    %Solving the uv equation
    for i = 2:N
        aN_uv(i) = (nu + nu_t_n(i)/sigma_k)/dy_n(i);
        aS_uv(i) = (nu + nu_t_s(i)/sigma_k)/dy_s(i);
        aP_uv(i) = aN_uv(i) + aS_uv(i);
        Su_uv(i) = (P_12(i) + Pi_12_s(i) + Pi_12_r(i) + Pi_12_s_p(i) + Pi_12_r_p(i))*dy(i);
    end

    %Gauss Seidel method
    for i = 2:N
        uv(i+1) = (aN_uv(i)*uv(i+2) + aS_uv(i)*uv(i) + Su_uv(i))/aP_uv(i);
    end
    uv(end) = 0; %Dirichilet BC at channel centerline
    uv = (1-urf)*uv + urf*uv_old; %Under Relaxation
    %Residual Calcuation
    for i = 2:N
        residual_uv = residual_uv + abs(aN_uv(i)*uv(i+2) + aS_uv(i)*uv(i) + Su_uv(i) - aP_uv(i)*uv(i+1));
    end

    %Solving the epsilon equation
    for i = 2:N
        aN_eps(i) = (nu + nu_t_n(i)/sigma_eps)/dy_n(i);
        aS_eps(i) = (nu + nu_t_s(i)/sigma_eps)/dy_s(i);
        aP_eps(i) = aN_eps(i) + aS_eps(i);
        Su_eps(i) = c1_eps*(eps_old(i+1)/k_old(i+1))*Pk(i)*dy(i);
        Sp_eps(i) = c2_eps*(eps_old(i+1)/k_old(i+1))*dy(i);
    end

    %Gauss Seidel method
    for i = 2:N
        eps(i+1) = (aN_eps(i)*eps(i+2) + aS_eps(i)*eps(i) + Su_eps(i))/(aP_eps(i) + Sp_eps(i));
    end
    eps(end) = eps(end-1); eps(1) = eps(2); %Neumann BC at wall and channel centerline
    eps = (1-urf)*eps + urf*eps_old; %Under Relaxation    
    %Residual Calcuation
    for i = 2:N
        residual_eps = residual_eps + abs(aN_eps(i)*eps(i+2) + aS_eps(i)*eps(i) + Su_eps(i) - (aP_eps(i) + Sp_eps(i))*eps(i+1));
    end

    %Calculating TKE
    k = 0.5*(u2 + v2 + w2);
    
    %Calculating the residual
    R = max([residual_u,residual_u2,residual_v2,residual_w2,residual_uv,residual_eps]);
end

%Calculating yplus
y_plus = ynode*ustar/nu;
fprintf('The value of ustar is: %f\n', ustar);
fprintf('Yplus for first node away from wall is : %f\n', y_plus(2));

%Loading DNS Data
load dns_data.dat
load u_dns.dat
load u2_dns.dat
load v2_dns.dat
load w2_dns.dat
load uv_dns.dat

%Calculating TKE from dns data
k_dns=0.5*(u2_dns+v2_dns+w2_dns);

%Epslion DNS
ustar_dns = 1;
epsilon_dns=dns_data(:,2)*ustar_dns^4/nu; 

%Plotting
%U velocity
figure(1)
plot(u_dns,y_dns,'LineWidth',1.25);
hold on
plot(u,ynode,'-x','LineWidth',1.25);
xlabel('u')
ylabel('y')
title('U velocity')
legend('DNS','RSM','Location','northwest');

%u2 Reynolds Stress
figure(2);
plot(y_dns,u2_dns,'LineWidth',1.25);
hold on
plot(ynode,u2,'-x','LineWidth',1.25);
xlabel('y')
ylabel('u2')
title('u2 Reynolds Stress')
legend('DNS','RSM');

%v2 Reynolds Stress
figure(3);
plot(y_dns,v2_dns,'LineWidth',1.25);
hold on
plot(ynode,v2,'-x','LineWidth',1.25);
xlabel('y')
ylabel('v2')
title('v2 Reynolds Stress')
legend('DNS','RSM');

%w2 Reynolds Stress
figure(4);
plot(y_dns,w2_dns,'LineWidth',1.25);
hold on
plot(ynode,w2,'-x','LineWidth',1.25);
xlabel('y')
ylabel('w2')
title('w2 Reynolds Stress')
legend('DNS','RSM');

%Reynolds Shear Stress
figure(5);
plot(y_dns,-uv_dns,'LineWidth',1.25);
hold on
plot(ynode,-uv,'-x','LineWidth',1.25);
xlabel('y')
ylabel('-uv')
title('Reynolds Shear Stress')
legend('DNS','RSM');

%Turbulence Kinetic Energy
figure(6);
plot(y_dns,k_dns,'LineWidth',1.25);
hold on
plot(ynode,k,'-x','LineWidth',1.25);
xlabel('y')
ylabel('-k')
title('Turbulence Kinetic Energy')
legend('DNS','RSM');

%Dissipation Rate of TKE
figure(7);
plot(y_dns,epsilon_dns,'LineWidth',1.25);
hold on
plot(ynode,eps,'-x','LineWidth',1.25);
xlabel('y')
ylabel('epsilon')
title('Dissipation Rate of TKE')
legend('DNS','RSM');

%Comparison with Log Law
figure(8)
yplus_dns = y_dns*ustar_dns / nu;
uplus_loglaw = (1/kappa)*log(yplus_dns) + B;

% Plotting log law from yplus = 50
semilogx(yplus_dns(32:end), uplus_loglaw(32:end), 'LineWidth', 1.25);
hold on

uplus = u/ustar;
semilogx(y_plus, uplus, '-x', 'LineWidth', 1.25)

xlabel('yplus')
ylabel('uplus')
title('Comparison of velocity profile with log law')
legend('Log Law', 'RSM','Location','northwest');
