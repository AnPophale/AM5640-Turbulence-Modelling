clc; clear all; close all;format shortg; format compact;

%Node locations
ynode = load("y_dns.dat");
N = numel(ynode) - 2; %Number of cells

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
    
%Physcial Properties
nu = 1/395;
rho = 1;
u_star = 1;

%Yplus Values
y_plus = ynode/(nu/u_star);

%Model Constants
C_mu = 0.09;
C1 = 1.44;
C2 = 1.92;
sigma_k = 1;
sigma_epsilon = 1.3;
kappa = 0.41;

% Residual Error Limit
residual_limit = 10^(-4);
residual = 1;

% Under Relaxation Factor
urf = 0.5;

%Initial Guess and BCs
u = ones(N+2,1); u(1) = 0; %Wall BC
k = ones(N+2,1); k(1) = 0; %Wall BC

epsilon = ones(N+2,1);
epsilon(2:8) = 2*nu*k(2:8)./(ynode(2:8).^2); %Wall BC
epsilon(1) = epsilon(2);

nu_t = zeros(N+2,1);
nu_t_n = zeros(N,1); nu_t_s = zeros(N,1);

aP_u = zeros(N,1); aN_u = zeros(N,1); aS_u = zeros(N,1); Su_u = zeros(N,1);
aP_k = zeros(N,1); aN_k = zeros(N,1); aS_k = zeros(N,1); Su_k = zeros(N,1); Sp_k = zeros(N,1);
aP_e = zeros(N,1); aN_e = zeros(N,1); aS_e = zeros(N,1); Su_e = zeros(N,1); Sp_e = zeros(N,1);

%Main interation loop
iter = 0;
while residual > residual_limit
    u_old = u;
    k_old = k;
    epsilon_old = epsilon;

    %Calculating Turbulence Viscosity at Nodes
    for i= 2:N+2
        nu_t(i) = C_mu*(k_old(i)^2)/epsilon_old(i);
    end
    nu_t(1) = 0; %Wall BC

    %Interpolating Turbulence Viscosity to Cell Faces
    for i = 1:N
        nu_t_n(i) = fn(i)*nu_t(i+2) + (1-fn(i))*nu_t(i+1);
        nu_t_s(i) = fs(i)*nu_t(i) + (1-fs(i))*nu_t(i+1);
    end
    nu_t_n(end) = nu_t(end);
    nu_t_s(1) = nu_t(1);
    
    %Calculating aP, aN, aS, Su for u equation
    for i = 1:N
        aN_u(i) = (nu + nu_t_n(i))/dy_n(i);
        aS_u(i) = (nu + nu_t_s(i))/dy_s(i);
        aP_u(i) = aN_u(i) + aS_u(i);
        Su_u(i) = dy(i)/rho;
    end
    
    %Solving for u
    for i = 1:N
        u(i+1) = (aN_u(i)*u(i+2) + aS_u(i)*u(i) + Su_u(i))/aP_u(i);
    end
    u(end) = u(end-1); %Neumann BC at channel centerline

    %Under Relaxation
    u = (1-urf)*u + urf*u_old;

    %Calculating aP, aN, aS, Su, Sp for k equation
    for i = 1:N
        aN_k(i) = (nu + (nu_t_n(i)/sigma_k))/dy_n(i);
        aS_k(i) = (nu + (nu_t_s(i)/sigma_k))/dy_s(i);
        aP_k(i) = aN_k(i) + aS_k(i);

        Pk = nu_t(i+1)*(((u_old(i+2)-u_old(i))/(dy_n(i) + dy_s(i)))^2); %Central difference on non uniform grid?
        Su_k(i) = Pk*dy(i);

        Sp_k(i) = epsilon_old(i+1)*dy(i)/k_old(i+1);
    end

    %Solving for k
    for i = 1:N
        k(i+1) = (aN_k(i)*k(i+2) + aS_k(i)*k(i) + Su_k(i))/(aP_k(i) + Sp_k(i));
    end
    k(end) = k(end-1); %Neumann BC at channel centerline

    %Under Relaxation
    k = (1-urf)*k + urf*k_old;

    %Calculating aP, aN, aS, Su, Sp for epsilon equation
    for i = 1:N
        aN_e(i) = (nu + (nu_t_n(i)/sigma_epsilon))/dy_n(i);
        aS_e(i) = (nu + (nu_t_s(i)/sigma_epsilon))/dy_s(i);
        aP_e(i) = aN_e(i) + aS_e(i);

        Pk = nu_t(i+1)*(((u_old(i+2)-u_old(i))/(dy_n(i) + dy_s(i)))^2);
        Su_e(i) = C1*Pk*epsilon_old(i+1)*dy(i)/k_old(i+1);

        Sp_e(i) = C2*epsilon_old(i+1)*dy(i)/k_old(i+1);
    end

    %Solving for epsilon
    for i = 1:N
        epsilon(i+1) = (aN_e(i)*epsilon(i+2) + aS_e(i)*epsilon(i) + Su_e(i))/(aP_e(i) + Sp_e(i));
    end
    epsilon(end) = epsilon(end-1); %Neumann BC at channel centerline
    epsilon(2:8) = 2*nu*k(2:8)./(ynode(2:8).^2); %Wall BC
    epsilon(1) = epsilon(2);

    %Under Relaxation
    epsilon = (1-urf)*epsilon + urf*epsilon_old;

    %Residual Calculation
    residual = max([norm(u - u_old), norm(k-k_old), norm(epsilon-epsilon_old)]);
    iter = iter+1;
end

%Plotting Results
%DNS Data
load dns_data.dat
load y_dns.dat
load u_dns.dat
load u2_dns.dat
load v2_dns.dat
load w2_dns.dat

%Calculating TKE from dns data
k_dns=0.5*(u2_dns+v2_dns+w2_dns);

%All terms in TKE equation are normalized by ustar^4/nu
epsilon_dns=dns_data(:,2)*u_star^4/nu; 
Pk_dns = dns_data(:,3)*u_star^4/nu;
Dk_pressure_dns = dns_data(:,4)*u_star^4/nu;
Dk_turbulent_dns = dns_data(:,5)*u_star^4/nu;
Dk_viscous_dns = dns_data(:,6)*u_star^4/nu;

%U Velocity
figure(1)
semilogx(y_dns/nu,u_dns,'-r','LineWidth',1.25);
hold on
semilogx(ynode/nu,u,'-b','LineWidth',1.25)
xlabel('y plus'); 
ylabel('U'); 
title('U-velocity');
legend('DNS','K-epsilon Model','Location','northwest')

%Turbulent Kinetic Energy
figure(2)
plot(y_plus,k_dns,'-r','LineWidth',1.25);
hold on
plot(y_plus,k,'-b','LineWidth',1.25);
xlabel('y plus'); 
ylabel('k'); 
title('Turbulence kinetic energy');
legend('DNS','K-epsilon Model')

%Dissipation rate of TKE
figure(3)
plot(y_plus,epsilon_dns,'-r','LineWidth',1.25)
hold on
plot(y_plus,epsilon,'-b','LineWidth',1.25)
xlabel('y plus'); 
ylabel('Epsilon'); 
title('Dissipation rate of k');
legend('DNS','K-epsilon Model');

%Production Rate of TKE
%Calculating Pk model
Pk = zeros(N+2,1);
for i = 1:N
    Pk(i+1) = nu_t(i+1)*(((u_old(i+2)-u_old(i))/(dy_n(i) + dy_s(i)))^2);
end
Pk(1) = 0; %Turbulence viscosity on wall is zero
Pk(end) = 0; %No velocity gradient at centerline

%Plotting
figure(4)
plot(y_plus,Pk_dns,'-r','LineWidth',1.25)
hold on
plot(y_plus,Pk,'-b','LineWidth',1.25)
xlabel('y plus')
ylabel('P_k')
title('Production rate of k')
legend('DNS','K-epsilon Model');

%Viscous diffusion of TKE
%Calculating nu*dk/dy
nudk_dy = zeros(N+2,1);
for i = 1:N
    nudk_dy(i+1) = nu*(k(i+2)-k(i))/(dy_n(i)+dy_s(i));
end
nudk_dy(1) = nu*(k(2)-k(1))/dy_s(1);
nudk_dy(N+2) = nu*(k(N+2)-k(N+1))/dy_n(end);

%Calculating d/dy(nu*dk/dy)) 
Dk_viscous = zeros(N+2,1);
for i = 1:N
    Dk_viscous(i+1) = (nudk_dy(i+2)-nudk_dy(i))/(dy_n(i)+dy_s(i));
end
Dk_viscous(1) = (nudk_dy(2)-nudk_dy(1))/dy_s(1);
Dk_viscous(N+2) = (nudk_dy(N+2)-nudk_dy(N+1))/dy_n(end);

%Plotting
figure(5)
plot(y_plus,Dk_viscous_dns,'-r','LineWidth',1.25)
hold on
plot(y_plus,Dk_viscous,'-b','LineWidth',1.25)
xlabel('y plus')
ylabel('D_k Viscous')
title('Viscous diffusion rate of k')
legend('DNS','K-epsilon Model');

%Turbulent diffusion of TKE
%Calculating nu_t*dk/dy
nu_tdk_dy = zeros(N+2,1);
for i = 1:N
    nu_tdk_dy(i+1) = nu_t(i+1)*(k(i+2)-k(i))/(dy_n(i)+dy_s(i));
end
nu_tdk_dy(1) = nu_t(1)*(k(2)-k(1))/dy_s(1);
nu_tdk_dy(N+2) = nu_t(N+2)*(k(N+2)-k(N+1))/dy_n(end);

%Calculating d/dy(nu_t*dk/dy)) 
Dk_turbulent = zeros(N+2,1);
for i = 1:N
    Dk_turbulent(i+1) = (nu_tdk_dy(i+2)-nu_tdk_dy(i))/(dy_n(i)+dy_s(i));
end
Dk_turbulent(1) = (nu_tdk_dy(2)-nu_tdk_dy(1))/dy_s(1);
Dk_turbulent(N+2) = (nu_tdk_dy(N+2)-nu_tdk_dy(N+1))/dy_n(end);

%Plotting
figure(6)
plot(y_plus,Dk_turbulent_dns,'-r','LineWidth',1.25)
hold on
plot(y_plus,Dk_turbulent,'-b','LineWidth',1.25)
xlabel('y plus')
ylabel('D_k Turbulent')
title('Turbulent diffusion rate of k')
legend('DNS','K-epsilon Model');

%Budget of TKE (DNS)
figure(7)
hold on
plot(y_plus,Pk_dns,'LineWidth',1.25)
plot(y_plus,epsilon_dns,'LineWidth',1.25)
plot(y_plus,Dk_turbulent_dns,'LineWidth',1.25)
plot(y_plus,Dk_viscous_dns,'LineWidth',1.25)
plot(y_plus,Dk_pressure_dns,'LineWidth',1.25)
xlabel('yplus')
title('Budget of TKE (DNS)')
legend('P_k','Epsilon','D_k Viscous','D_k Turbulent','D_k Pressure');

%Budget of TKE (Model)
Dk_pressure = zeros(N+2,1); %Pressure diffusion neglected in k epsilon model
figure(8)
hold on
plot(y_plus,Pk,'LineWidth',1.25)
plot(y_plus,epsilon,'LineWidth',1.25)
plot(y_plus,Dk_turbulent,'LineWidth',1.25)
plot(y_plus,Dk_viscous,'LineWidth',1.25)
plot(y_plus,Dk_pressure,'LineWidth',1.25)
xlabel('yplus')
title('Budget of TKE (k-epsilon model)')
legend('P_k','Epsilon','D_k Viscous','D_k Turbulent','D_k Pressure');