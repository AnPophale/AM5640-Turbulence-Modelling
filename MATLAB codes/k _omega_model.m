clc; clear all; close all; format shortg; format compact;

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
beta_star = 0.09;
C1 = 5/9;
C2 = 3/40;
sigma_k = 2;
sigma_omega = 2;
kappa = 0.41;

% Residual Error Limit
residual_limit = 10^(-4);
residual = 1;

% Under Relaxation Factor
urf = 0.5;

%Initial Guess and BCs
u = ones(N+2,1); u(1) = 0; %Wall BC
k = ones(N+2,1); k(1) = 0; %Wall BC

omega = ones(N+2,1);
omega(2:8) = (1/C2)*6*nu./(ynode(2:8).^2); %Wall BC
omega(1) = omega(2);

nu_t = zeros(N+2,1);
nu_t_n = zeros(N,1); nu_t_s = zeros(N,1);

aP_u = zeros(N,1); aN_u = zeros(N,1); aS_u = zeros(N,1); Su_u = zeros(N,1);
aP_k = zeros(N,1); aN_k = zeros(N,1); aS_k = zeros(N,1); Su_k = zeros(N,1); Sp_k = zeros(N,1);
aP_omega = zeros(N,1); aN_omega = zeros(N,1); aS_omega = zeros(N,1); Su_omega = zeros(N,1); Sp_omega = zeros(N,1);

%Main interation loop
iter = 1;
while residual > residual_limit
    residual_u = 0;
    residual_k = 0;
    residual_omega = 0;

    u_old = u;
    k_old = k;
    omega_old = omega;

    %Calculating Turbulence Viscosity at Nodes
    for i= 1:N+2
        nu_t(i) = k_old(i)/omega_old(i);
    end

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
    
    %Solving for u using Gauss Seidel Method
    for i = 1:N
        u(i+1) = (aN_u(i)*u(i+2) + aS_u(i)*u(i) + Su_u(i))/aP_u(i);
    end
    u(end) = u(end-1); %Neumann BC at channel centerline
    
    %Under Relaxation
    u = (1-urf)*u + urf*u_old;

    %Residual for u
    for i = 1:N
        residual_u = residual_u + abs(aN_u(i)*u(i+2) + aS_u(i)*u(i) + Su_u(i) - aP_u(i)*u(i+1));
    end
    residual_u_vec(iter) = residual_u;

    %Calculating aP, aN, aS, Su, Sp for k equation
    for i = 1:N
        aN_k(i) = (nu + (nu_t_n(i)/sigma_k))/dy_n(i);
        aS_k(i) = (nu + (nu_t_s(i)/sigma_k))/dy_s(i);
        aP_k(i) = aN_k(i) + aS_k(i);

        Pk = nu_t(i+1)*(((u_old(i+2)-u_old(i))/(dy_n(i) + dy_s(i)))^2);
        Su_k(i) = Pk*dy(i);

        Sp_k(i) = beta_star*omega_old(i+1)*dy(i);
    end

    %Solving for k using Gauss Seidel Method
    for i = 1:N
        k(i+1) = (aN_k(i)*k(i+2) + aS_k(i)*k(i) + Su_k(i))/(aP_k(i) + Sp_k(i));
    end
    k(end) = k(end-1); %Neumann BC at channel centerline

    %Under Relaxation
    k = (1-urf)*k + urf*k_old;

    %Residual for k
    for i = 1:N
        residual_k = residual_k + abs(aN_k(i)*k(i+2) + aS_k(i)*k(i) + Su_k(i) - (aP_k(i) + Sp_k(i))*k(i+1));
    end
    residual_k_vec(iter) = residual_k;

    %Calculating aP, aN, aS, Su, Sp for omega equation
    for i = 1:N
        aN_omega(i) = (nu + (nu_t_n(i)/sigma_omega))/dy_n(i);
        aS_omega(i) = (nu + (nu_t_s(i)/sigma_omega))/dy_s(i);
        aP_omega(i) = aN_omega(i) + aS_omega(i);

        Pk = nu_t(i+1)*(((u_old(i+2)-u_old(i))/(dy_n(i) + dy_s(i)))^2);
        Su_omega(i) = C1*Pk*omega_old(i+1)*dy(i)/k_old(i+1);

        Sp_omega(i) = C2*dy(i)*omega_old(i+1);
    end

    %Solving for omega using Gauss Seidel Method 
    %(Omega is not solved for first 8 nodes, it comes from the boundary condition) 
    for i = 8:N
        omega(i+1) = (aN_omega(i)*omega(i+2) + aS_omega(i)*omega(i) + Su_omega(i))/(aP_omega(i) + Sp_omega(i));
    end
    omega(end) = omega(end-1); %Neumann BC at channel centerline
    omega(2:8) = (1/C2)*6*nu./(ynode(2:8).^2); %Wall BC
    omega(1) = omega(2);

    %Under Relaxation
    omega = (1-urf)*omega + urf*omega_old;

    %Residual for omega
    for i = 8:N
        residual_omega = residual_omega + abs(aN_omega(i)*omega(i+2) + aS_omega(i)*omega(i) + Su_omega(i) - (aP_omega(i) + Sp_omega(i))*omega(i+1));
    end
    residual_omega_vec(iter) = residual_omega;

    %Residual Calculation
    residual = max([residual_u, residual_k, residual_omega]);
    iter = iter+1;
end

%Calculating epsilon from k and omega
epsilon = beta_star*k.*omega;

%Calculating Wall Shear Stress
tau_wall = nu*rho*(u(2)-u(1))/dy_s(1)

%Plotting the residual
figure(1)
semilogy(1:1:iter-1,residual_u_vec,'-r','LineWidth',1.25)
hold on 
semilogy(1:1:iter-1,residual_k_vec,'-k','LineWidth',1.25)
hold on
semilogy(1:1:iter-1,residual_omega_vec,'-b','LineWidth',1.25)
xlabel('Iteration Number')
ylabel('Residual')
legend('u','k','omega')
title('Residual Plot')

%Plotting Results
%Loading DNS Data
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
figure(2)
semilogx(y_dns/nu,u_dns,'-r','LineWidth',1.25);
hold on
semilogx(ynode/nu,u,'-k','LineWidth',1.25)
xlabel('y plus'); 
ylabel('U'); 
title('U-velocity');
legend('DNS','k-omega Model','Location','northwest')

%Turbulent Kinetic Energy
figure(3)
plot(y_plus,k_dns,'-r','LineWidth',1.25);
hold on
plot(y_plus,k,'-k','LineWidth',1.25);
xlabel('y plus'); 
ylabel('k'); 
title('Turbulence kinetic energy');
legend('DNS','k-omega Model')

%Dissipation rate of TKE
figure(4)
plot(y_plus,epsilon_dns,'-r','LineWidth',1.25)
hold on
plot(y_plus,epsilon,'-k','LineWidth',1.25)
xlabel('y plus'); 
ylabel('Epsilon'); 
title('Dissipation rate of k');
legend('DNS','k-omega Model');

%Reynolds Shear Stress
%Calculating nu_t*dk/dy
nu_tdu_dy = zeros(N+2,1);
for i = 1:N
    nu_tdu_dy(i+1) = nu_t(i+1)*(u(i+2)-u(i))/(dy_n(i)+dy_s(i));
end
nu_tdu_dy(1) = nu_t(1)*(u(2)-u(1))/dy_s(1);
nu_tdu_dy(N+2) = nu_t(N+2)*(u(N+2)-u(N+1))/dy_n(end);

%Plotting
uv_dns = load("uv_dns.dat");
figure(5)
plot(y_plus,nu_tdu_dy,'-r','LineWidth',1.25)
hold on
plot(y_plus,-uv_dns,'-k','LineWidth',1.25)
xlabel('y plus')
ylabel('-uv')
title('Reynolds Shear Stress')
legend('DNS','k-omega Model');

%Turbulence Viscosity
figure(6)
plot(y_plus,nu_t,'-r','LineWidth',1.25)
xlabel('y plus')
ylabel('nu_t')
title('Turbulence Viscosity')
legend('k-omega Model','Location','northwest')

%Production Rate of TKE
%Calculating Pk model
Pk = zeros(N+2,1);
for i = 1:N
    Pk(i+1) = nu_t(i+1)*(((u_old(i+2)-u_old(i))/(dy_n(i) + dy_s(i)))^2);
end
Pk(1) = 0; %Turbulence viscosity on wall is zero
Pk(end) = 0; %No velocity gradient at centerline

%Plotting
figure(7)
plot(y_plus,Pk_dns,'-r','LineWidth',1.25)
hold on
plot(y_plus,Pk,'-k','LineWidth',1.25)
xlabel('y plus')
ylabel('P_k')
title('Production rate of k')
legend('DNS','k-omega Model');

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
figure(8)
plot(y_plus,Dk_viscous_dns,'-r','LineWidth',1.25)
hold on
plot(y_plus,Dk_viscous,'-k','LineWidth',1.25)
xlabel('y plus')
ylabel('D_k Viscous')
title('Viscous diffusion rate of k')
legend('DNS','k-omega Model');

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
figure(9)
plot(y_plus,Dk_turbulent_dns,'-r','LineWidth',1.25)
hold on
plot(y_plus,Dk_turbulent,'-k','LineWidth',1.25)
xlabel('y plus')
ylabel('D_k Turbulent')
title('Turbulent diffusion rate of k')
legend('DNS','k-omega Model');

%Pressure Diffusion of TKE
figure(10)
plot(y_plus,zeros(N+2,1),'-r','LineWidth',1.25) %Pressure diffusion neglected in k omega model
hold on
plot(y_plus,Dk_pressure_dns,'-k','LineWidth',1.25)
xlabel('y plus')
ylabel('D_k Pressure')
title('Pressure diffusion rate of k')
legend('DNS','k-omega Model');

%Budget of TKE (DNS)
figure(11)
hold on
plot(y_plus,Pk_dns,'LineWidth',1.25)
plot(y_plus,epsilon_dns,'LineWidth',1.25)
plot(y_plus,Dk_viscous_dns,'LineWidth',1.25)
plot(y_plus,Dk_turbulent_dns,'LineWidth',1.25)
plot(y_plus,Dk_pressure_dns,'LineWidth',1.25)
xlabel('yplus')
title('Budget of TKE (DNS)')
legend('P_k','Epsilon','D_k Viscous','D_k Turbulent','D_k Pressure');

%Budget of TKE (Model)
Dk_pressure = zeros(N+2,1); %Pressure diffusion neglected in k omega model
figure(12)
hold on
plot(y_plus,Pk,'LineWidth',1.25)
plot(y_plus,epsilon,'LineWidth',1.25)
plot(y_plus,Dk_viscous,'LineWidth',1.25)
plot(y_plus,Dk_turbulent,'LineWidth',1.25)
plot(y_plus,Dk_pressure,'LineWidth',1.25)
xlabel('yplus')
title('Budget of TKE (k-omega model)')
legend('P_k','Epsilon','D_k Viscous','D_k Turbulent','D_k Pressure');