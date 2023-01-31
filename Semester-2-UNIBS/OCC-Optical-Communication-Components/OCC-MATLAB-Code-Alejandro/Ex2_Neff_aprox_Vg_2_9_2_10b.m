% Refractive index and group velocity approximation for LP01 mode 
close all
clear all
clc

%% Variables declaration
mu0 = 4*pi*1e-7;              %%% Magnetic permeability
eps0 = 8.8542e-12;            %%% Electric permitivity
c = 1/sqrt(mu0*eps0);         %%% Light speed
a = 4.5e-6;                   %%% Fiber core radius
n1 = 1.447;                   %%% Fiber core refractive index
n2 = 1.443;                   %%% Fiber cladding refractive index
Vmonomode = 2.405;            %%% V for monomode regime
lambda_monomode = 2*pi*a*sqrt((n1^2 - n2^2))/Vmonomode;

%%%%%%%%%%%%%%%%%%%%%%%%%%%
Points_X = 20000;
Points_L = 60; %110;
l0 = 0;                                     %%% Bessel function order l
lambdamax = 4.5e-6;                         %%% Maximum wavelenght
neff = zeros(1,Points_L);                     %%% Neff empty variable
beta = zeros(1,Points_L);                     %%% Beta'' empty variable
neff_approx = zeros(1,Points_L);              %%% Neff approx empty variable
lambda = linspace(1e-8,lambdamax,Points_L);   %%% Wavelenghts under study

%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i = 1:Points_L
    %% RHS and LHS calculation
    V = 2*pi*a*sqrt((n1^2 - n2^2))/lambda(i);   %%% V calculation
    X = linspace(0,V,Points_X);                   %%% Fixed variable to act as x axis X=kt*a              
    Y = real(sqrt(V.^2-X.^2));                  %%% Relationship X,Y
    LHS0 = -X.*besselj(l0-1,X)./besselj(l0,X);  %%% Left hand side of dispersion relation 
    RHS0 = Y.*besselk(l0-1,Y)./besselk(l0,Y);   %%% Right hand side of dispersion relation 

    %% Find LP01 solution
    minDistance = inf;
    distances = zeros(1,length(LHS0));
    for ni = 1 : length(LHS0)
        distances(ni) = abs(LHS0(ni) - RHS0(ni));
        if distances(ni) < minDistance
          minXindex = ni;
          minDistance = distances(ni);
        else
          break
        end
    end
  
    %% Propagation parameters calculation
    XLP01 = X(minXindex);
    YLP01 = Y(minXindex);
    wLP01 = sqrt((XLP01^2 + YLP01^2)/(mu0*eps0*(n1^2-n2^2)*a^2));
    ktLP01 = XLP01/a;
    gammaLP01 = YLP01/a;
    beta0 = (wLP01/c*n1 + wLP01/c*n2)/2;
    delta_beta = (n1^2*wLP01^2*mu0*eps0-beta0^2-ktLP01^2)/(2*beta0);
    beta(i) = beta0+delta_beta;
    neff(i) = beta(i)/(wLP01*sqrt(mu0*eps0));

    %% Effective refractive index approximation calculation
    neff_approx(i) = n2+(n1-n2)*(1.1428-0.996/V)^2;
end

window = ((lambda >= 1e-6)&(lambda <= 1.6e-6));

lambda_vg = nonzeros(lambda.*window)';
omega_vg = 2*pi*c./lambda_vg;

beta_vg = nonzeros(beta.*window)';

vg = diff(omega_vg)./diff(beta_vg);

%% Plots
%%%%%% Neff real and approximation %%%%%%
figure('Renderer', 'painters', 'Position', [400 50 900 700])
hold on
plot(lambda,neff,'-g','LineWidth',2.5)
plot(lambda,neff_approx,'--r','LineWidth',2.5)
xline(lambda_monomode,'-.','LineWidth',2);

%%% Draw arrow
x = [0.4 0.7];
y = [0.7 0.7];
annotation('textarrow',x,y)

%%% Insert text
xclad = 1.5e-6;
yclad = 1.446;
str = 'Monomodal region';
text(xclad,yclad,str)

%%% Other plot settings
title('Neff(lambda)');
xlabel('Wavelenght(m)');
ylabel('Neff');
axis([0 lambdamax 1.443 1.447])
legend('Real', 'Gaussian')
ax = gca;
ax.XAxis.Exponent = -6;
set(findall(gcf,'type','text'),'FontSize',30);
set(gca,'FontSize',30);
grid;

%%%%%% Group velocity %%%%%%
figure('Renderer', 'painters', 'Position', [400 50 900 700])
hold on
plot(lambda_vg(1:end-1),vg,'-g','LineWidth',2.5)
xline(lambda_monomode,'-.','LineWidth',2);

%%% Draw arrow
x = [0.53 0.83];
y = [0.7 0.7];
annotation('textarrow',x,y)

%%% Insert text
xclad = 1.3e-6;
yclad = 2.0714e8;
str = 'Monomodal region';
text(xclad,yclad,str)

%%% Other plot settings
title('Vg(lambda)');
xlabel('Wavelenght(m)');
ylabel('Vg');
legend('Real', 'Gaussian')
axis([lambda_vg(1) lambda_vg(end) 2.0708*1e8 2.0716*1e8])
ax = gca;
ax.XAxis.Exponent = -6;
ax.YAxis.Exponent = 8;
set(findall(gcf,'type','text'),'FontSize',30);
set(gca,'FontSize',30);
grid;
