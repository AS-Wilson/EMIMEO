% Exercise 2 - Figure 2.9

close all
clear all
clc

% Variables
mu0 = 4*pi*1e-7;              %%% Magnetic permeability
eps0 = 8.8542e-12;            %%% Electric permitivity
c = 1/sqrt(mu0*eps0);         %%% Light speed
a = 4.5e-6;                   %%% Fiber core radius
n1 = 1.447;                   %%% Fiber core refractive index
n2 = 1.443;                   %%% Fiber cladding refractive index
V_monomode = 2.045;%0.6756;%2.405;            %%% V for monomode regime
lambda_monomode = 2*pi*a*sqrt((n1^2 - n2^2))/V_monomode;

% Calculation to get the solution of the trascendental equation
num_points = 20000;
X_max = 10; % maximum of the x axis
epsilon = 0.001; %accuracy


x = linspace(0,X_max,num_points);
J_left = -x.*besselj(-1,x)./besselj(0,x);

y = real(sqrt(V_monomode^2-x.^2));
J_right = y.*besselk(-1,y)./besselk(0,y);

for i = 1:length(x)
    distance = abs((J_left(i) - J_right(i)));
    if (J_left(i) < 0 || J_left(i) > J_right(i))
        continue
    elseif distance < epsilon
        min_distance = distance;
        x_sol = x(i);
        index = i;
        kt = x(i)/a;
        gamma = real(sqrt(V_monomode^2-x(i).^2))/a;
        break
    end
 end

%kt = zeros(1,length(v));
%gamma = zeros(1,length(v));
lambda_min = 1e-8;
lambda_max = 4.5e-6;
lambda_points = 10000;
lambda = linspace(lambda_min,lambda_max,lambda_points);

w_freq = 2*pi*c./lambda;
beta0 = w_freq*(n1+n2)/(2*c);
beta_l = (n1*n1*w_freq.*w_freq*eps0*mu0-kt*kt)./(2*beta0) + beta0/2;

n_eff = beta_l./(w_freq*sqrt(eps0*mu0));


%%%%%
% plots
plot(lambda*1e6,n_eff, '-','LineWidth',2)

%xlim([0 5]) %xlim([0 3])
ylim([1.440 1.45]) %ylim([-1 3])
xlabel('$\lambda (\mu m)$','interpreter','latex');
ylabel('$n_{eff}$','interpreter','latex')
title('Refractive index for V = 2.045','interpreter','latex');
xline(lambda_monomode*1e6,'-.','LineWidth',2) % line separating core and cladding

% ESTA MAL. Lo que se hace es calcular para cada lambda, se calcula el V y
% ver si LHS-RHS da cero. Cuando eso ocurra se obtiene el correspondiente V
% el correspondiente X, el correspondiente kt, el correspondiente gamma,
% el correspondiente w, el correspondiente beta0, delta_beta y por ultimo
% el correspondiente beta_l

