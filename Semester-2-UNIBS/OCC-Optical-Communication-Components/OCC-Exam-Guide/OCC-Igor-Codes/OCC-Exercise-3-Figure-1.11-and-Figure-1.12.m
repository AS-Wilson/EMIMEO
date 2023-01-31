% Exercise 3 - Figure 2.11 and 2.12

close all
clear all
clc

% Variables
mu0 = 4*pi*1e-7;              %%% Magnetic permeability
eps0 = 8.8542e-12;            %%% Electric permitivity
c = 1/sqrt(mu0*eps0);         %%% Light speed
Ai_SiO2 = [0.6961663, 0.4079426, 0.8974994];        %%% 100% SiO2 parameters
Ai_GeO2_SiO2 = [0.7136824, 0.4254807, 0.8964226];   %%% 8% GeO2, 92% SiO2 parameters
Lambdai_SiO2 = [0.0684043, 0.1162414, 9.8961609]*1e-6;      %%% 100% SiO2 wavelengths
Lambdai_GeO2_SiO2 = [0.0617167, 0.1270814, 9.8961614]*1e-6; %%% 8% GeO2, 92% SiO2 wavelengths

% Vectors for calculations
lambda_in = 1e-6;
lambda_fin = 1.8e-6;
points = 1000;
lambda = linspace(lambda_in, lambda_fin, points);
omega = 2*pi*c./lambda;

n_SiO2 = zeros(1,points);
n_GeO2_SiO2 = zeros(1,points);


% Calculation of n
for i=1:points
    sum_SiO2 = 0;
    sum_GeO2_SiO2 = 0;
    % Calculating the sumatory
    for j=1:3
        sum_SiO2 = sum_SiO2 + (Ai_SiO2(j)*lambda(i)^2)/(lambda(i)^2-Lambdai_SiO2(j)^2);
        sum_GeO2_SiO2 = sum_GeO2_SiO2 + (Ai_GeO2_SiO2(j)*lambda(i)^2)/(lambda(i)^2-Lambdai_GeO2_SiO2(j)^2);
    end

    n_SiO2(i) = sqrt(sum_SiO2 + 1);
    n_GeO2_SiO2(i) = sqrt(sum_GeO2_SiO2 + 1);
end

% Calculation of betas (prop. constant)
beta_SiO2 = 2*pi./lambda.*n_SiO2;
beta_GeO2_SiO2 = 2*pi./lambda.*n_GeO2_SiO2;

% Calculation of group velocities
Vg_SiO2 = diff(omega)./diff(beta_SiO2);
Vg_GeO2_SiO2 = diff(omega)./diff(beta_GeO2_SiO2);


% Plot of n(lambda)
plot(lambda*1e6,n_SiO2,'-',lambda*1e6,n_GeO2_SiO2,'-','LineWidth',2)
xlabel('$\lambda$($\mu$m)','interpreter','latex');
ylabel('n','interpreter','latex')
title('Refractive index with Sellmeier formula','interpreter','latex');
legend('100% SiO2','8% GeO2')

% Plot of Vg(lambda)
plot(lambda(2:end)*1e6,Vg_SiO2*1e-8,'-',lambda(2:end)*1e6,Vg_GeO2_SiO2*1e-8,'-','LineWidth',2)
xlabel('$\lambda$($\mu$m)','interpreter','latex');
ylabel('$V_{g}$ [$10^{8}$ m/s]','interpreter','latex')
title('Refractive index with Sellmeier formula','interpreter','latex');
legend('100% SiO2','8% GeO2')




