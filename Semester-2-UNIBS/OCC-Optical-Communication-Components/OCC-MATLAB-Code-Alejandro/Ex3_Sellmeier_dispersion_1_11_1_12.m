% Refractive index and group velocity calculation from Sellmeier equation
close all
clear all
clc

%% Variables declaration
mu0 = 4*pi*1e-7;                                    %%% Magnetic permeability
eps0 = 8.8542e-12;                                  %%% Electric permitivity
c=sqrt(mu0*eps0)^-1;                                %%% Light speed
Ai_SiO2 = [0.6961663, 0.4079426, 0.8974994];        %%% 100% SiO2 parameters
Ai_GeO2_SiO2 = [0.7136824, 0.4254807, 0.8964226];   %%% 8% GeO2, 92% SiO2 parameters
Lambdai_SiO2 = [0.0684043, 0.1162414, 9.8961609]*1e-6;      %%% 100% SiO2 wavelengths
Lambdai_GeO2_SiO2 = [0.0617167, 0.1270814, 9.8961614]*1e-6; %%% 8% GeO2, 92% SiO2 wavelengths

%%%%%%%%%%%%%%%%%%%%%%%%%%%
Points = 1000;
lambdamin = 1e-6;                                   %%% Minimum wavelenght
lambdamax = 1.6e-6;                                 %%% Maximum wavelenght
nSiO2 = zeros(1,Points);                            %%% NSiO2 empty variable
nGeO2_SiO2 = zeros(1,Points);                       %%% NGeO2_SiO2 empty variable
lambda = linspace(lambdamin,lambdamax,Points);      %%% Wavelenghts under study
omega = 2*pi*c./lambda;                             %%% Wavelenghts under study

%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Refractive index calculation

for i = 1:Points
    
    sum_1 = 0;
    sum_2 = 0;
    for j = 1:length(Ai_SiO2)
        sum_1 = sum_1 + Ai_SiO2(j)*...
            lambda(i)^2/...
            (lambda(i)^2 - Lambdai_SiO2(j)^2);
        
        sum_2 = sum_2 + Ai_GeO2_SiO2(j)*...
            lambda(i)^2/...
            (lambda(i)^2 - Lambdai_GeO2_SiO2(j)^2);
    end
    nSiO2(i) = sqrt(sum_1+1);
    nGeO2_SiO2(i) = sqrt(sum_2+1);

end
beta_SiO2 = 2*pi./lambda.*nSiO2;
beta_GeO2_SiO2 = 2*pi./lambda.*nGeO2_SiO2;

vgSiO2 = diff(omega)./diff(beta_SiO2);
vgGeO2_SiO2 = diff(omega)./diff(beta_GeO2_SiO2);

% k0_2_SiO2 = diff(omega,2)./diff(beta_SiO2,2);
% k0_2_GeO2_SiO2 = diff(omega,2)./diff(beta_GeO2_SiO2,2);
% 
% k0_3_SiO2 = diff(omega,3)./diff(beta_SiO2,3);
% k0_3_GeO2_SiO2 = diff(omega,3)./diff(beta_GeO2_SiO2,3);

%% Plots
%%%%%% Neff pure SiO2 and Ge-doped SiO2 %%%%%%
figure('Renderer', 'painters', 'Position', [400 50 900 700])
hold on
plot(lambda,nSiO2,'-b','LineWidth',2.5)
plot(lambda,nGeO2_SiO2,'-r','LineWidth',2.5)

%%% Other plot settings
title('N Sellmeier(lambda)');
xlabel('Wavelenght(um)');
ylabel('N');
axis([lambdamin lambdamax 1.44 1.465])
legend('100% SiO_2', '8% GeO_2, 92% SiO_2')
ax = gca;
ax.XAxis.Exponent = -6;
set(findall(gcf,'type','text'),'FontSize',30);
set(gca,'FontSize',30);
grid;
hold off;

%%%%%% Group velocity pure SiO2 and Ge-doped SiO2 %%%%%%
figure('Renderer', 'painters', 'Position', [400 50 900 700])
hold on
plot(lambda(2:end),vgSiO2,'-b','LineWidth',2.5)
plot(lambda(2:end),vgGeO2_SiO2,'-r','LineWidth',2.5)

%%% Other plot settings
title('Vgroup(lambda)');
xlabel('Wavelenght(um)');
ylabel('Vg');
axis([lambdamin lambdamax 2.03e8 2.055e8])
legend('100% SiO_2', '8% GeO_2, 92% SiO_2', 'Location', 'east')
ax = gca;
ax.XAxis.Exponent = -6;
ax.YAxis.Exponent = 8;
set(findall(gcf,'type','text'),'FontSize',30);
set(gca,'FontSize',30);
grid;
hold off;

% %%%%%% Second order dispersion pure SiO2 and Ge-doped SiO2 %%%%%%
% figure('Renderer', 'painters', 'Position', [400 50 900 700])
% hold on
% plot(lambda(2:end-1),k0_2_SiO2,'-b','LineWidth',2.5)
% plot(lambda(2:end-1),k0_2_GeO2_SiO2,'-r','LineWidth',2.5)
% 
% %%% Other plot settings
% title('Second order Dispersion');
% xlabel('Wavelenght(um)');
% ylabel('K_0´´');
% legend('100% SiO_2', '8% GeO_2, 92% SiO_2')
% ax = gca;
% set(findall(gcf,'type','text'),'FontSize',30);
% set(gca,'FontSize',30);
% grid;
% hold off;
% 
% %%%%%% Third´ order dispersion pure SiO2 and Ge-doped SiO2 %%%%%%
% figure('Renderer', 'painters', 'Position', [400 50 900 700])
% hold on
% plot(lambda(2:end-2),k0_3_SiO2,'-b','LineWidth',2.5)
% plot(lambda(2:end-2),k0_3_GeO2_SiO2,'-r','LineWidth',2.5)
% 
% %%% Other plot settings
% title('Third order Dispersion');
% xlabel('Wavelenght(um)');
% ylabel('K_0´´´');
% legend('100% SiO_2', '8% GeO_2, 92% SiO_2')
% ax = gca;
% set(findall(gcf,'type','text'),'FontSize',30);
% set(gca,'FontSize',30);
% grid;
% hold off;