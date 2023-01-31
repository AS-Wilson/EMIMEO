% Effective refractive index variation with respect to normalised frequency
close all
clear all
clc

%% Variables declaration
mu0 = 4*pi*1e-7;                                    %%% Magnetic permeability
eps0 = 8.8542e-12;                                  %%% Electric permitivity
c = 1/sqrt(mu0*eps0);                               %%% Light speed
a = 4.5e-6;                                         %%% Fiber core radius
l0=0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%
Points = 1000;                                          %%% Bessel function order l
neff = zeros(1,Points);                          %%% Neff empty variable
nSiO2 = 1.443;                            %%% NSiO2 empty variable
nGeO2_SiO2 = 1.447;                       %%% NGeO2_SiO2 empty variable
Vmax = 4;
V=linspace(0.01,Vmax,Points);
%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i = 1:Points
    %% RHS and LHS calculation
    X = linspace(0,V(i),Points);                                   %%% Fixed variable to act as x axis X=kt*a              
    Y = real(sqrt(V(i).^2-X.^2));                                  %%% Relationship X,Y
    LHS0 = -X.*besselj(l0-1,X)./besselj(l0,X);                  %%% Left hand side of dispersion relation 
    RHS0 = Y.*besselk(l0-1,Y)./besselk(l0,Y);                   %%% Right hand side of dispersion relation 

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
    wLP01 = sqrt((XLP01^2 + YLP01^2)/(mu0*eps0*(nGeO2_SiO2^2-nSiO2^2)*a^2));
    ktLP01 = XLP01/a;
    gammaLP01 = YLP01/a;
    beta0 = (wLP01/c*nGeO2_SiO2 + wLP01/c*nSiO2)/2;
    delta_beta = (nGeO2_SiO2^2*wLP01^2*mu0*eps0-beta0^2-ktLP01^2)/(2*beta0);
    beta = beta0+delta_beta;
    neff(i) = beta/(wLP01*sqrt(mu0*eps0));

end

%% Plots
%%%%%% Neff with respect to lambda %%%%%%
figure('Renderer', 'painters', 'Position', [400 50 900 700])
hold on
plot(V,neff,'-g','LineWidth',2.5)

%%% Other plot settings
title('Neff LP01 vs V');
xlabel('V');
ylabel('Neff');
set(findall(gcf,'type','text'),'FontSize',30);
set(gca,'FontSize',30);
grid;
