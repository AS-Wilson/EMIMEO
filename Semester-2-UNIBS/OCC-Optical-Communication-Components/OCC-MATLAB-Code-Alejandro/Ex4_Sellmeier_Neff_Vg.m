% Refractive index and group velocity NSellmeier for LP01 mode 
close all
clear all
clc

%% Variables declaration
mu0 = 4*pi*1e-7;                                    %%% Magnetic permeability
eps0 = 8.8542e-12;                                  %%% Electric permitivity
c = 1/sqrt(mu0*eps0);                               %%% Light speed
a = 4.5e-6;                                         %%% Fiber core radius
Ai_SiO2 = [0.6961663, 0.4079426, 0.8974994];        %%% 100% SiO2 parameters
Ai_GeO2_SiO2 = [0.7136824, 0.4254807, 0.8964226];   %%% 8% GeO2, 92% SiO2 parameters
Lambdai_SiO2 = [0.0684043, 0.1162414, 9.8961609]*1e-6;      %%% 100% SiO2 wavelengths
Lambdai_GeO2_SiO2 = [0.0617167, 0.1270814, 9.8961614]*1e-6; %%% 8% GeO2, 92% SiO2 wavelengths


%%%%%%%%%%%%%%%%%%%%%%%%%%%
Points_X = 25000;
Points_L = 40;
l0 = 0;                                          %%% Bessel function order l
lambdamin = 1e-6;                                %%% Minimum wavelenght
lambdamax = 1.6e-6;                              %%% Maximum wavelenght
neff = zeros(1,Points_L);                          %%% Neff empty variable
beta = zeros(1,Points_L);                          %%% Beta'' empty variable
lambda = linspace(lambdamin,lambdamax,Points_L);   %%% Wavelenghts under study
nSiO2 = zeros(1,Points_L);                            %%% NSiO2 empty variable
nGeO2_SiO2 = zeros(1,Points_L);                       %%% NGeO2_SiO2 empty variable

%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i = 1:Points_L
    %% N1 and N2 calculations
    for j = 1:Points_L
        sum_1 = 0;
        sum_2 = 0;
        for k = 1:length(Ai_SiO2)
            sum_1 = sum_1 + Ai_SiO2(k)*...
                lambda(j)^2/...
                (lambda(j)^2 - Lambdai_SiO2(k)^2);

            sum_2 = sum_2 + Ai_GeO2_SiO2(k)*...
                lambda(j)^2/...
                (lambda(j)^2 - Lambdai_GeO2_SiO2(k)^2);
        end
        nSiO2(j) = sqrt(sum_1+1);
        nGeO2_SiO2(j) = sqrt(sum_2+1);
    end
    
    %% RHS and LHS calculation
    V = 2*pi*a*sqrt((nGeO2_SiO2(i)^2 - nSiO2(i)^2))/lambda(i);  %%% V calculation
    X = linspace(0,V,Points_X);                                   %%% Fixed variable to act as x axis X=kt*a              
    Y = real(sqrt(V.^2-X.^2));                                  %%% Relationship X,Y
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
    wLP01 = sqrt((XLP01^2 + YLP01^2)/(mu0*eps0*(nGeO2_SiO2(i)^2-nSiO2(i)^2)*a^2));
    ktLP01 = XLP01/a;
    gammaLP01 = YLP01/a;
    beta0 = (wLP01/c*nGeO2_SiO2(i) + wLP01/c*nSiO2(i))/2;
    delta_beta = (nGeO2_SiO2(i)^2*wLP01^2*mu0*eps0-beta0^2-ktLP01^2)/(2*beta0);
    beta(i) = beta0+delta_beta;
    neff(i) = beta(i)/(wLP01*sqrt(mu0*eps0));

end

omega = 2*pi*c./lambda;
vg = diff(omega)./diff(beta);

%% Plots
%%%%%% Neff real and approximation %%%%%%
figure('Renderer', 'painters', 'Position', [400 50 900 700])
hold on
plot(lambda,neff,'-g','LineWidth',2.5)

%%% Other plot settings
title('Neff(lambda)');
xlabel('Wavelenght(m)');
ylabel('Neff');
axis([lambdamin lambdamax 1.445 1.465])
ax = gca;
ax.XAxis.Exponent = -6;
set(findall(gcf,'type','text'),'FontSize',30);
set(gca,'FontSize',30);
grid;

%%%%%% Group velocity %%%%%%
figure('Renderer', 'painters', 'Position', [400 50 900 700])
hold on
plot(lambda(1:end-1),vg,'-g','LineWidth',2.5)
[VG,IN] = max(vg);
xline(lambda(IN),'-.','LineWidth',2);

%%% Insert text
xclad = lambda(IN);
yclad = 2.035e8;
str = 'Beta2 = 0';
text(xclad,yclad,str)

%%% Other plot settings
title('Vg(lambda)');
xlabel('Wavelenght(m)');
ylabel('Vg');
axis([lambdamin lambdamax 2.02*1e8 2.04*1e8])
ax = gca;
ax.XAxis.Exponent = -6;
ax.YAxis.Exponent = 8;
set(findall(gcf,'type','text'),'FontSize',30);
set(gca,'FontSize',30);
grid;