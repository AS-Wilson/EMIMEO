% Ricardo Escobar
% Chapter 1

clc; clear all; close all;

li = 1e-6; lf = 2e-6; nl = 15;
lambda = linspace(li,lf,nl);
mu0 = 4*pi*1e-7;                                  
eps0 = 8.8542e-12;  
c = 3e8; % m/s

a = 2.15e-6; %Fiber core radius

%Data from Table 1.4: Sellmeier parameters

A_si = [0.6961663, 0.4079426, 0.8974994];
A_ge = [0.7136824, 0.4254807, 0.8964226];
lambda_si = [0.0684043, 0.1162414, 9.8961609]*1e-6; 
lambda_ge = [0.0617167, 0.1270814, 9.8961614]*1e-6;
A_BK7 = [1.03961212, 0.231792344, 1.01046945]; % BK7 parameters
lambda_BK7 = [sqrt(6.00069867e-3), sqrt(2.00179144e-2), sqrt(103.560653)]*1e-6; % BK7 wavelengths

n_si = zeros(1,nl);
sum_si = 0; 
sum_ge = 0; n_ge = n_si;
sum_BK7 = 0; n_BK7 = n_si;
beta_G = n_si; neff_G = n_si;

% Now we introduce the fiber calculations
nX = 50000;
l = 0;

for i=1:nl
    for j = 1:nl
        sum_si = 0; sum_ge = 0; sum_BK7 = 0;
        for k = 1:length(A_si)
            sum_si = sum_si + A_si(k)*lambda(j)^2/(lambda(j)^2-lambda_si(k)^2);
            sum_ge = sum_ge + A_ge(k)*((lambda(j)^2)/((lambda(j).^2)-(lambda_ge(k)^2)));
            sum_BK7 = sum_BK7 + A_BK7(k)*((lambda(j)^2)/((lambda(j).^2)-(lambda_BK7(k)^2)));
        end
        n_si(j)=sqrt(1+sum_si);
        n_ge(j)=sqrt(1+sum_ge);
        n_BK7(j)=sqrt(1+sum_BK7);
    end
    
    V = 2*pi*a*sqrt((n_ge(i)^2 - n_si(i)^2))/lambda(i);
    X = linspace(0, V, nX);
    Y = real(sqrt(V.^2-X.^2));

    LHS = -X.*besselj(l-1, X)./besselj(l, X);
    RHS = Y.*besselk(l-1,Y)./besselk(l,Y);
    
    minDist = inf;
    minXindex = 0;
    dist = zeros(1,length(LHS));
    for j = 1:length(LHS)
        dist(j) = abs(LHS(j) - RHS(j));
        if dist(j) < minDist
          minXindex = j;
          minDist = dist(j);
        else
          break
        end
    end
    Xsol = X(minXindex);
    Ysol = Y(minXindex);
    w_G = sqrt((Xsol^2 + Ysol^2)/(mu0*eps0*(n_ge(i)^2-n_si(i)^2)*a^2));
    kt_G = Xsol/a;
    gamma_G = Ysol/a;
    beta0 = (w_G/c*n_ge(i) + w_G/c*n_si(i))/2;
    dBeta = (n_ge(i)^2*w_G^2*mu0*eps0-beta0^2-kt_G^2)/(2*beta0);
    beta_G(i) = beta0+dBeta;
    neff_G(i) = beta_G(i)/(w_G*sqrt(mu0*eps0));
end

w_si = 2*pi*c./lambda;
beta_si = w_si.*n_si/c;
beta1_si = diff(beta_si)./diff(w_si);
beta2_si = diff(beta1_si)./diff(w_si(1:end-1));
beta3_si = diff(beta2_si)./diff(w_si(1:end-2));
Vg_si = diff(w_si)./diff(beta_si);

w_ge = 2*pi*c./lambda;
beta_ge = w_ge.*n_ge/c;
beta1_ge = diff(beta_ge)./diff(w_ge);
beta2_ge = diff(beta1_ge)./diff(w_ge(1:end-1));
beta3_ge = diff(beta2_ge)./diff(w_ge(1:end-2));
Vg_ge = diff(w_ge)./diff(beta_ge);

w_BK7 = 2*pi*c./lambda;
beta_BK7 = w_BK7.*n_BK7/c;
beta1_BK7 = diff(beta_BK7)./diff(w_BK7);
beta2_BK7 = diff(beta1_BK7)./diff(w_BK7(1:end-1));
beta3_BK7 = diff(beta2_BK7)./diff(w_BK7(1:end-2));
Vg_BK7 = diff(w_BK7)./diff(beta_BK7);

omega_G = 2*pi*c./lambda;
beta1_G = diff(beta_G)./diff(omega_G);
beta2_G = diff(beta1_G)./diff(omega_G(1:end-1));
beta3_G = diff(beta2_G)./diff(omega_G(1:end-2));
Vg_G = diff(omega_G)./diff(beta_G);

%% 1.11

figure
grid on; box on; hold on;
plot(lambda*1e6, n_si, 'LineWidth', 2)
plot(lambda*1e6, n_ge, 'LineWidth', 2)
plot(lambda*1e6, n_BK7, 'LineWidth', 2)
% xlim([1,1.6]); ylim([1.44,1.465]);
xlabel('Wavelength [\mu m]'); ylabel('n');
set(gca,'FontSize',16);
legend('100% SiO2','8% GeO2, 92% SiO2','BK7');

figure
grid on; box on; hold on;
plot(lambda(2:end)*1e6, Vg_si, 'LineWidth', 2)
plot(lambda(2:end)*1e6, Vg_ge, 'LineWidth', 2)
plot(lambda(2:end)*1e6, Vg_BK7, 'LineWidth', 2)
% xlim([1, 1.6]); ylim([2.03, 2.055]);
xlabel('Wavelength [\mu m]'); ylabel('Vg [m/s]');
set(gca,'FontSize',16);
legend('100% SiO2','8% GeO2, 92% SiO2','BK7');

%% 1.12

figure
grid on; box on; hold on;
plot(lambda(1:end-2)*1e6, beta2_si*1e27, 'LineWidth', 2)
plot(lambda(1:end-2)*1e6, beta2_ge*1e27, 'LineWidth', 2) 
plot(lambda(1:end-2)*1e6, beta2_BK7*1e27, 'LineWidth', 2)
plot(lambda(1:end-2)*1e6, beta2_G*1e27, 'LineWidth', 2)
xlim([1,1.6]);
xlabel('Wavelength [\mu m]'); ylabel('K'''' [ps^2/km]');
set(gca,'FontSize',16);
legend({'100% SiO2','8% GeO2, 92% SiO2','BK7', 'Guided'}, 'Location','southwest');

figure
grid on; box on; hold on;
plot(lambda(1:end-3)*1e6, beta3_si*1e40, 'LineWidth', 2)
plot(lambda(1:end-3)*1e6, beta3_ge*1e40, 'LineWidth', 2) 
plot(lambda(1:end-3)*1e6, beta3_BK7*1e40, 'LineWidth', 2)
% xlim([1,1.6]); %ylim([0,0.2]);
xlabel('Wavelength [\mu m]'); ylabel('K'''''' (ps^3/km)');
set(gca,'FontSize',16);
legend({'100% SiO2','8% GeO2, 92% SiO2', 'BK7'},'Location','northwest');
