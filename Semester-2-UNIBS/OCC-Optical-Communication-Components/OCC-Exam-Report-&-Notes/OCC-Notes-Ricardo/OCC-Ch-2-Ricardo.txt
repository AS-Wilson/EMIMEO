% Ricardo Escobar
% Chapter 2
 
% Normalized frequencies V = [1, 2]

clear; clc; close all;

% set(0,'defaultfigurecolor',[1 1 1])

a = 4.5e-6; % Core radius
n1 = 1.447; % Core refractive index
n2 = 1.443; % Cladding refractive index
c = 3e8;

% For M(r) we have
% M_l(r) = A1J_l(k_t*r) in the core
% M_l(r) = A2K_l(gamma*r) in the cladding

% For LPlm mode we will solve the RHS and LHS of the equation:
% -[X J_{l-1}(X)/J_{l}(X)] = [Y K_{l-1}(Y)/K_{l}(Y)]
% With V^2 = X^2 + Y^2

Xi = 0; Xf = 7; nX = 10000;
ri = 0; rf = 4*a; nr = 2*500;
li = 50e-9; lf = 5e-6; nl = 50;
lambda = linspace(li,lf,nl);
l=0; m=1; % LPlm mode
V = 2;

%% Calculations

% Dispersion relation definitions from equations
X = linspace(Xi, Xf, nX);
r = linspace(0, rf, nr); 
rIndex = find(abs(r-4.5e-6) < 0.01e-6); rCore = r(1:rIndex); rClad = r(rIndex+1:1000);
Y = real(sqrt(V^2-X.^2));
w = (V*c/a)*sqrt(1/(n1^2-n2^2));

% Bessel functions: J = besselj(NU,Z) is the Bessel function of the first kind, J_nu(Z).

RHS = -X.*besselj(l-1, X)./besselj(l, X);
LHS = Y.*besselk(l-1,Y)./besselk(l,Y);

% Xsol = 0.9793 % For V = 1
% Xsol = 1.5282 % For V = 2 
% This may fail
index = find(abs(RHS-LHS)<0.1); 
Xsol = X(index(floor(0.5*length(index))));
Ysol = real(sqrt(V^2-Xsol^2));


% By solving it graphically
beta0 = w*(n1+n2)/(2*c); % (w/c)*n2 < beta0 < (w/c)*n1
kt = Xsol/a;
gamma = Ysol/a;
beta1 = (n1*n1*w*w-((kt*c)^2))/(2*beta0*c*c) + beta0/2 ;
ne = beta1*c/w;
fprintf('The initial value beta_0 is %4.2e and the final value of beta is %4.2e resulting in an effective refractive index of %4.2e\n', beta0, beta1, ne);

% Now we can obtain M(r)
A1 = 1;
A2 = A1*besselj(l,kt*a)./besselk(l,gamma*a);
Mcore = A1*besselj(l,kt*rCore);
Mclad = A2*besselk(l,gamma*rClad);
M_r = [Mcore Mclad];

%Gaussian approximation
W0 = a*(0.65+1.619*V^(-3/2)+2.879*V^(-6));
M_G = exp(-((r./W0).^2));

% Because I did not like to plot points
q = isnan(LHS);
var1 = 1;
for i = 1:nX
    if(q(i) == 1)
        LHS(i) = 0;
    end
    if var1 == 0
        RHS(i) = NaN;
    end
    if(RHS(i)>10 && var1 == 1)
        var1 = 0;
    end
end
% End of useless stuff

% 2.10 n_e calculation

RHS_G = RHS;
ne = zeros(1,length(lambda));
ne_G = ne; beta_G = ne; w_G = ne;
j = 1;

for index = lambda
    V_G = (2*pi*a/index)*sqrt(n1^2-n2^2);
    Y_G = real(sqrt(V_G^2-X.^2));
    LHS_G = Y_G.*besselk(l-1,Y_G)./besselk(l,Y_G);
    i = 2;
    delta = RHS_G - LHS_G;
    min = (delta).*(circshift(delta,1));
    while  min(i)>0
        i = i+1 ;   
    end
    Z_G = X(i); 
    w_G(j) = (V_G*c/a)*sqrt(1/(n1^2-n2^2));
%     beta0_G = (w_G/(2*c))*(n1+n2); % (w/c)*n2 < beta0 < (w/c)*n1
    beta0_G = (w_G(j)/c)*n1 + ((w_G(j)/c)*n2-(w_G(j)/c)*n1)*rand(1,1);
    beta_G(j) = (n1*n1*w_G(j)*w_G(j)-((Z_G*c/a)^2))/(2*c*c*beta0_G)+beta0_G/2;   
    ne_G(j) = n2+(n1-n2)*((1.1428-0.996/V_G)^2); % Gaussian approximation
    j=j+1;
end
ne = beta_G*c./w_G;
Vg=diff(w_G)./diff(beta_G);

%% 2.8

figure
grid on; box on; hold on;
plot(X, RHS, 'LineWidth', 2)
plot(X, LHS, 'LineWidth', 2)
xlabel('X'); title('l=0');
axis([0 2.5 -1 3])
set(gca,'FontSize',16);
legend('LHS','RHS');

% figure
% grid on; box on; hold on;
% plot(r*1e6, M_r, 'LineWidth', 2)
% xlabel('r [\mu m]'); ylabel('M(r) [m^{-1}]'); title('l=0');
% axis([0 4*a*1e6 0 1])
% set(gca,'FontSize',16);
% legend('M(r)');

%% 2.10

% Change V = 2

figure
grid on; box on; hold on;
plot(r*1e6, M_r, 'LineWidth', 2)
plot(r*1e6, M_G, '--', 'LineWidth', 2)
xlabel('r [\mu m]'); ylabel('M(r) [m^{-1}]'); title('l=0');
axis([0 4*a*1e6 0 1])
set(gca,'FontSize',16);
legend('M(r) analytical', 'Gaussian approximation');

figure
grid on; box on; hold on;
plot(lambda*1e6,ne, 'LineWidth', 2)
plot(lambda*1e6,ne_G,'--', 'LineWidth', 2);
xlabel('wavelength(\mu m)');
ylabel('n_e');
xlim([0,4.5]);
ylim([1.443,1.447]);
set(gca,'FontSize',16);
legend('Analytical','Gaussian approximation');

%% 2.9

figure
grid on; box on; hold on;
plot(lambda*1e6,ne, 'LineWidth', 2)
plot(lambda*1e6,ne_G,'--', 'LineWidth', 2);
xlabel('wavelength [\mu m]');
ylabel('n_e');
xlim([0,4.5]);
ylim([1.443,1.447]);
set(gca,'FontSize',16);
legend('Analytical','Gaussian approximation');

figure
grid on; box on; hold on;
plot(lambda(2:end)*1e6,Vg*1e-8, 'LineWidth', 2); 
xlabel('Wavelength [\mu m]');
ylabel('Vg [10^8 m/s]');
xlim([1,1.6]);
% ylim([2.0708,2.0716]);
set(gca,'FontSize',16);
