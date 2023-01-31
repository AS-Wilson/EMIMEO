% Sellmeier equations -> ne(lambda) and vg(lambda)
close all
clear all
clc

%% Variables declaration
c = 3.0e8;                    %%% Light speed
mu0 = 4*pi*1e-7;              %%% Magnetic permeability
eps0 = 8.8542e-12;            %%% Electric permitivity
a = 4.5e-6;                   %%% Fiber core radius
n1 = 1.447;                   %%% Fiber core refractive index
n2 = 1.443;                   %%% Fiber cladding refractive index
Vmax = 2;                     %%% Maximum normalised frequency
l0 = 0;                       %%% Bessel function order l
%%%%%%%%%%%%%%%%%%%%%%%%%%%
Points = 6000;
r = linspace(0,4*a,Points);
M_real = zeros(Vmax,Points);
leg = [];

for V = 1:Vmax
    %% RHS and LHS calculation
    X = linspace(0,V,Points);  %%% Fixed variable to act as x axis X=kt*a              
    Y = real(sqrt(V.^2-X.^2));          %%% Relationship X,Y
    LHS0 = -X.*besselj(l0-1,X)./besselj(l0,X);   %%% Left hand side of dispersion relation 
    RHS0 = Y.*besselk(l0-1,Y)./besselk(l0,Y);    %%% Right hand side of dispersion relation 

    %% Find LP01 solution
    minDistance = inf;
    distances = zeros(1,length(LHS0));
    for ni = 1 : length(LHS0)
        distances(ni) = abs(LHS0(ni) - RHS0(ni));
        if distances(ni) < minDistance
          minXindex = ni;
          minDistance = distances(ni);
        end
    end

%     % Report to command window:
%     fprintf('Closest distance is %f which occurs in X = %.4d\n',...
%       minDistance, X(minXindex));

    %% Propagation parameters calculation
    XLP01 = X(minXindex);
    YLP01 = Y(minXindex);
    wLP01 = sqrt((XLP01^2 + YLP01^2)/(mu0*eps0*(n1^2-n2^2)*a^2));
    ktLP01 = XLP01/a;
    gammaLP01 = YLP01/a;
    beta0 = (wLP01/c*n1 + wLP01/c*n2)/2;
    delta_beta = (n1^2*wLP01^2*mu0*eps0-beta0^2-ktLP01^2)/(2*beta0);
    beta = beta0+delta_beta;
    neff = beta*c/wLP01;

    %% M(r) calculation
    for i = 1:length(r)
        if r(i) < a
            M_real(V,i) = besselj(l0,ktLP01*r(i));
        else
            M_real(V,i) = besselj(l0,ktLP01*a)/besselk(l0,gammaLP01*a)*...
                besselk(l0,gammaLP01*r(i));
        end
    end
    entry = strcat('V=', string(V));
    leg = [leg, entry];
end

%% Plots
%%%%%% M(r) for LP01 %%%%%%

figure('Renderer', 'painters', 'Position', [400 50 900 700])
hold on
for v=1:Vmax
    plot(r,M_real(v,:), 'LineWidth', 2);
end
xline(a,'-.','LineWidth',2);

xcore = 0.4*a;
ycore = 1.3;
str = 'Core';
text(xcore,ycore,str)
xclad = 1.1*a;
yclad = 1.3;
str = 'Cladding';
text(xclad,yclad,str)

title('M(r)')
xlabel('r(m)');
ylabel('E');
axis([0 3*a -0.5 1.5])
legend(leg)
ax = gca;
ax.XAxis.Exponent = -6;
set(findall(gcf,'type','text'),'FontSize',30);
set(gca,'FontSize',30);
grid;

