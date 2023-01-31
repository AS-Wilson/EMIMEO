% Gaussian approximation for LP01 mode 
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%
Points = 3000;
r = linspace(0,4*a,Points);
l0 = 0;                       %%% Bessel function order l
Vmax = 2;                     %%% Maximum normalised frequency

%%%%%%%%%%%%%%%%%%%%%%%%%%%
for V = 1:Vmax
    %% RHS and LHS calculation
    X = linspace(0,V,Points);                   %%% Fixed variable to act as x axis X=kt*a              
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
        end
    end

    % Report to command window:
    fprintf('Closest distance is %f which occurs in X = %.4d\n',...
      minDistance, X(minXindex));

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

    M_real = zeros(1,Points);
    for i = 1:length(r)
        if r(i) < a
            M_real(i) = besselj(l0,ktLP01*r(i));
        else
            M_real(i) = besselj(l0,ktLP01*a)/besselk(l0,gammaLP01*a)*...
                besselk(l0,gammaLP01*r(i));
        end
    end

    %% Gaussian approximation calculation
    wo = a*(0.65+1.619*V^(-3/2)+2.879*V^(-6));
    M_gaussian = exp(-r.^2/wo^2);
    neff_approx = n2+(n1-n2)*(1.1428-0.996/V)^2;

    %% Plots

    %%%%%% M(r) for LP01 and Gaussian approximation %%%%%%
    figure('Renderer', 'painters', 'Position', [400 50 900 700])
    hold on
    plot(r,M_real,'-g','LineWidth',2.5)
    plot(r,M_gaussian,'--r','LineWidth',2.5)
    xline(a,'-.','LineWidth',2);
    
    xcore = 0.4*a;
    ycore = 1.1;
    str = 'Core';
    text(xcore,ycore,str)
    xclad = 1.1*a;
    yclad = 1.1;
    str = 'Cladding';
    text(xclad,yclad,str)

    title(sprintf('M(r) for V = %i',V))
    xlabel('r(m)');
    ylabel('|M|');
    axis([0 3*a -0.25 1.25])
    legend('Real', 'Gaussian')
    ax = gca;
    ax.XAxis.Exponent = -6;
    set(findall(gcf,'type','text'),'FontSize',30);
    set(gca,'FontSize',30);
    grid;
end
