% Exercise 1 - Figure 2.8

clc
clear all %close all the variables
close all %close the figures

% Constants
c = 3e8; % [m/s]
eps0 = 8.854e-12; % [F/m]
mu0 = pi*4e-7; % [H/m]

% Fiber optic characteristics
n1 = 1.447;
n2 = 1.443;
a = 4.5e-6; % [m]

% Variables
num_points = 10000000;
X_max = 10; % maximum of the x axis
epsilon = 0.00001; %accuracy

x = linspace(0,X_max,num_points);
J_left = -x.*besselj(-1,x)./besselj(0,x);

v = [1 2 6];

y = zeros(length(v),length(x));
J_right = zeros(length(v),length(x));

kt = zeros(1,length(v));
gamma = zeros(1,length(v));



% Calculation
leg = ["LHS"];

for i = 1:1:length(v)

    entry = strcat('RHS, V=', string(v(i))); % decoratives for the plot
    leg = [leg, entry]; % decoratives for the plot


    sol_vector = [v(i)];
    y(i,:) = real(sqrt(v(i)^2-x.^2));
    J_right(i,:) = y(i,:).*besselk(-1,y(i,:))./besselk(0,y(i,:));


    for j = 1:1:length(x)
        if J_left(j) > 100
            J_left(j) = nan;
        elseif (J_left(j) < 0 || J_left(j) > J_right(i,1))
            continue
        elseif abs((J_left(j) - J_right(i,j))) < epsilon
            if abs(sol_vector(end) - J_left(j)) < 0.0001
                continue
            else
                % in sol_vector we introduce:
                % [Value of V, Value of X, Value of Y]
                sol_vector(end+1) = x(j);
                sol_vector(end+1) = J_left(j);

                % calculate kt and gamma
                kt(i) = x(j)/a;
                %gamma(i) = J_left(j)/a;
                gamma(i) = real(sqrt(v(i)^2-x(j).^2))/a;
            end
        end
        solutions{i} = sol_vector;
    end
end


figure(1)
plot(x,J_left,'-',x,J_right)
xlim([0 8]) %xlim([0 3])
ylim([-1 8]) %ylim([-1 3])
xlabel('X','interpreter','latex');
title('$LP_{0m}$','interpreter','latex');
legend(leg);

% kt values (limited to one cut of the graphics otherwise it does not work
% x_axis = linspace(0,1,10);
% figure(2)
% for i=length(v)
%     y_axis = kt(i)*ones(size(x_axis));
%     plot(x_axis, y_axis, 'r')
%     hold on
% end


% Kt only for V = 1,2
figure(2)
yline(kt(1),'-','V = 1 and kt='+string(kt(1)))

figure(3)
yline(kt(2),'-','V = 2 and kt='+string(kt(2)))


% Gamma only for V = 1,2
figure(4)
yline(gamma(1),'-','V = 1 and gamma='+string(gamma(1)))

figure(5)
yline(gamma(2),'-','V = 2 and gamma='+string(gamma(2)))




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
% 
% We calculate beta0. It has to be between the calculation with n1 and n2
% So: (w/c)*n2 < beta0 < (w/c)*n1. We use n = (n1+n2)/2

% kt and gamma used
% kt_used = kt(1);
% gamma_used = gamma(1);
% For each V we get frequency, beta0 and beta_l
%V = 1;
w = zeros(1,length(v));
beta0 = zeros(1,length(v));
beta_l = zeros(1,length(v));

for i=1:length(v)
    w(i) = (v(i)*c/a)*sqrt(1/(n1^2-n2^2));
    beta0(i) = w(i)*(n1+n2)/(2*c);
    beta_l(i) = (n1*n1*w(i)*w(i)*eps0*mu0-kt(i)*kt(i))/(2*beta0(i)) + beta0(i)/2;
end

% % calculation of beta0
% beta0 = w*(n1+n2)/(2*c);
% % calculation of beta_l
% beta_l = (n1*n1*w*w*eps0*mu0-kt_used*kt_used)/(2*beta0) + beta0/2;
% figure(6)
% title('$\beta_{l}$','interpreter','latex')
% yline(beta_l,'-','value='+string(beta_l))


% Now we want to obtain M(r)
% Now we can obtain M(r)
n_points = 6000;
r = linspace(0,5*a,n_points);
M_real = zeros(length(v),n_points);
leg = [];

A1 = 1;
for i=1:length(v)
    for j=1:length(r)
        A2 = A1*besselj(0,kt(i)*a)./besselk(0,gamma(i)*a);
        if r(j) < a
            M_real(i,j) = A1*besselj(0,kt(i)*r(j));
        else
            M_real(i,j) = A2*besselk(0,gamma(i)*r(j));
        end
    end
    entry = strcat('V=', string(v(i)));
    leg = [leg, entry];
end


figure(6) % adding a lot of decoration
hold on
% 
% for i = length(v)
%     plot(r*1e6,M_real(i,:), 'LineWidth', 2);
% end

plot(r*1e6,M_real, '-','LineWidth',2)
legend(leg,'AutoUpdate','off');

xlim([0 20]) %xlim([0 3])
ylim([0 1.2]) %ylim([-1 3])
xlabel('r($\mu m$)','interpreter','latex');
ylabel('M(r)[$m^{-1}$]','interpreter','latex')
title('$LP_{0m}$','interpreter','latex');


xline(a*1e6,'-.','LineWidth',2) % line separating core and cladding

% next code is to put text for CORE and CLADDING
xcore = 0.4*a*1e6;
ycore = 1.1;
str = 'Core';
text(xcore,ycore,str)
xclad = 1.1*a*1e6;
yclad = 1.1;
str = 'Cladding';
text(xclad,yclad,str)








