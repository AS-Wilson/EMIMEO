% LPlm modes plot
close all
clear all
clc

%% Variables declaration
c=3.0e8;                    %%% Light Speed
mu0=4*pi*1e-7;              %%% Magnetic permeability
eps0=8.8542e-12;            %%% Electric permitivity
a=4.5e-6;                   %%% Fiber core radius
n1=1.447;                   %%% Fiber core refractive index
n2=1.443;                   %%% Fiber cladding refractive index
Vmax = 6;                   %%% Normalised frequency max
l=0;                        %%% Bessel function order l
Points=3000;
RHS = zeros(Vmax,Points);

%%%%%%%%%%%%%%%%%%%%%%%%%%
xmax=Vmax;
X=linspace(0,xmax,Points);    %%% Fixed variable to act as x axis X=kt*a
leg = ["LHS"];

%% RHS and LHS calculation
for v=1:Vmax
    Y=real(sqrt(v.^2-X.^2));                         %%% Relationship X,Y
    LHS = -X.*besselj(l-1,X)./besselj(l,X);        %%% Left hand side of dispersion relation for l=0
    RHS(v,:) = Y.*besselk(l-1,Y)./besselk(l,Y);    %%% Right hand side of dispersion relation 
    entry = strcat('RHS, V=', string(v));
    leg = [leg, entry];
end
%% Plots
%%%%%% LP0m %%%%%%

figure('Renderer', 'painters', 'Position', [400 50 900 700])
plot(X,LHS,'.b', 'LineWidth', 2);
hold on
for v=1:Vmax
    plot(X,RHS(v,:), 'LineWidth', 2);
end
axis([0 xmax+1 -l-1 Vmax+1])
title(strcat('LP', string(l),'m'))
legend(leg)
xlabel('X');
ylabel('F');
set(findall(gcf,'type','text'),'FontSize',30);
set(gca,'FontSize',30);
grid;

%%%%%% LP1m %%%%%%
% figure
% plot(X,LHS1,X,RHS1,'.');
% axis([0 xmax+1 -4 10])
% title('LP1m')
% xlabel('X');
% ylabel('F');
% set(findall(gcf,'type','text'),'FontSize',30);
% set(gca,'FontSize',30);
% grid;
