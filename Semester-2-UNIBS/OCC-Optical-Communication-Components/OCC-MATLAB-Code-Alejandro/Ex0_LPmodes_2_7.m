% LPlm modes plot
close all
clear all
clc

%% Variables declaration
V = 6;                      %%% Normalised frequency max
lmax=1;                     %%% Bessel function order l

Points=2000;
l=0:lmax;
xmax=V;
X=linspace(0,xmax,Points);  %%% Fixed variable to act as x axis X=kt*a

for i=1:length(l)
    %% RHS and LHS calculation
    Y=real(sqrt(V^2-X.^2));                        %%% Relationship X,Y
    RHS = Y.*besselk(l(i)-1,Y)./besselk(l(i),Y);    %%% Right hand side of dispersion relation 
    LHS = -X.*besselj(l(i)-1,X)./besselj(l(i),X);   %%% Left hand side of dispersion relation for l=0
    
    %% Plots
    figure('Renderer', 'painters', 'Position', [400 50 900 700])
    plot(X,LHS,'.b', 'LineWidth', 2);
    hold on
    plot(X,RHS,'-r', 'LineWidth', 2);
    axis([0 xmax+1 -1 V+1])
    title(strcat('LP', string(l(i)),'m'))
    legend('LHS', strcat('RHS, V=',string(V)))
    xlabel('X');
    ylabel('F');
    set(findall(gcf,'type','text'),'FontSize',30);
    set(gca,'FontSize',30);
    grid;
end