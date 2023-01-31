% Calculos basicos

clc
clear all %close all the variables
close all %close the figures

% Constantes universales
c = 3e8; %m/s
epsilon0 = 8.854e-12; %F/m
mu0 = pi*4e-7; %H/m

n1 = 1.447;
n2 = 1.443;
a = 4.5e-6; %metros


% Vamos a ver un valor de beta en valores de V de: 2 < V < 2.3
% Cojo V = 2.15 que está en la mitad
V = 2.15;
lambda_prueba = sqrt((2*pi*a/V)*(n1^2-n2^2));
omega_prueba = 2*pi*c/lambda_prueba;
beta_1_prueba = n1*omega_prueba*sqrt(epsilon0*mu0);
beta_2_prueba = n2*omega_prueba*sqrt(epsilon0*mu0);


% Comprobacion de los resultados

% Para V = 1
V_1 = 1;
X_1 = 0.9793;
Y_1 = sqrt(V_1^2-X_1^2);
kt_1 = X_1/a;
gamma_1 = Y_1/a;

% Para V = 2
V_2 = 2;
X_2 = 1.5282;
Y_2 = sqrt(V_2^2-X_2^2);
kt_2 = X_2/a;
gamma_2 = Y_2/a;

% Para V = 6
V_6 = 6;
X_6a = 2.405;
Y_6a = sqrt(V_6^2-X_6a^2);
kt_6a = X_6a/a;
gamma_6a = Y_6a/a;

X_6b = 5.52;
Y_6b = sqrt(V_6^2-X_6b^2);
kt_6b = X_6b/a;
gamma_6b = Y_6b/a;


