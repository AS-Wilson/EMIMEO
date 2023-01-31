% Plasmonic biosensor

clear all
close all
clc

n_input=1.5;
n_output=1;
lambda=532e-9;
n_met=sqrt(-10+i*1);
d=40e-9;
n_f=1.33;

% angle incidence: theta_inc
d_f = 30e-9;% thickness of bio-nanofilm: d_f


%step 1

n_layers=[n_met, n_f];
d_layers=[d, d_f]; % 0 is the d_f
N_angles=270;
angles=linspace(0,89,N_angles);
for ii=1:N_angles
    T_R_A=tmm(lambda,angles(ii),n_layers,d_layers,n_input,n_output,1);
    T(ii)=T_R_A(1);
    R(ii)=T_R_A(2);
    A(ii)=T_R_A(3);
end

plot(angles,R,'LineWidth',3)
 

%step 3: Linearity

d_efe=[0 10 20 30 40 50 60];
theta_dip = [44.87 46.26 47.97 49.92 52.2 54.6 59.98];
figure
plot(d_efe, theta_dip, 'LineWidth',3)





