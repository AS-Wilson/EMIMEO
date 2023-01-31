% Excitation of a surface plasmon polariton (SSP)
clear all
close all
clc


n_input=1.5;
n_output=1;
lambda=532;
eps_relative=sqrt(-10+i*1); %different convention for physics
d=70;

N_angles=270;
angles=linspace(0,89,N_angles);
for ii=1:N_angles
    T_R_A=tmm(lambda,angles(ii),eps_relative,d,n_input,n_output,1);
    T(ii)=T_R_A(1);
    R(ii)=T_R_A(2);
    A(ii)=T_R_A(3);
end
figure; plot(angles, R,'LineWidth',3);
xlabel('Angle of Incidence (deg)', 'FontName', 'Times New Roman', 'FontSize', 20);
ylabel('Reflectance')
set(gca, 'FontName','Time New Roman','FontSize',20)


% now we calculate things
for ii=1:1:7
    d_new=ii*10;
    for jj=1:N_angles
    T_R_A=tmm(lambda,angles(jj),eps_relative,d_new,n_input,n_output,1);
    T(jj)=T_R_A(1);
    R(jj)=T_R_A(2);
    A(jj)=T_R_A(3);
    end
    
    [R_dip(ii), index]=min(R);
    theta_dip(ii)=angles(1,index)

end

% 
% n_sp=sqrt((-10+i*1)/(-10+i*1+1))
% theta_sp=asin(1.035/1.5)*180/pi;
% k_sp=2