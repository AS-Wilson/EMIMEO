close all
clear all
clc

eps_inf=9.5;
omega_p=1.36e16;
gamma=1.05e14;
lam=linspace(300e-9,2000e-9,10000);
freq=3e8./lam;
omega=2*pi*freq;


%step 1
eps_gold=eps_inf-omega_p^2./(omega.^2-i.*omega.*gamma);


figure(1)
%plot(lam, eps_gold, 'LineWidth',3)

plot(lam,real(eps_gold), lam, imag(eps_gold))
%plot(lam,real(eps_gold))
%plot(lam,imag(eps_gold))
axis([300e-9 2000e-9 -200 100]);


%step 2
eps_d= 1; %air
n_sp=sqrt((eps_gold.*eps_d)./(eps_gold+eps_d));
k_sp=omega./3e8.*n_sp;

figure(2)
plot(real(k_sp),omega)