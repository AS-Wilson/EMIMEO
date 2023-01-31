% Lab 5 - Fabry-Perot resonator

clear all
close all
clc

N_wavelenghts=100;
wavelenghts=linspace(400e-9, 700e-9,N_wavelenghts);
for ii=1:N_wavelenghts
    T_R_A=tmm(wavelenghts(ii),0,3+1i*0.2,100e-9,1.5,1,0);
    T(ii)=T_R_A(1);
    R(ii)=T_R_A(2);
    A(ii)=T_R_A(3);
end
figure(1)
plot(wavelenghts,T,wavelenghts,R)



% lab

freq=linspace(6e13,6e14,1000);
lam=3e8./freq;

%% d = 200 nm
for ii=1:length(lam)
    T_R_A=tmm(lam(ii),0,3,200e-9,1,1,0);
    T(ii)=T_R_A(1);
    R(ii)=T_R_A(2);
    A(ii)=T_R_A(3);
end

figure()
plot(freq,T,freq,R)

figure
plot(lam,T,lam,R)


%% d = 500 nm
for ii=1:length(lam)
    T_R_A=tmm(lam(ii),0,3,500e-9,1,1,0);
    T(ii)=T_R_A(1);
    R(ii)=T_R_A(2);
    A(ii)=T_R_A(3);
end

figure()
plot(lam,T,lam,R)

%% d = 800 nm
for ii=1:length(lam)
    T_R_A=tmm(lam(ii),0,3,800e-9,1,1,0);
    T(ii)=T_R_A(1);
    R(ii)=T_R_A(2);
    A(ii)=T_R_A(3);
end


figure()
plot(lam,T,lam,R)
%%




% kf*d=N*pi -> 2*pi*nf*d/lambda = N*pi

% d = 500nm -> 1st resonance @ 3 um
% d = 200nm -> 1st resonance @ 1.2 um
% d = 800nm -> 1st resonance @ 4.8 um


clear all

freq=linspace(6e13,6e14,1000);
lam=3e8./freq;
thicknesses=linspace(200e-9,800e-9,100);
for jj=1:length(thicknesses)
    for ii=1:length(lam)
    T_R_A=tmm(lam(ii),0,3,thicknesses(jj),1,1,0);
    T(ii,jj)=T_R_A(1);
    R(ii,jj)=T_R_A(2);
    A(ii,jj)=T_R_A(3);
    end
end
figure()
surf(lam*1e6,thicknesses*1e9,T');colorbar;shading interp;view(2);colormap jet

figure()
surf(lam*1e6,thicknesses*1e9,R');colorbar;shading interp;view(2);colormap jet

figure()
surf(lam*1e6,thicknesses*1e9,A');colorbar;shading interp;view(2);colormap jet
