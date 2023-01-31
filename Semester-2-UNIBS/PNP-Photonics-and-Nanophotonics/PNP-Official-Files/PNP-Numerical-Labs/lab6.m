% Lab 6 - A dielectric mirror with high reflectivity

% lam_bragg=3um
% d1=0.25um
% d2=0.5um
% n1=3
% n2=1.5
clear all
close all
clc

lam=linspace(1e-6,8e-6,1000);
n1=3;
n2=1.5;
d1=0.25e-6;
d2=0.5e-5;
n=[n1 n2];
d=[d1 d2];
for ii=1:length(lam)
    T_R_A=tmm(lam(ii),0,n,d,1,1,0);
    T(ii)=T_R_A(1);
    R(ii)=T_R_A(2);
    A(ii)=T_R_A(3);
end

figure()
plot(lam,R)


n=[n1 n2 n1 n2];
d=[d1 d2 d1 d2];
for ii=1:length(lam)
    T_R_A=tmm(lam(ii),0,n,d,1,1,0);
    T(ii)=T_R_A(1);
    R(ii)=T_R_A(2);
    A(ii)=T_R_A(3);
end

figure()
plot(lam*1e6,R)

% function of freq
figure()
plot(3e8./lam*1e6,R)

%% more layers
n=[n1 n2 n1 n2 n1 n2 n1 n2 n1 n2 n1 n2];
d=[d1 d2 d1 d2 d1 d2 d1 d2 d1 d2 d1 d2];
for ii=1:length(lam)
    T_R_A=tmm(lam(ii),0,n,d,1,1,0);
    T(ii)=T_R_A(1);
    R(ii)=T_R_A(2);
    A(ii)=T_R_A(3);
end

figure()
plot(lam*1e6,R)

% function of freq
figure()
plot(3e8./lam*1e6,R)