
%
% PROPAGATION METHOD FOR THE SOLUTION OF NLSE
%
% i Fz - beta2/2*Ftt + gamma * |F|^2 F=0
%

clc
clear all
close all

% MATERIAL PROPERTIES

beta2=0;
gamma=1;

% TEMPORAL COORDINATE

t0=-40;
t1=40;
nt=501;
t=linspace(t0,t1,nt);
deltat=(t1-t0)/(nt-1);

% SPATIAL COORDINATE
z0=0;
z1=10;
nz=101;
z=linspace(z0,z1,nz);
deltaz=(z1-z0)/(nz-1);

h=deltaz;

% FREQUENCY COORDINATE
indfreq=-nt/2:1:nt/2-1;
omega=pi/t1.*indfreq;


%INPUT ENVELOPE
amp=1;
tau0=1;
SHIFT=5;

FIN=amp*exp(-t.^2/(2*tau0^2)); %Gaussian Pulse


FFIN=deltat*fftshift(fft(FIN)); %SPECTRAL PROFILE

% figure(1)
% plot(t,abs(FIN))
% xlabel('t'); ylabel('|F|')
% set(findall(gcf,'type','text'),'FontSize',20);
% set(gca,'FontSize',20);
% grid;

% figure(2)
% plot(t,unwrap(angle(FIN)))
% xlabel('t'); ylabel('angle F')
% set(findall(gcf,'type','text'),'FontSize',30);
% set(gca,'FontSize',30);
% grid;
% 
% figure(3)
% plot(omega,abs(FFIN))
% xlabel('\omega'); ylabel('|FFIN|')
% grid
% set(findall(gcf,'type','text'),'FontSize',30);
% set(gca,'FontSize',30);
% grid;



%%%%%%%%%%%%%%
% CORE OF THE PROGRAM
%%%%%%%%%%%%%%

F=zeros(ceil(nt),floor(nz));
FF=zeros(ceil(nt),floor(nz));

F(:,1)=FIN;
FF(:,1)=FFIN;

q=FIN;

for loop_step=2:1:nz

    %LINEAR DISPERSIVE STEP
    qs=deltat*fftshift(fft(q)); %FFT
    qs_old=qs;
    
    prop=beta2/2*omega.^2;
    fact=j*prop*h;
    qs=qs_old.*exp(fact); %calculation of the envelope, undere
    %the effect of propagation beta2
    
    q=(1/deltat)*ifft(ifftshift(qs)); %coming back in the time domain
    
       
    %NONLINEAR CHI3 EFFECT  STEP
    
    q_old=q;
    q=q_old.*exp(j*gamma*abs(q_old).^2*h);
    
    %SAVE DATA EVERY NZ
    
    F(:,loop_step)=q;
    FF(:,loop_step)=deltat*fftshift(fft(q));
    
    
end

save dati


figure
mesh(t,z,abs(F)')
xlabel('t')
ylabel('z')
zlabel('|F|')
set(findall(gcf,'type','text'),'FontSize',30);
set(gca,'FontSize',30);
grid;


figure
plot(t,abs(F(:,1)),t,abs(F(:,end)),'r')
xlabel('t')
ylabel('|F|')
set(findall(gcf,'type','text'),'FontSize',30);
set(gca,'FontSize',30);
grid;


figure
mesh(t,z,unwrap(angle((F)))')
xlabel('t')
ylabel('z')
zlabel('angle F')
set(findall(gcf,'type','text'),'FontSize',30);
set(gca,'FontSize',30);
grid;

figure
mesh(omega,z,abs(FF)')
xlabel('\omega')
ylabel('z')
zlabel('|FF|')
set(findall(gcf,'type','text'),'FontSize',30);
set(gca,'FontSize',30);
grid;

figure
plot(omega,abs(FF(:,1)),omega,abs(FF(:,end)),'r')
xlabel('\omega')
ylabel('|F|')
set(findall(gcf,'type','text'),'FontSize',30);
set(gca,'FontSize',30);
grid;