% PROPAGATION METHOD FOR THE SOLUTION OF NLSE
%
% i Fz - beta2/2*Ftt + gamma * |F|^2 F=0
%

clc
clear
close all

% MATERIAL PROPERTIES

beta2=-1;
gamma=1;

% TEMPORAL COORDINATE

t0=-40;
t1=40;
nt=2001;
t=linspace(t0,t1,nt);
deltat=(t1-t0)/(nt-1);

% SPATIAL COORDINATE
z0=0;
z1=10;
nz=501;
z=linspace(z0,z1,nz);
deltaz=(z1-z0)/(nz-1);

h=deltaz;

% FREQUENCY COORDINATE
indfreq=-nt/2:1:nt/2-1;
omega=pi/t1.*indfreq;


%INPUT ENVELOPE

FIN=sech(t);


FFIN=deltat*fftshift(fft(FIN)); %SPECTRAL PROFILE


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
    fact=1i*prop*h;
    qs=qs_old.*exp(fact); %calculation of the envelope, undere
    %the effect of propagation beta2
    
    q=(1/deltat)*ifft(ifftshift(qs)); %coming back in the time domain
    
       
    %NONLINEAR CHI3 EFFECT  STEP
    
    q_old=q;
    q=q_old.*exp(1i*gamma*abs(q_old).^2*h);
    
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
xlim([-4 4])
ylabel('|F|')
set(findall(gcf,'type','text'),'FontSize',30);
set(gca,'FontSize',30);
grid;