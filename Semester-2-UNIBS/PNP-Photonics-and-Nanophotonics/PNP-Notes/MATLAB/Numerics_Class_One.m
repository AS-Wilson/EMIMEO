%% Part One %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Lambda = 600e-6;

TE_Polarisation = 0;
TM_Polarisation = 1;

n_water = 1.33;
n_output_1 = 1;

samples = 1000;


theta = linspace(0, 89, samples);


%% Scenario One
for angle = 1:samples
    TRA_Array(angle,:) = tmm(Lambda,theta(angle),1,0,n_water,n_output_1,TE_Polarisation);
end


figure(1)

title('Scenario One')

hold on

a(1) = plot(theta,TRA_Array(:,1));
a(2) = plot(theta,TRA_Array(:,2));
a(3) = plot(theta,TRA_Array(:,3));

hold off

legend(a, 'Transmission','Reflection', 'Absorption')
grid on
xlim([0 89])
ylim([-0.01 1.01])
xlabel('Angle');
ylabel('T, R, A');




%% Scenario Two
for angle = 1:samples
    TRA_Array(angle,:) = tmm(Lambda,theta(angle),1,0,n_water,n_output_1,TM_Polarisation);
end


figure(2)

title('Scenario Two')

hold on

b(1) = plot(theta,TRA_Array(:,1));
b(2) = plot(theta,TRA_Array(:,2));
b(3) = plot(theta,TRA_Array(:,3));

hold off

legend(b, 'Transmission','Reflection', 'Absorption')
grid on
xlim([0 89])
ylim([-0.01 1.01])
xlabel('Angle');
ylabel('T, R, A');




%% Scenario Three
n_air = 1;

for angle = 1:samples
    TRA_Array(angle,:) = tmm(Lambda,theta(angle),1,0,n_air,n_water,TM_Polarisation);
end


figure(3)

title('Scenario Two')

hold on

c(1) = plot(theta,TRA_Array(:,1));
c(2) = plot(theta,TRA_Array(:,2));
c(3) = plot(theta,TRA_Array(:,3));

hold off

legend(c, 'Transmission','Reflection', 'Absorption')
grid on
xlim([0 89])
ylim([-0.01 1.01])
xlabel('Angle');
ylabel('T, R, A');




%% Part Two - Step One %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
celerity = 3e8;

samples = 10000;
Lambda = linspace(300e-9, 2000e-9, samples);

freq = celerity./Lambda;

Epsilon_Infinity = 9.5;
Gamma = 1.05e14;
omega_p = 1.36e16;
omega = 2*pi*freq;


Epsilon_Gold = Epsilon_Infinity - (omega_p.^2 ./ (omega.^2 -1i.*omega.* Gamma));


figure(4)

title('Lab Exercise, Step One')

hold on
d(1) = plot(Lambda, real(Epsilon_Gold));
d(2) = plot(Lambda, imag(Epsilon_Gold));
hold off

legend(d, 'Real\{Epsilon\}','Complex\{Epsilon\}')
grid on
xlim([300e-9 2000e-9])
ylim([-200 100])
xlabel('Wavelength');
ylabel('Permittivity of Gold');




%% Part Two - Step Two %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


Surface_Plasmon = omega./celerity .* sqrt((Epsilon_Gold.*1) ./ (Epsilon_Gold + 1));

figure(5)

title('Lab Exercise, Step Two')

hold on
e(1) = plot(real(Surface_Plasmon), omega);
hold off

legend(e, 'Real\{Wavenumber Surface Plasmon\}')
grid on
%xlim([300e-9 2000e-9])
%ylim([-200 100])
xlabel('Wavenumber');
ylabel('Angular Frequency');

