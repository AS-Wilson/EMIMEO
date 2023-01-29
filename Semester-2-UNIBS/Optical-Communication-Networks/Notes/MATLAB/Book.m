% SSFM code for solving the normalized NLS equation 
% $ By G. P. Agrawal for the 6th edition of NLFO book
fiblen = 5; % fiber length (in units of L_D)
beta2 = -1; % sign of GVD parameter beta_2
N = 1; % soliton order

%---set simulation parameters
nt = 1024; Tmax = 32; % FFT points and window size 
step_num = round(20*fiblen*N^2); % No. of z steps
deltaz = fiblen/step_num; % step size in z
dtau = (2*Tmax)/nt; % step size in tau

tau = (-nt/2:nt/2-1)*dtau; % time array
omega = fftshift(-nt/2:nt/2-1)*(pi/Tmax); % omega array 
uu = sech(tau); % sech pulse shape (can be modified)

%---Plot input pulse shape and spectrum
temp = fftshift(ifft(uu)); % Fourier transform 
spect = abs(temp).^2; % input spectrum
spect = spect./max(spect); % normalize
freq = fftshift(omega)/(2*pi); % freq. array

figure(1)
subplot(2,1,1);
    plot(tau, abs(uu).^2, '--k'); hold on; axis([-5 5 0 inf]);
    xlabel('Normalized Time'); ylabel('Normalized Power'); 
subplot(2,1,2);
    plot(freq, spect, '--k'); hold on; axis([-.5 .5 0 inf]); 
    xlabel('Normalized Frequency'); ylabel('Spectral Power');


%---store dispersive phase shifts to speedup code
dispersion = exp(0.5i*beta2*omega.^2*deltaz); % phase factor 
hhz = 1i*N^2*deltaz; % nonlinear phase factor

%*********[ Beginning of MAIN Loop]***********
% scheme: 1/2N -> D -> 1/2N; first half step nonlinear
temp = uu.*exp(abs(uu).^2.*hhz/2); % note hhz/2
for n=1:step_num
    f_temp = ifft(temp).*dispersion; 
    uu = fft(f_temp);
    temp = uu.*exp(abs(uu).^2.*hhz);
end
uu = temp.*exp(-abs(uu).^2.*hhz/2); % Final field
%***************[ End of MAIN Loop ]**************

%----Plot output pulse shape and spectrum
temp = fftshift(ifft(uu)); % Fourier transform
spect = abs(temp).^2; % output spectrum
spect = spect./max(spect); % normalize

figure(2)
subplot(2,1,1)
    plot(tau, abs(uu).^2, '-b'); hold off;
    xlabel('Normalized Time'); ylabel('Normalized Power');
subplot(2,1,2)
    plot(freq, spect, '-b'); hold off;
    xlabel('Normalized Frequency'); ylabel('Spectral Power');



