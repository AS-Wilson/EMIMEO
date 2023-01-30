% SSFM code for solving the normalized NLS equation $ By G. P. Agrawal for the 6th edition of NLFO book
fiblen = 5; beta2 = -1; N = 1;
% fiber length (in units of L_D)
% sign of
% soliton
GVD parameter beta_2
order
%---set simulation parameters
nt = 1024; Tmax = 32; % FFT points and window size step_num = round(20*fiblen*N^2); % No. of z steps
deltaz = fiblen/step_num; dtau = (2*Tmax)/nt;
tau = (-nt/2:nt/2-1)*dtau;
omega = fftshift(-nt/2:nt/2-1)*(pi/Tmax); % omega array uu = sech(tau); % sech pulse shape (can be modified)
%---Plot input pulse shape and spectrum
temp = fftshift(ifft(uu)); spect = abs(temp).^2;
spect = spect./max(spect); freq = fftshift(omega)/(2*pi); subplot(2,1,1);
% Fourier transform % input spectrum
% normalize
% freq. array
% step size in z
% step size in tau
plot(tau, abs(uu).^2, ’--k’); hold on; axis([-5 5 0 inf]);
% time array
693694 Numerical code for the NLS equation
xlabel(’Normalized Time’); ylabel(’Normalized Power’); subplot(2,1,2);
plot(freq, spect, ’--k’); hold on; axis([-.5 .5 0 inf]); xlabel(’Normalized Frequency’); ylabel(’Spectral Power’);
%---store dispersive phase shifts to speedup code
dispersion = exp(0.5i*beta2*omega.^2*deltaz); % phase factor hhz = 1i*N^2*deltaz; % nonlinear phase factor
%*********[ Beginning of MAIN Loop]***********
% scheme: 1/2N -> D -> 1/2N; first half step nonlinear
temp = uu.*exp(abs(uu).^2.*hhz/2); for n=1:step_num
f_temp = ifft(temp).*dispersion; uu = fft(f_temp);
temp = uu.*exp(abs(uu).^2.*hhz);
end
uu = temp.*exp(-abs(uu).^2.*hhz/2); %***************[ End of MAIN Loop ]**************
%----Plot output pulse shape and spectrum
temp = fftshift(ifft(uu)); spect = abs(temp).^2; spect = spect./max(spect); subplot(2,1,1)
% Fourier transform
% output spectrum
% normalize
plot(tau, abs(uu).^2, ’-b’); hold off; subplot(2,1,2)
plot(freq, spect, ’-b’); hold off;