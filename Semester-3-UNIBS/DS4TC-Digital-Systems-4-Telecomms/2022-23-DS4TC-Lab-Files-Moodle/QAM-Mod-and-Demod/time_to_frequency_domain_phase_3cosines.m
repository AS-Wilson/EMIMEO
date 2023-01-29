% Plot sum-of-three-sines in the time and frequency domains (complex)

% set parameters
fs = 10240;                     % sampling rate
T_max = 1;                      % sim end time
t = (0:1/fs:(T_max-1/fs))';     % time vector
f1 = 100;                       % frequency of 1st tone 
f2 = 200;                       % frequency of 2nd tone
f3 = 300;
A1 = 10;                        % amplitudes of tones
A2 = 1;
A3 = 4;

% create sum of cosines signal
s_t = A1*cos((2*pi*f1*t)+pi/4) + A2*cos((2*pi*f2*t)+pi/6) + + A3*cos(2*pi*f3*t);      % cosines


% create very small impulse and add to signal (helps to ensure correctness 
% of phase spectrum)
r_t = zeros(length(s_t),1);   
r_t(1) = 0.01;                % try commenting out this line - see what happens!
y_t = s_t + r_t;

% perform FFT using rectangular window... 
Nfft = 1024;                         % number of points in FFT
X_1 = (1/Nfft) * fftshift(fft(y_t,Nfft));        % calculate FFT
f = ((-Nfft/2):(Nfft/2)-1)/Nfft*fs;              % create frequency scale

% plot time domain waveform
Scale = 1;
figure(101); 
plot(t,s_t,'b-','LineWidth',2);
axis([0 0.03 -20 20]);
xlabel('time (s)');
ylabel('amplitude');
title('Sum of Cosines');
grid on;

% plot figure (frequency domain, 2 sided)
figure(102); 
stem(f,Scale*real(X_1),'r-o','MarkerSize',5,'MarkerFaceColor','r','LineWidth',2);
axis([-400 400 0 10]);
xlabel('frequency (Hz)');
ylabel('real part');
title('FFT with rectangular window');
grid on;

% plot figure (frequency domain, 2 sided)
figure(103); 
stem(f,Scale*imag(X_1),'r-o','MarkerSize',5,'MarkerFaceColor','r','LineWidth',2);
axis([-400 400 -10 10]);
xlabel('frequency (Hz)');
ylabel('imag part');
title('FFT with rectangular window');
grid on;

% plot FFT (frequency spectrum, linear)
figure(104); 
stem(f,Scale*abs(X_1),'r-o','MarkerSize',5,'MarkerFaceColor','r','LineWidth',2);
axis([-400 400 0 10]);
xlabel('frequency (Hz)');
ylabel('magnitude');
title('FFT with rectangular window');
grid on;

% plot FFT (frequency spectrum, linear)
figure(105); 
stem(f,rem(angle(X_1),pi),'r-o','MarkerSize',5,'MarkerFaceColor','r','LineWidth',2);
axis([-400 400 -pi pi]);
xlabel('frequency (Hz)');
ylabel('phase (radians)');
title('FFT with rectangular window');
grid on;