clc; clear all; clear time_plot; clear fft_plot; clear functions; file = 1;

% Specify plot parameters
txtsize=20;
ltxtsize=9;
pwidth=4; %4;
pheight=4; %4;
pxoffset=1; %0.65;
pyoffset=1; %0.5;
markersize=5;

% number filter taps
taps = 128;
% where to plot from
start = int32((taps/100)+1)*100;
inc = 50;

Fs1 = 10*2*pi;
% Create deterministic and stochastic digital data streams




figure(1);
plot(n(start:start+inc)*Fs1, sin_wave(start:start+inc), 'o', ...
    n(start:start+inc)*Fs1, sin_wave(start:start+inc));
time_plot(start, start+inc, txtsize, ltxtsize, pwidth, pheight, pxoffset, pyoffset, ...
    markersize, {'Original', 'Sinusoidal Signal : Time domain'});

figure(2);
plot(n(start:start+inc)*Fs1, random(start:start+inc),'o', ...
    n(start:start+inc)*Fs1, random(start:start+inc));
time_plot(start, start+inc, txtsize, ltxtsize, pwidth, pheight, pxoffset, pyoffset, ...
    markersize, {'Original', 'Random Binary Signal : Time domain'});

figure(3);
fft_plot(sin_wave(start:end), 1, txtsize, ltxtsize, pwidth, pheight, ...
    pxoffset,  pyoffset, markersize, ...
    {'Original', 'Sinusoidal Signal : Fourier domain'}, '$f_{s1}$');

figure(4);
fft_plot(random(start:end), 1, txtsize, ltxtsize, pwidth, pheight, ...
    pxoffset, pyoffset, markersize, ...
    {'Original', 'Random Binary Signal : Fourier domain'}, '$f_{s1}$');

%% Create lowpass filter and apply it to both data streams
% b = firls(n,f,a),
%     n is the FIR filter order
%     f is a vector of pairs of frequency points,
%     a is a vector containing the desired amplitude at the points in f




figure(5)
plot(n(start:start+inc)*Fs1, sin_bwlimited(start:start+inc), '-o');
time_plot(start, start+inc, txtsize, ltxtsize, pwidth, pheight, pxoffset, ...
    pyoffset, markersize, ...
    {'Band Limited', 'Sinusoidal Signal : Time domain'});
figure(6)
plot(n(start:start+inc)*Fs1, random_bwlimited(start:start+inc),'o');
time_plot(start, start+inc, txtsize, ltxtsize, pwidth, pheight, pxoffset, ...
    pyoffset, markersize, ...
    {'Band Limited', 'Random Binary Signal : Time domain'});

figure(7);
fft_plot(sin_bwlimited(start:end), 1, txtsize, ltxtsize, pwidth, pheight, ...
    pxoffset, pyoffset, markersize, ...
    {'Band Limited', 'Sinusoidal Signal : Fourier domain'}, '$f_{s1}$');
figure(8);
fft_plot(random_bwlimited(start:end), 1, txtsize, ltxtsize, pwidth, ...
    pheight, pxoffset, pyoffset, markersize, ...
    {'Band Limited', 'Random Binary Signal : Fourier domain'}, '$f_{s1}$');

%% y = upsample(x,n)
%     increases the sampling rate of x by inserting (n-1) zeros
%     between samples.





figure(9)
plot(n(start*N:(start+inc)*N)*Fs1, sin_up(start*N:(start+inc)*N), 'o');
time_plot(start*N, (start+inc)*N, txtsize, ltxtsize, pwidth, pheight, ...
    pxoffset, pyoffset, markersize, ...
    {sprintf('Improperly Upsampled (%d)', N), 'Sinusoidal Signal : Time domain'});
figure(10)
plot(n(start*N:(start+inc)*N)*Fs1, random_up(start*N:(start+inc)*N),'o');
time_plot(start*N, (start+inc)*N, txtsize, ltxtsize, pwidth, pheight, ...
    pxoffset, pyoffset, markersize, ...
    {sprintf('Improperly Upsampled (%d)', N),  'Random Binary Signal : Time domain'});

figure(11);
fft_plot(sin_up(start:end), N, txtsize, ltxtsize, pwidth, pheight, ...
    pxoffset, pyoffset, markersize, ...
    {sprintf('Improperly Upsampled (%d)', N), 'Sinusoidal Signal : Fourier domain'}, '$f_{s2}$');
figure(12);
fft_plot(random_up(start:end), N, txtsize, ltxtsize, pwidth, pheight, pxoffset, ...
    pyoffset, markersize, ...
    {sprintf('Improperly Upsampled (%d)', N), 'Random Binary Signal : Fourier domain'}, '$f_{s2}$');

%% Attempt to downsampling by M without filtering
% This is incorrect, but is instructive to show what artifacts occur





figure(13)
plot(n(start*N/M:(start+inc)*N/M)*Fs1, sin_up_down(start*N/M:(start+inc)*N/M), 'o');
time_plot(start*N/M, (start+inc)*N/M, txtsize, ltxtsize, pwidth, pheight, pxoffset, pyoffset, ...
    markersize, {sprintf('Improperly Upsampled (%d) then Downsampled (%d)', N, M), 'Sinusoidal Signal : Time domain'});
figure(14)
plot(n(start*N/M:(start+inc)*N/M)*Fs1, random_up_down(start*N/M:(start+inc)*N/M),'o');
time_plot(start*N/M, (start+inc)*N/M, txtsize, ltxtsize, pwidth, pheight, pxoffset, pyoffset, ...
    markersize, {sprintf('Improperly Upsampled (%d) then Downsampled (%d)', N, M), 'Random Binary Signal : Time domain'});

figure(15);
fft_plot(sin_up_down(start:end), N, txtsize, ltxtsize, pwidth, pheight, pxoffset, ...
    pyoffset, markersize, ...
    {sprintf('Improperly Upsampled (%d) then Downsampled (%d)', N, M),  'Sinusoidal Signal : Fourier domain'}, '$f_{s3}$');
figure(16);
fft_plot(random_up_down(start:end), N, txtsize, ltxtsize, pwidth, pheight, ...
    pxoffset, pyoffset, markersize, ...
    {sprintf('Improperly Upsampled (%d) then Downsampled (%d)', N, M), 'Random Binary Signal : Fourier domain'}, '$f_{s3}$');

%% Lowpass filtering of baseband periodic replica followed by downsampling
% (correct approach)





start = start + int32(taps/(2*N));
figure(17)
plot(n(start*N:(start+inc)*N)*Fs1, sin_up_filtered(start*N:(start+inc)*N), 'o');
time_plot(start*N, (start+inc)*N, txtsize, ltxtsize, pwidth, pheight, pxoffset, pyoffset, ...
    markersize, {sprintf('Properly Upsampled (%d) and Filtered', N), 'Sinusoidal Signal : Time domain'});
figure(18)
plot(n(start*N:(start+inc)*N)*Fs1, random_up_filtered(start*N:(start+inc)*N),'o');
time_plot(start*N, (start+inc)*N, txtsize, ltxtsize, pwidth, pheight, pxoffset, pyoffset, ...
    markersize, {sprintf('Properly Upsampled (%d) and Filtered', N), 'Random Binary Signal : Time domain'});

figure(19);
fft_plot(sin_up_filtered(start:end), 1, txtsize, ltxtsize, pwidth, pheight, pxoffset, ...
    pyoffset, markersize, ...
    {sprintf('Properly Upsampled (%d) and Filtered', N),  'Sinusoidal Signal : Fourier domain'}, '$f_{s2}$');
figure(20);
fft_plot(random_up_filtered(start:end), 1, txtsize, ltxtsize, pwidth, pheight, ...
    pxoffset, pyoffset, markersize, ...
    {sprintf('Properly Upsampled (%d) and Filtered', N), 'Random Binary Signal : Fourier domain'}, '$f_{s2}$');

figure(21)
plot(n(start*N/M:(start+inc)*N/M)*Fs1, sin_up_filtered_down(start*N/M:(start+inc)*N/M), 'o');
time_plot(start*N/M, (start+inc)*N/M, txtsize, ltxtsize, pwidth, pheight, pxoffset, pyoffset, ...
    markersize, {sprintf('Properly Upsampled (%d), Filtered, and Downsampled (%d)', N, M), 'Sinusoidal Signal : Time domain'});
figure(22)
plot(n(start*N/M:(start+inc)*N/M)*Fs1, random_up_filtered_down(start*N/M:(start+inc)*N/M),'o');
time_plot(start*N/M, (start+inc)*N/M, txtsize, ltxtsize, pwidth, pheight, pxoffset, pyoffset, ...
    markersize, {sprintf('Properly Upsampled (%d), Filtered, and Downsampled (%d)', N, M), 'Random Binary Signal : Time domain'});

figure(23);
fft_plot(sin_up_filtered_down(start:end), 1, txtsize, ltxtsize, pwidth, pheight, pxoffset, ...
    pyoffset, markersize, ...
    {sprintf('Properly Upsampled (%d), Filtered and Downsampled (%d)', N, M),  'Sinusoidal Signal : Fourier domain'}, '$f_{s3}$');
figure(24);
fft_plot(random_up_filtered_down(start:end), 1, txtsize, ltxtsize, pwidth, pheight, ...
    pxoffset, pyoffset, markersize, ...
    {sprintf('Properly Upsampled (%d), Filtered and Downsampled (%d)', N, M), 'Random Binary Signal : Fourier domain'}, '$f_{s3}$');

clear time_plot; clear fft_plot; clear functions;
Vars=whos;
PersistentVars=Vars([Vars.persistent]);
PersistentVarNames={PersistentVars.name};
clear(PersistentVarNames{:});
