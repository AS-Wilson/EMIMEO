% takes in of time domain vector, and provides
% output of FFT in dBFS (full scale)
function out = do_fft(in)
    L = length (in); % Window Length of FFT
    in_HannWnd = in' .* hanning(L ,'periodic');
    out = 20*(log10(abs(fftshift(fft(in_HannWnd,L))/(L/2))));
end
