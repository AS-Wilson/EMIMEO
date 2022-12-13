function fft_plot(data, points, txtsize, ltxtsize, pwidth, pheight, pxoffset, ...
    pyoffset, markersize, titlestr, units)
    persistent file;
    fft_data = do_fft(data);
    fft_axis=do_fs(data);
    [psor,lsor] = findpeaks(fft_data(int32(length(data)/2):length(data)), ...
        fft_axis(int32(length(data)/2):length(data)), ...
    'NPeaks', points, 'SortStr', 'descend');
    if mean(psor) > -23
        [M,I] = min(lsor);
        plot(fft_axis, fft_data, lsor(I), psor(I), 'o');
        xlabel('$\frac{f_s}{2}$','interpreter','latex');
        str = sprintf('%2.1fdB @ %.3f', psor(I), lsor(I));
        text(lsor(I)+.02,psor(I), ...
            strcat(str, '$\frac{f_s}{2}$'), ...
            'interpreter','latex');
    else
        plot(fft_axis, fft_data);
    end
    xlabel(strcat({'Frequency'}, {' '}, {units}), 'interpreter','latex');
    ylabel('Signal Amplitude (dB)');
    ylim([-200 0]);
    xlim([-0.5 0.5]);
    %xticks([-.5 -.25  0  .25  .5]);
    grid on;
    title(titlestr);
end
