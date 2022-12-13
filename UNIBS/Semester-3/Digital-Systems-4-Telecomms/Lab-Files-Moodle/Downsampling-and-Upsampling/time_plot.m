function time_plot(x1, x2, txtsize, ltxtsize, pwidth, pheight, pxoffset, ...
    pyoffset, markersize, titlestr)
    persistent file;
    xlabel('Discrete Time (n)');ylabel('Signal Amplitude');
    xlim([x1 x2]);
    title(titlestr);
end
