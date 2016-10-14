function [opt, plotopt] = DefaultOptions()

    opt.filt.window = 'butterworth'; % window shape
    opt.filt.halfdef = 'power'; % how is filt.N defined?
    opt.filt.N = NaN; % will be set based on cutoff later.
    opt.filt.cutoff = 2./[0.14 0.16]; % (days) band pass filter windows
    opt.filt.debugflag = 0; % debugging spectrum plots in BandPass()

    opt.name = 'bc2m1';
    opt.nmode = 2; % what mode am I look for? used for normalization
    opt.filter_temp = 1;
    opt.debug = 0;  % debugging spectrum plots

    plotopt.nmode = [2]; % which theoretical mode am I looking for?
    plotopt.plotcorr = 1;
    plotopt.plotstd = 0;
    plotopt.ploterr = 1;
    plotopt.plotOLS = 1;
    plotopt.plotWTLS = 1;

    plotopt.window = opt.filt.window;
    plotopt.name = opt.name;
end