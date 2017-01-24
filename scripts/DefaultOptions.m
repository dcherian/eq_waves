function [opt, plotopt] = DefaultOptions()

    opt.filt.window = 'butterworth'; % window shape
    opt.filt.halfdef = 'power'; % how is filt.N defined?
    opt.filt.N = NaN; % will be set based on cutoff later.
    opt.filt.cutoff = 2./[0.14 0.16]; % (days) band pass filter windows
    opt.filt.debugflag = 0; % debugging spectrum plots in BandPass()

    % Save butterworth filter parameters so as to save time later.
    % butterworth requires (desired freq)/(sampling freq/2)
    % the factor of 4 is because
    % a) 1./opt.cutoff = (desired freq)/2
    %      as designed for other windows.
    % b) (sampling freq/2) = 1/2 cpd.
    [opt.filt.b, opt.filt.a] = ...
        butter(1, sort(2./opt.filt.cutoff/(1/2)), 'bandpass');

    opt.name = 'bc2m1';
    opt.nmode = 2; % what mode am I look for? used for normalization
    opt.filter_temp = 1;
    opt.debug = 0;  % debugging spectrum plots
    opt.debugRegression = 0; % debug regression
    opt.TagainstDHT = 1; % T on y-axis of regression?
    opt.InterpGapLength = 5; % gap length over which to linearly interpolate

    opt.corr_sig = 0.08; % correlation significance level
    opt.numMC = 5e3; % number of iterations when estimating error bounds

    plotopt.nmode = [2]; % which theoretical mode am I looking for?
    plotopt.plotcorr = 1;
    plotopt.plotstd = 0;
    plotopt.ploterr = 1;
    plotopt.plotOLS = 1;
    plotopt.plotWTLS = 0;
    plotopt.MarkWaterDepth = 1;
    plotopt.plotPhaseLag = 0;

    plotopt.window = opt.filt.window;
    plotopt.name = opt.name;
end