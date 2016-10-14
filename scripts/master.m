%% first figure out modes

[opt, plotopt] = DefaultOptions;
plotopt.plotWTLS = 0;
plotopt.ploterr = 0;

InferModeShape(opt);
PlotModeMap(plotopt);
export_fig -r300 -nofontswap -depsc images/10-14-bc2m1-butter.png

%% compare windows
plotopt.window = 'rect';
PlotModeMap(plotopt);
export_fig -r300 images/09-28-bc2m1-rect.png

plotopt.window = 'gauss';
PlotModeMap(plotopt);
export_fig -r300 images/09-28-bc2m1-gauss.png

ylim([-500 0]);
export_fig -r300 images/09-28-bc2m1-top500.png

%% bc2m1 with different filtering
clear opt

% main options
opt.debug = 0; % debugging spectrum plots
opt.filter_temp = 1; % filter temperature also?

opt.filt.halfdef = 'power'; % how is filt.N defined?
opt.filt.N = NaN; % will be set based on cutoff later.
opt.filt.debugflag = 0; % debugging spectrum plots in BandPass()
opt.name = 'bc2m1';
opt.filt.cutoff = 2./[0.135 0.155]; % (days) band pass filter windows

opt.filt.window = 'rect';
InferModeShape(opt);

opt.filt.window = 'gauss';
InferModeShape(opt);

opt.filt.window = 'butterworth';
InferModeShape(opt);

%% bc2m2

opt.name = 'bc2m2';
opt.filt.cutoff = 2./[0.19 0.23]; % (days) band pass filter windows
InferModeShape(opt);

plotopt.name = opt.name;
plotopt.nmode = [1 2]; % which theoretical mode am I looking for?
plotopt.plotcorr = 0;
plotopt.plotstd = 0;
PlotModeMap(plotopt);

export_fig -r300 images/09-24-bc2m2.png
