%% first figure out modes

clear opt

% main options
opt.debug = 0; % debugging spectrum plots
opt.filter_temp = 1; % filter temperature also?
opt.n_modes = 3; % number of modes to calculate

opt.filt.window = 'gauss'; % window shape
opt.filt.halfdef = 'power'; % how is filt.N defined?
opt.filt.N = NaN; % will be set based on cutoff later.
opt.filt.debugflag = 0; % debugging spectrum plots in BandPass()

%% bc2m1

opt.name = 'bc2m1';
opt.n_mode = 2; % which theoretical mode am I looking for?
opt.filt.cutoff = 2./[0.135 0.155]; % (days) band pass filter windows
InferModeShape(opt);
PlotModeMap(opt.name);

export_fig -r300 images/09-21-bc2m1.png

%% bc2m2

opt.name = 'bc2m2';
opt.n_mode = 2;  % which theoretical mode am I looking for?
opt.filt.cutoff = 2./[0.19 0.21]; % (days) band pass filter windows
InferModeShape(opt);
PlotModeMap(opt.name);

export_fig -r300 images/09-21-bc2m2.png

%% bc1m1

opt.name = 'bc1m1';
opt.n_mode = 1;  % which theoretical mode am I looking for?
opt.filt.cutoff = 2./[0.19 0.21]; % (days) band pass filter windows
InferModeShape(opt);
PlotModeMap(opt.name);

export_fig -r300 images/09-21-bc1m1.png
