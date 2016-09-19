%% first figure out modes

clear opt

% main options
opt.debug = 0; % debugging spectrum plots
opt.filter_temp = 1; % filter temperature also?
opt.name = 'bc2m1';
opt.n_modes = 3; % number of modes to calculate
opt.n_mode = 2; % which theoretical mode am I looking for?

opt.filt.cutoff = 2./[0.135 0.155]; % (days) band pass filter windows
opt.filt.window = 'rect'; % window shape
opt.filt.halfdef = 'power'; % how is filt.N defined?
opt.filt.N = NaN; % will be set based on cutoff later.
opt.filt.debugflag = 0; % debugging spectrum plots in BandPass()

InferModeShape(opt);

%% Plot bc2m1 shape

PlotModeMap('bc2m1');
