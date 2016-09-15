%% first figure out modes

% main options
opt.debug = 0; % debugging spectrum plots
opt.butterworth = 1; % butterworth filter if 1, else running mean
opt.filter_temp = 1; % filter temperature also?
opt.name = 'bc2m1';
opt.windows = [5 9]; % (days) band pass filter windows
opt.n_modes = 3; % number of modes to calculate
opt.n_mode = 2; % which theoretical mode am I looking for?

InferModeShape(opt);

%% Plot bc2m1 shape

PlotModeMap('bc2m1');
