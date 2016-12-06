% driver script for figuring out null hypothesis rejection limits
% using a Monte Carlo approach.

clear

numMC = 1e3; % number of monte carlo iterations
tlen = 2e4; % length of synthetic temp time series
SpectralSlope = -2; % spectral slope of synthetic temp time series

opt = DefaultOptions;
mc = ReadTaoTriton;
mc.timetemp = mc.timedht;
opt.name = 'montecarlo';

% generate numMC time series.
% dimension 1 of syntemp == depth == 1 Monte Carlo iteration.
syntemp = synthetic_timeseries_known_spectrum([tlen numMC], 1, 1, SpectralSlope)';

for mm=1:length(mc.lon)
    for nn=1:length(mc.lat)
        disp([mm nn])
        mc.T{mm,nn} = syntemp;
        mc.depth{mm,nn} = [1:numMC]';

        modes = InferModeShape(opt, mc, mm, nn);

        if isfield(modes, 'corr')
            r(mm,nn,:) = modes.corr{mm,nn};
            m(mm,nn,:) = modes.InferredModeOLS{mm,nn};
            se(mm,nn,:) = modes.StdErrorOLS{mm,nn};
        else
            r(mm,nn,:) = NaN;
            m(mm,nn,:) = NaN;
            se(mm,nn,:) = NaN;
        end

        % save memory
        mc.temp{mm,nn} = [];
    end
end

rvec = cut_nan(r(:));
mvec = cut_nan(m(:));
sevec = cut_nan(se(:));

% fisher transformation of corr. coeff.
w = 1/2*log((1+rvec)./(1-rvec));
rfit = fitdist(w, 'normal')

% transform reg slope m. this should be t-distributed for white noise
t = mvec./sevec;
slfit = fitdist(t, 'normal')