% this script generates fake temp & dynamic height time
% series and use some of the code in InferModeShape to see how
% well I do at getting things back out.
clear;

plotfigures = 0;

[opt, plotopt] = DefaultOptions;

opt.rednoise = 1;
opt.filter_temp = 1;
plotopt.plotWTLS = 0;

% synthetic time series parameters
ntimes = 8000; % number of time steps
zfull = -5000:25:0; % full depth z-vector
zsamp = -[500 450 400 300 250 200 150 100 75 50 25]; % z vector for measurements
Tmode(:,1) = -sin(pi*zfull/max(abs(zfull/16))  + 0* 2*pi * rand(size(zfull))); % ideal mode-1 shape
Tmode(:,2) = -sin(2*pi*zfull/max(abs(zfull/16)) + 0* 2*pi * rand(size(zfull))); % ideal mode-2 shape

% generate a wave amplitude time series
% one for each wave vertical mode
% (of course these have to be at different frequencies)
rng('shuffle');
tseries = zeros([ntimes size(Tmode, 2)]);
freq = [0.08 0.1 0.145 0.2];

opt.filt.cutoff = 2./[0.14 0.16];

rednoiseamp = 2e2;
if opt.rednoise
    for ff=freq([2 4]) % mode 1
        tseries(:,1) = tseries(:,1) + ...
            rednoiseamp * rednoise([ntimes 1]) + ...
            1 * rednoise([ntimes 1]).*sin(2*pi*ff*[0:ntimes-1]');
    end
    for ff=freq([1 3]) % mode 2
        tseries(:,2) = tseries(:,2) + ...
            rednoiseamp * rednoise([ntimes 1]) + ...
            1 * rednoise([ntimes 1]).*sin(2*pi*ff*[0:ntimes-1]');
    end
else
    for ff=freq([2 4]) % mode 1
        tseries(:,1) = tseries(:,1) + ...
            rand([ntimes 1]) + ...
            2 * rand([ntimes 1]).*sin(2*pi*ff*[0:ntimes-1]');
    end
    for ff=freq([1 3]) % mode 2
        tseries(:,2) = tseries(:,2) + ...
            rand([ntimes 1]) + ...
            2 * rand([ntimes 1]).*sin(2*pi*ff*[0:ntimes-1]');
    end
end
tseries = 4 * tseries./max(tseries(:));

% add gaps
gapstart = randi([40, 100], 5);
gaplength = randi([200, 1500], 5);

for ii=1:length(gapstart)
    tseries(gapstart(ii):gapstart(ii)+gaplength(ii)) = NaN;
end

% Use amplitude and mode shape to generate T,S time series
Tsamp = nan([ntimes length(zsamp)]);
for zz=1:length(zsamp)
    Tsamp(:,zz) = 20 + ...
        sum(bsxfun(@times, tseries, Tmode(zfull == zsamp(zz),:)), 2);
    Ssamp(:,zz) = 35 + ...
        sum(bsxfun(@times, tseries, Tmode(zfull == zsamp(zz),:)), 2);
end

% calculate dynamic height time series.
dynht = trapz(zsamp, sw_svan(Ssamp, Tsamp, sw_pres(zsamp, 0)), 2);

if plotfigures
    figure('Position', [360 156 800 600]);
    hax(1) = subplot(121);
    subplot(hax(1));
    PlotSpectrum(dynht);
    PlotSpectrum(Tsamp(:,5));
end

%% BandPass filtering
dynht = BandPass(dynht, opt.filt);
if opt.filter_temp
    for zz=1:length(zsamp)
        Tsamp(:,zz) = BandPass(Tsamp(:,zz), opt.filt);
    end
end

if plotfigures
    subplot(hax(1));
    PlotSpectrum(dynht);
    PlotSpectrum(Tsamp(:,5));
    ylim([1e-30 1e3])
end
%% total least squares and plot
[infer_mode, infer_mode_error, corrcoeff, dof] = ...
    DoRegression(dynht', Tsamp', opt);

% get into structure so that PlotMode can be reused
[mmax, imax] = max(abs(infer_mode(:,1)));
imnorm = sign(infer_mode(end-1,1)) * mmax;
modes.InferredModeOLS{1,1} = infer_mode(:,1)./imnorm;
modes.InferredModeErrorOLS{1,1} = infer_mode_error(:,1)./imnorm;

[mmax, imax] = max(abs(infer_mode(:,2)));
imnorm = sign(infer_mode(end-1,2)) * mmax;
modes.InferredModeWTLS{1,1} = infer_mode(:,2)./imnorm;
modes.InferredModeErrorWTLS{1,1} = infer_mode_error(:,2)./imnorm;

modes.corr{1,1} = corrcoeff./sign(imnorm);
modes.depth{1,1} = -zsamp;
flatbot.IdealTempMode(1,1,:,:) = Tmode; %./Tmode(zfull == zsamp(1));
flatbot.zTmode = -zfull;

% Plot mode structure
if plotfigures
    hax(2) = subplot(122);
    handles = PlotMode({modes; flatbot}, 1, 1, plotopt, hax(2));
    title('Test with synthetic time series');
    handles.hleg = legend('Location', 'SouthEast');
    handles.hleg.Box = 'off';
    xlim([-1.2 1.2]);
    ylim([-750 0]);
end
