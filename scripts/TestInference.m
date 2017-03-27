% this script generates fake temp & dynamic height time
% series and use some of the code in InferModeShape to see how
% well algorithm infers the signal.
% if opt.nullhyp = 1, dynht is repaced by rednoise and you should get no significant slope

function TestInference()

    plotfigures = 1;

    [opt, plotopt] = DefaultOptions;

    opt.nullhyp = 1; % if 1, dynht is substituted with red noise
    opt.filter_temp = 1;
    opt.numMC = 1000;
    plotopt.plotWTLS = 0;

    freq = [0.08 0.1 0.145 0.2]; % frequency of spectral peaks
    opt.filt.cutoff = 2./[0.14 0.16]; % filter cutoffs

    % synthetic time series parameters
    ntimes = 8000; % number of time steps
    zfull = -5000:25:0; % full depth z-vector
    zsamp = -[500 450 400 300 250 200 150 100 75 50 25]; % z vector for measurements
    nz = length(zsamp);

    [dynht, Tsamp, Tmode] = CreateData(freq, ntimes, zfull, zsamp);

    if plotfigures
        figure('Position', [360 156 800 600]);
        hax(1) = subplot(121);
        subplot(hax(1));
        PlotSpectrum(dynht);
        PlotSpectrum(Tsamp(:,5));
        linex(freq);
    end

    if opt.nullhyp
        [NoiseAmpNull, NoiseSlopeNull] =  EstimateNoiseSpectrum(dynht, opt);

        dynht = synthetic_timeseries_known_spectrum(length(Tsamp), 1, ...
                                                    NoiseAmpNull, NoiseSlopeNull);
        if plotfigures
            PlotSpectrum(dynht);
            title('Dyn. height is red noise!');
        end

        % Monte-carlo regression slope etc. with peaked time-series (temp.)
        % and random time series (dyn ht.)
        % (slope, intercept, corrcoeff) all span 0 as expected.
        % TestMC(NoiseSlopeNull, 1, 1, ntimes, Tsamp(:,6)', NoiseAmpNull);
    end

    %% This must happen before BandPass filtering.

    % "background" spectrum properties.
    if opt.TagainstDHT
        for ii=1:size(Tsamp, 2)
            [noise.amp(ii), noise.slope(ii)] = ...
                EstimateNoiseSpectrum(Tsamp(:,ii), opt, 0);
        end
    end

    [sigslope.m, ~, sigslope.mdist] = CalcSignificanceBounds(dynht, Tsamp, opt, 0);

    %% BandPass filtering
    dynht = BandPass(dynht, opt.filt);
    if opt.filter_temp
        for zz=1:length(zsamp)
            Tsamp(:,zz) = BandPass(Tsamp(:,zz), opt.filt);
        end
    end

    if plotfigures
        axes(hax(1));
        PlotSpectrum(dynht);
        PlotSpectrum(Tsamp(:,5));
    end

    %% least squares and plot
    [infer_mode, infer_mode_error, corrcoeff, dof] = ...
        DoRegression(dynht', Tsamp', noise, sigslope, opt);

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
    flatbot.etopoDepth = 1;

    % Plot mode structure
    if plotfigures
        hax(2) = subplot(122);cla;
        try
            handles = PlotMode({modes; flatbot}, 1, 1, plotopt, hax(2));
            title('Test with synthetic time series');
            handles.hleg = legend('Location', 'SouthEast');
            handles.hleg.Box = 'off';
            xlim([-1.2 1.2]);
            ylim([-750 0]);
        catch ME
            title('No significant slopes found.');
        end
    end
end

function [dynht, Tsamp, Tmode] = CreateData(freq, ntimes, zfull, zsamp);

% ideal mode-1 shape
    Tmode(:,1) = -sin(pi*zfull/max(abs(zfull/16))  + 0* 2*pi * rand(size(zfull)));
    % ideal mode-2 shape
    Tmode(:,2) = -sin(2*pi*zfull/max(abs(zfull/16)) + 0* 2*pi * rand(size(zfull)));

    % generate a wave amplitude time series
    % one for each wave vertical mode
    % (of course these have to be at different frequencies)
    rng('shuffle');
    tseries = zeros([ntimes size(Tmode, 2)]);

    snr = sqrt(0.01);
    noiseSlope0 = -2;
    for ff=freq([2 4]) % mode 1
        rand1 = synthetic_timeseries_known_spectrum(ntimes, 1, 1, noiseSlope0);
        rand2 = synthetic_timeseries_known_spectrum(ntimes, 1, 1, noiseSlope0);

        tseries(:,1) = tseries(:,1) + rand1 + ...
            snr*max(rand1(:))./max(rand2(:)) .* rand2 .*sin(2*pi*ff*[0:ntimes-1]');
    end
    for ff=freq([1 3]) % mode 2
        rand1 = synthetic_timeseries_known_spectrum(ntimes, 1, 1, noiseSlope0);
        rand2 = synthetic_timeseries_known_spectrum(ntimes, 1, 1, noiseSlope0);

        tseries(:,2) = tseries(:,2) + rand1 + ...
            snr*max(rand1(:))./max(rand2(:)) .* rand2 .*sin(2*pi*ff*[0:ntimes-1]');
    end

    % Normalize so that temperatures are realistic enough that seawater toolbox works.
    tseries = bsxfun(@rdivide, tseries, max(tseries, [], 1));

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
    assert(isreal(dynht), 'Dynamic height is complex! Check temp amplitude');

end