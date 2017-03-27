function [infer_mode, infer_mode_error, corrcoeff, intercept, stdError] = ...
        DoRegression(dhtinput, Tinput, noise, sigslope, opt)

    if ~exist('opt', 'var') | ~isfield(opt, 'debug')
        opt.debug = 0;
    end

    nz = size(Tinput, 1);

    infer_mode = nan([nz 2]);
    infer_mode_error = nan([nz 2]);
    corrcoeff = nan([nz 1]);
    dof = nan([nz 1]);
    stdError = dof;

    if opt.debugNullRegression, hnull = figure; end

    % 1. generate numMC noise time series
    if ~opt.TagainstDHT
        noisemat = synthetic_timeseries_known_spectrum( ...
            [length(dhtinput) opt.numMC], 1, noise.amp, noise.slope);
    end

    for ii = 1:nz
        dht = dhtinput - nanmean(dhtinput);
        T = Tinput(ii,:) - nanmean(Tinput(ii,:));

        % find all nan's in both datasets
        mask = ~(isnan(dht) | isnan(T));
        T(~mask) = NaN;
        dht(~mask) = NaN;

        if all(mask == 0) % all NaNs
            continue;
        end

        if opt.debugRegression
            figure;
            plot(dht(mask)); hold on;
            plot(T(mask));
            title(['filtered time series for regression' ...
                   ' | depth level = ' num2str(ii)]);
        end

        % use 0.5cm error on dynamic height (units cm)
        % use 0.01 C error on temperature
        % (Tom's suggestion).
        if opt.TagainstDHT
            rrwtls = mf_wtls(dht', T', 0.1, 0.01, 0);

            [rrols([3 1]), rrols([4 2]), ~, stderr] = ...
                dcregress(dht', T', NaN, 0, opt.debugRegression, 0);

            if abs(rrols(1)) > sigslope.m(ii)
                noisemat = synthetic_timeseries_known_spectrum( ...
                    [length(dhtinput) opt.numMC], 1, ...
                    noise.amp(ii), noise.slope(ii));

                m = MonteCarloRegressionSlopeDistrib( ...
                    dht, T, rrols(1), rrols(3), noisemat, opt);

                rrols(1) = nanmean(m);
            else
                rrols(1) = NaN;
                m = [1 1]*NaN;
            end
        else
            rrwtls = mf_wtls(T', dht', 0.01, 0.5, 0);

            [rrols([3 1]), rrols([4 2]), ~, stderr] = ...
                dcregress(T', dht', [], ...
                          0, opt.debugRegression, 0);

            if abs(rrols(1)) > sigslope.m(ii)
                noisemat = synthetic_timeseries_known_spectrum( ...
                    [length(dhtinput) opt.numMC], 1, ...
                    noise.amp(ii), noise.slope(ii));

                m = MonteCarloRegressionSlopeDistrib( ...
                    T, dht, rrols(1), rrols(3), noisemat, opt);
                rrols(1) = nanmean(m);
            else
                rrols(1) = NaN;
                m = [1 1]*NaN;
            end
        end

        infer_mode(ii,1:2) = [rrols(1) rrwtls(1)];
        % symmetric 95% confidence interval from distribution of
        % regression slopes from Monte Carlo simulation.
        if isnan(rrols(1))
            infer_mode_error(ii,1) = NaN;
        else
            infer_mode_error(ii,1) = mean(abs(rrols(1) - calc95(m)));
        end
        infer_mode_err(ii,2) = rrwtls(2);
        intercept(ii,:) = [rrols(3) rrols(4)];
        stdError(ii) = stderr(2);

        % "Since the slope from the GMFR is simply a ratio of
        % variances, it is ``transparent'' to the
        % determination of correlation coefficients or
        % coefficients of determination. It is these
        % correlations, not the slope of the line that test
        % the strength of the linear relationship between the
        % two variables" - Emery & Thompson (2001), pg. 249
        corrcoeff(ii) = min(min( ...
            corrcoef(dht(mask)', T(mask)')));

        if ~strcmpi(opt.name, 'montecarlo')
            if abs(corrcoeff(ii)) <= opt.corr_sig
                % corrcoeff = 0 means insignificant
                % corrcoeff = NaN means no data.
                corrcoeff(ii) = 0;
                infer_mode(ii, 2) = NaN;
                infer_mode_error(ii, 2) = NaN;
            end
        end

        if opt.debugNullRegression
            figure(hnull);
            subplot(2, ceil(nz/2), ii);
            histogram(sigslope.mdist{ii});
            if ~isnan(rrols(1));
                linex(rrols(1));
                calc95(m);
            end
            title(num2str(ii));
        end
    end
end

function [m] = MonteCarloRegressionSlopeDistrib(x, y0, slope, intercept, noisemat, opt)
% Monte Carlo estimation of bounds on regression slope
% presumes that provided slope is significant.

    ticstart = tic;
    disp('Monte Carlo error estimation');
    m = nan([1 size(noisemat,2)]);
    parfor mc=1:size(noisemat, 2)
        noisevec = BandPass(noisemat(:,mc), opt.filt);
        y =  slope*x' + intercept + noisevec;
        coeff = dcregress(x, y, NaN, 0, 0, 0, 0);
        c(mc) = coeff(1);
        m(mc) = coeff(2);
    end
    toc(ticstart);

    % figure;
    % hold on; histogram(m); linex(slope);

    % keyboard;
    % figure;
    % PlotSpectrum(slope*x+intercept);
    % PlotSpectrum(noisemat(:,mc));
    % PlotSpectrum(BandPass(noisemat(:,mc), opt.filt));
    % PlotSpectrum(y);
    % PlotSpectrum(BandPass(y0, opt.filt));
    % legend('slope*x', 'noise', 'noise filt', 'y', 'T');
end