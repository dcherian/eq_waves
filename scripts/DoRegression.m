function [infer_mode, infer_mode_error, corrcoeff, intercept, stdError] = ...
        DoRegression(dhtinput, Tinput, opt)

    if ~exist('opt', 'var') | ~isfield(opt, 'debug')
        opt.debug = 0;
    end

    nz = size(Tinput, 1);

    infer_mode = nan([nz 2]);
    infer_mode_error = nan([nz 2]);
    corrcoeff = nan([nz 1]);
    dof = nan([nz 1]);
    stdError = dof;

    for ii = 1:nz
        dht = dhtinput;
        T = Tinput(ii,:);

        % find all nan's in both datasets
        mask = ~(isnan(dht) | isnan(T));
        T(~mask) = NaN;
        dht(~mask) = NaN;

        if all(mask == 0) % all NaNs
            continue;
        end

        if opt.debug
            figure(opt.hdbg);
            hdbgax2 = subplot(212);
            plot(dht(mask)); hold on;
            plot(T(mask));
            title('filtered time series for regression');
            keyboard;
        end

        % use 0.5cm error on dynamic height (units cm)
        % use 0.01 C error on temperature
        % (Tom's suggestion).
        if opt.TagainstDHT
            rrwtls = mf_wtls(dht', T', 0.1, 0.01, 0);

            [rrols([3 1]), rrols([4 2]), ~, stderr] = ...
                dcregress(dht', T'-nanmean(T), [], 0,0,0);
        else
            rrwtls = mf_wtls(T', dht', 0.01, 0.5, 0);

            [rrols([3 1]), rrols([4 2]), ~] = ...
                dcregress(T'-nanmean(T), dht', [], 0,0,0);
        end

        infer_mode(ii,1:2) = [rrols(1) rrwtls(1)];
        infer_mode_error(ii,1:2) = [1.96*stderr(2)*opt.SlopeSigma rrwtls(2)];
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
    end
end
