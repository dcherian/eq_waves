function [infer_mode, infer_mode_error, corrcoeff, dof] = ...
        DoRegression(dhtinput, Tinput, opt)

    if ~exist('opt', 'var') | ~isfield(opt, 'debug')
        opt.debug = 0;
    end

    nz = size(Tinput, 1);

    infer_mode = nan([nz 1]);
    infer_mode_error = infer_mode;
    corrcoeff = infer_mode;
    dof = infer_mode;

    for ii = 1:nz
        dht = dhtinput;
        T = Tinput(ii,:);

        % find all nan's in both datasets
        mask = ~(isnan(dht) | isnan(T));
        T(~mask) = NaN;
        dht(~mask) = NaN;

        if all(mask == 0)
            infer_mode(ii) = NaN;
            infer_mode_error(ii,1) = NaN;
            continue;
        end

        % regress to find mode shape
        % infer_mode(ii) = dhtfilt(mask)' \ T(mask)';
        % [infer_mode(ii), bint] = ...
        %     regress(T(mask)', dhtfilt(mask)');
        % infer_mode_error(ii,1) = bint(2) - infer_mode(ii);

        if opt.debug
            figure(opt.hdbg);
            hdbgax2 = subplot(212);
            plot(dht(mask)); hold on;
            plot(T(mask));
            title('filtered time series for regression');
            keyboard;
        end

        if opt.filter_temp
            dof(ii) = floor(min([calcdof(dht) ...
                                calcdof(T)]) * 2/pi);
        else
            dof(ii) = calcdof(dht) * 2/pi;
        end

        Nsamp = ceil(length(mask) / dof(ii)) + 1;

        if length(cut_nan(dht(1:Nsamp:end))) <= 2 | ...
                length(cut_nan(isnan(T(1:Nsamp:end)))) <= 2
            continue;
        end

        rrwtls = mf_wtls(dht(1:Nsamp:end)', T(1:Nsamp:end)', ...
                         5e-2, nanstd(T), 0);

        [rrols([3 1]), rrols([4 2]), ~] = ...
            dcregress(dht', T'-nanmean(T), dof(ii), [], 0);

        % "Since the slope from the GMFR is simply a ratio of
        % variances, it is ``transparent'' to the
        % determination of correlation coefficients or
        % coefficients of determination. It is these
        % correlations, not the slope of the line that test
        % the strength of the linear relationship between the
        % two variables" - Emery & Thompson (2001), pg. 249
        corrcoeff(ii) = min(min( ...
            corrcoef(dht(mask)', T(mask)')));

        if 0 %abs(corrcoeff(ii)) <= corr_sig(dof(ii)-2, 0.95)
            % 0 means insignificant, NaN means no data.
            corrcoeff(ii) = 0;
            infer_mode(ii) = NaN;
            infer_mode_error(ii) = NaN;
        else
            infer_mode(ii,1:2) = [rrols(1) rrwtls(1)];
            infer_mode_error(ii,1:2) = [rrols(2) rrwtls(2)];
        end
    end
end
