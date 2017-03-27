%% Monte Carlo regression Slope

function [mdist, m, r, rdist] = TestMC(spectralSlope, DoBandPass, ...
                                       makePlot, len, yin, spectralAmp)

    numMC = 5e3;

    if ~exist('spectralSlope', 'var'), spectralSlope = -3; end
    if ~exist('DoBandPass', 'var'), DoBandPass = 1; end
    if ~exist('makePlot', 'var'), makePlot = 0; end
    if ~exist('len', 'var'), len = 5e3; end
    if ~exist('spectralAmp', 'var'), spectralAmp = 1; end
    if ~exist('yin', 'var')
        yin = synthetic_timeseries_known_spectrum(len, 1, spectralAmp, spectralSlope);
    end

    opt = DefaultOptions;
    if DoBandPass
        warning(['Running BandPass filter with default options.', ...
                 'cutoff = [' num2str(opt.filt.cutoff, '%.2f ') ']']);
        if makePlot
            figure;
            subplot(211); plot(yin);
        end
        yin = BandPass(yin, opt.filt);
        if makePlot
            subplot(212); plot(yin);
        end
    end

    tic;
    clear c m r

    c = nan([numMC 1]);
    m = c; r = c; mdist = c; rdist = c;

    if all(isnan(yin)), return; end

    xxmat = synthetic_timeseries_known_spectrum([len numMC], 1, ...
                                                spectralAmp, spectralSlope);

    c = nan([numMC 1]);
    m = c; mdist = c; r = c;

    parfor ii=1:numMC
        % xnoise = 0;
        % xx = slope * (yin + xnoise) + intercept;
        % ynoise = 0.2*max(abs(xx))*whitenoise(size(xx));
        % xx = xx + ynoise;

        xx = xxmat(:,ii);
        if size(xx) ~= size(yin), xx = xx'; end

        if DoBandPass
            xx = BandPass(xx, opt.filt);
        end

        if size(yin) ~= size(xx)
            xx = xx';
        end

        [coeff,~,~,err] = dcregress(xx, yin, [], 0, 0, 0);
        c(ii) = coeff(1);
        m(ii) = coeff(2);

        % t hypothesis test
        % (b_1 - β_{10})/se(b_1)
        % Draper & Smith pg. 36 eqn. (1.4.9)
        mdist(ii) = (m(ii)-0)/err(2);

        mask = isnan(yin) | isnan(xx);
        rmat = corrcoef(yin(~mask), xx(~mask));
        r(ii)= rmat(1,2);
    end

    % the variable rdist == w (bendat piersol eqn. 4.58) is
    % normally distributed with mean = 0 and σ²=1/(dof-3)
    rdist = 1/2*log((1+r)./(1-r));

    if makePlot
        nbins = 40;
        figure;
        subplot(221)
        plot(xx,yin, '*');
        title(['Spectral slope = ' num2str(spectralSlope)]);
        beautify;

        subplot(222)
        histogram(r, nbins, 'Normalization', 'countdensity', ...
                  'EdgeColor', 'none');
        linex(calc95(r), [], 'r');
        title('corr. coeff');
        beautify;

        subplot(223);
        hsl = histogram(m, nbins, 'Normalization', 'countdensity', ...
                        'EdgeColor', 'none');
        xbin = hsl.BinEdges;
        %linex(slope, [], 'k')
        %linex([coeff(2) coeff(2)-conf(2) coeff(2)+conf(2)]);
        linex(calc95(m), [], 'r');
        title('slope');
        beautify;

        subplot(224);
        histogram(c, nbins, 'Normalization', 'countdensity', 'EdgeColor', 'none');
        %linex(intercept, [], 'k');
        %linex([coeff(1) coeff(1)-conf(1) coeff(1)+conf(1)]);
        linex(calc95(c), [], 'r');
        title('intercept');
        beautify;

        FitDistrib = allfitdist(mdist, 'PDF');
        title(['Spectral Slope = ' num2str(spectralSlope) ' | DoBandPass ' ...
               '= ' num2str(DoBandPass)]);
    end

    % filename = ['MC-slope-' num2str(spectralSlope)];
    % if DoBandPass
    %     filename = [filename '-BandPass'];
    % end

    % save(filename);

    toc;
end