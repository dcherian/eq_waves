%% Monte Carlo regression Slope

%slope = 0.15;
%intercept = 0.2;

function [mdist, m, r, rdist] = TestMC(spectralSlope, DoBandPass, ...
                                       makePlot, len, xx, spectralAmp)

    numMC = 5e3;

    if ~exist('spectralSlope', 'var'), spectralSlope = -3; end
    if ~exist('DoBandPass', 'var'), DoBandPass = 1; end
    if ~exist('makePlot', 'var'), makePlot = 0; end
    if ~exist('len', 'var'), len = 5e3; end
    if ~exist('spectralAmp', 'var'), spectralAmp = 1; end
    if ~exist('xx', 'var')
        xx = synthetic_timeseries_known_spectrum(len, 1, spectralAmp, spectralSlope);
    end

    opt = DefaultOptions;
    if DoBandPass
        if makePlot
            figure;
            subplot(211); plot(xx);
        end
        xx = BandPass(xx, opt.filt);
        if makePlot
            subplot(212); plot(xx);
        end
    end

    tic;
    clear c m r

    c = nan([numMC 1]);
    m = c; r = c; mdist = c; rdist = c;
    for ii=1:numMC
        % xnoise = 0;
        % yy = slope * (xx + xnoise) + intercept;
        % ynoise = 0.2*max(abs(yy))*whitenoise(size(yy));
        % yy = yy + ynoise;

        yy = synthetic_timeseries_known_spectrum(len, 1, ...
                                                 spectralAmp, spectralSlope);
        if size(yy) ~= size(xx), yy = yy'; end

        if DoBandPass
            yy = BandPass(yy, opt.filt);
        end

        if size(xx) ~= size(yy)
            yy = yy';
        end

        [coeff,~,~,err] = dcregress(xx, yy, [], 0, 0, 0);
        c(ii) = coeff(1);
        m(ii) = coeff(2);

        % t hypothesis test
        % (b_1 - β_{10})/se(b_1)
        % Draper & Smith pg. 36 eqn. (1.4.9)
        mdist(ii) = (m(ii)-0)/err(2);

        mask = isnan(xx) | isnan(yy);
        rmat = corrcoef(xx(~mask), yy(~mask));
        r(ii)= rmat(1,2);
    end

    % the variable rdist == w (bendat piersol eqn. 4.58) is
    % normally distributed with mean = 0 and σ²=1/(dof-3)
    rdist = 1/2*log((1+r)./(1-r));

    %[coeff, conf, dof] = dcregress(xx, yy, length(yy)-2, 0, 0);
    % disp(['actual: slope = ' num2str(slope, '%.3f') ' | intercept = ' ...
    %       num2str(intercept, '%.3f')]);
    % disp(['dcregress: slope = [' num2str([coeff(2)-conf(2) coeff(2)+conf(2)], '%.3f ') ...
    %       '] | intercept = [' num2str([coeff(1)-conf(1) coeff(1)+conf(1)], ...
    %                                   '%.3f ') ']']);

    if makePlot
        nbins = 40;
        figure;
        subplot(221)
        plot(xx,yy, '*');
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