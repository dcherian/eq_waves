
%% figure out data frequency slopes
% tao = ReadTaoTriton(11,4);

% figure;
% subplot(121)
% hspec = PlotSpectrum(tao.dht{11,4});
% PlotSpectrum(synthetic_timeseries_known_spectrum(length(tao.dht{11,4}),1,0.05,-2.4));
% uistack(hspec, 'top');

% subplot(122);
% hspec = PlotSpectrum(tao.dht{11,4});
% PlotSpectrum(synthetic_timeseries_known_spectrum(length(tao.dht{11,4}),1,0.05,-2));
% uistack(hspec, 'top');

%% Monte Carlo regression Slope

%slope = 0.15;
%intercept = 0.2;

function [mdist, m] = TestMC(spectralSlope, DoBandPass, makePlot, len)

    numMC = 1e4;

    if ~exist('spectralSlope', 'var'), spectralSlope = -3; end
    if ~exist('DoBandPass', 'var'), DoBandPass = 1; end
    if ~exist('makePlot', 'var'), makePlot = 0; end
    if ~exist('len', 'var'), len = 5e3; end

    xx = synthetic_timeseries_known_spectrum(len, 1, 1, spectralSlope);
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
        xx = cut_nan(xx - nanmean(xx));
    end

    tic;
    clear c m r
    for ii=1:numMC
        % xnoise = 0;
        % yy = slope * (xx + xnoise) + intercept;
        % ynoise = 0.2*max(abs(yy))*whitenoise(size(yy));
        % yy = yy + ynoise;

        yy = synthetic_timeseries_known_spectrum(len, ...
                                                 1, 1, spectralSlope);
        if DoBandPass
            yy = BandPass(yy, opt.filt);
            yy = cut_nan(yy - nanmean(yy));
        end

        [coeff,~,~,err] = dcregress(xx, yy, [], 0, 0, 0);
        c(ii) = coeff(1);
        m(ii) = coeff(2);

        % t hypothesis test
        % (b_1 - β_{10})/se(b_1)
        % Draper & Smith pg. 36 eqn. (1.4.9)
        mdist(ii) = (m(ii)-0)/err(2);

        rmat = corrcoef(cut_nan(xx), cut_nan(yy));
        r(ii)= rmat(1,2);
    end
    toc;

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
    else
        FitDistrib = allfitdist(mdist);
    end

    % filename = ['MC-slope-' num2str(spectralSlope)];
    % if DoBandPass
    %     filename = [filename '-BandPass'];
    % end


    % save(filename);