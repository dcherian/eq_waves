function [NoiseAmp, NoiseSlope] = EstimateNoiseSpectrum(in, opt, ...
                                                      plotflag, hax)
    % Estimate noise spectrum as S = (NoiseAmp) * (freq)^(NoiseSlope)
    % Do this by fitting a line in
    %    (freq < low-frequency cutoff).
    % Optionally, use frequency information in
    %    (freq > high frequency cutoff).

    useHiFreq = 0;
    if ~exist('plotflag', 'var'), plotflag = 0; end
    if ~exist('opt', 'var'), opt = []; end

    SubsetLength = 256;
    [S, freq] = GappySpectrum(in, SubsetLength);

    % if opt is not provided, use full spectrum for fit
    if isempty(opt)
        indfreqlo = length(S);
        indfreqhi = length(S);
        freqreg = freq;
    else
        indfreqlo = find(freq > min(2./opt.filt.cutoff), 1, 'first')-1;
        if useHiFreq
            indfreqhi = find(freq < max(2./opt.filt.cutoff), 1, 'last')+1;
        else
            indfreqhi = length(S);
        end
        assert(freq(indfreqlo) <= min(2./opt.filt.cutoff));
        assert(freq(indfreqhi) >= max(2./opt.filt.cutoff));

        S(indfreqlo:indfreqhi) = NaN;
        freqreg = freq;
        freqreg(indfreqlo:indfreqhi) = NaN;
        % S(1:4) = NaN;
        % freqreg(1:4) = NaN;
    end

    [coeff] = dcregress(log(freqreg), log(S), NaN, 0, 0, 0, 0);
    NoiseAmp = exp(coeff(1));
    NoiseSlope = coeff(2);

    if plotflag
        xline = log(freq);
        yline = coeff(2)*xline + coeff(1);

        tseries = synthetic_timeseries_known_spectrum(...
            length(in), 1, NoiseAmp, NoiseSlope);
        % insert gaps at the same time instants as in input so that
        % PlotSpectrum gets quivalent input to work with
        tseries(isnan(in)) = NaN;

        if ~exist('hax', 'var')
            figure;
        else
            axes(hax);
        end
        hh = PlotSpectrum(in, SubsetLength);
        hh.DisplayName = 'Input';

        if ~isempty(opt)
            hh = PlotSpectrum(BandPass(in, opt.filt), SubsetLength);
            hh.DisplayName = 'Filtered input';
        end

        hplt = PlotSpectrum(tseries, SubsetLength);
        hplt.DisplayName = 'Generated noise';

        plot(exp(xline(1:indfreqlo)), exp(yline(1:indfreqlo)), 'k', ...
             'DisplayName', 'Straight line fit');
        plot(exp(xline(indfreqlo+1:indfreqhi)), exp(yline(indfreqlo+1:indfreqhi)), ...
             'k--', 'HandleVisibility', 'off');
        plot(exp(xline(indfreqhi+1:end)), exp(yline(indfreqhi+1:end)), ...
             'k--', 'HandleVisibility', 'off');
        limy = ylim;
        if ~isempty(opt)
            hpt = patch(2./[opt.filt.cutoff(2) opt.filt.cutoff(2) ...
                            opt.filt.cutoff(1) opt.filt.cutoff(1) ...
                            opt.filt.cutoff(2)], ...
                        [min(ylim) max(ylim) max(ylim) ...
                         min(ylim) min(ylim)], ...
                        'k', 'EdgeColor', 'none', 'FaceAlpha', 0.1, ...
                        'DisplayName', 'Filter band');
        end
        legend('Location', 'NorthEast');
        axis('square');
        beautify;
        grid on;
    end
end