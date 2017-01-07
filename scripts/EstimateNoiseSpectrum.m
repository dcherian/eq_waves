function [NoiseAmp, NoiseSlope] = EstimateNoiseSpectrum(in, opt, ...
                                                      plotflag, hax)
    % Estimate noise spectrum as S = (NoiseAmp) * (freq)^(NoiseSlope)
    % Do this by fitting a line in frequency band up to low
    % frequency cutoff

    if ~exist('plotflag', 'var'), plotflag = 0; end

    [S, freq] = GappySpectrum(in, 365);
    indfreq = find(freq > min(2./opt.filt.cutoff), 1, 'first')-1;
    assert(freq(indfreq) < min(2./opt.filt.cutoff));

    [coeff] = dcregress(log(freq(1:indfreq)), ...
                        log(S(1:indfreq)), NaN, 0, 0, 0, 0);
    NoiseAmp = exp(coeff(1));
    NoiseSlope = coeff(2);

    if plotflag
        xline = log(freq);
        yline = coeff(2)*xline + coeff(1);

        tseries = synthetic_timeseries_known_spectrum(...
            365, 1, NoiseAmp, NoiseSlope);
        if ~exist('hax', 'var')
            figure;
        else
            axes(hax);
        end
        PlotSpectrum(in, 365);
        hplt = PlotSpectrum(tseries); %hplt.Color = [1 1 1]*0.25;
        plot(exp(xline(1:indfreq)), exp(yline(1:indfreq)), 'k');
        plot(exp(xline(indfreq+1:end)), exp(yline(indfreq+1:end)), 'k--');
        limy = ylim;
        hpt = patch(2./[opt.filt.cutoff(2) opt.filt.cutoff(2) ...
                        opt.filt.cutoff(1) opt.filt.cutoff(1) ...
                        opt.filt.cutoff(2)], ...
                    [min(ylim) max(ylim) max(ylim) ...
                     min(ylim) min(ylim)], ...
                    'k', 'EdgeColor', 'none', 'FaceAlpha', 0.1);
        legend('Input', 'Generated noise', ...
               'Straight line fit', 'Straight line extrapolated', ...
               'Filter band', 'Location', 'SouthWest');
        axis('square');
        beautify;
        grid on;
    end
end