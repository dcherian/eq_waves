% Test effect of filtering on correlation
function [] = CorrFilter(spectralSlope, hax)

    if ~exist('spectralSlope', 'var'), spectralSlope = 0; end
    if ~exist('hax', 'var'), hax = []; end

    opt = DefaultOptions;
    xx = synthetic_timeseries_known_spectrum(200, 1, 1, spectralSlope);

    [c, lags] = xcorr(xx, 'coeff');

    if isempty(hax)
        figure;
        hax = gca;
    else
        axes(hax);
    end

    plot(lags, c)
    hold on;
    [c, lags] = xcorr(cut_nan(BandPass(xx, opt.filt)), 'coeff');
    plot(lags, c)
    linex(0); liney(0);
    title(['SpectralSlope = ' num2str(spectralSlope)]);
    xlabel('Lag (points)');
    ylabel('Correlation coeff');
    legend('Unfiltered', 'filtered', 'Location', 'NorthWest');
    beautify;