function [hplt] = PlotSpectrum(in, hax)

    if ~exist('hax', 'var'), hax = gca; end
    if size(in, 1) == 1, in = in'; end

    data = ~isnan(in);
    cc = bwconncomp(data);

    maxdataidx = 1;
    for ii=2:cc.NumObjects
        if length(cc.PixelIdxList{ii}) > ...
                length(cc.PixelIdxList{maxdataidx})
            maxdataidx = ii;
        end
    end

    axes(hax); hold on;
    [S, F] = spectrum_band_avg(in(cc.PixelIdxList{maxdataidx}), ...
                              1, 5, 'hann', 0);
    %[F,S] = mspec(in, [], 'cyclic');

    hplt = plot(F,S);
    hax.XScale = 'log';
    hax.YScale = 'log';
