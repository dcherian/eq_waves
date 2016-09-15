function [] = PlotSpectrum(in, hax)

    if ~exist('hax', 'var'), hax = gca; end
    if size(in, 1) == 1, in = in'; end

    axes(hax); hold on;
    [F,S] = mspec(in, [], 'cyclic');
    plot(F,S);
    hax.XScale = 'log';
    hax.YScale = 'log';
