function [hplt] = PlotSpectrum(in, SubsetLength, hax)

    if ~exist('hax', 'var'), hax = gca; end
    if ~exist('SubsetLength', 'var'), SubsetLength = []; end

    [S, F] = GappySpectrum(in, SubsetLength);
    axes(hax); hold on; length(F)
    hplt = plot(F, S, '.-');
    hax.XScale = 'log';
    hax.YScale = 'log';

    xlabel('Frequency (cpd)')
    ylabel('Spectral density (m^2/s^2/cpd)')

    %[F,S] = mspec(in, [], 'cyclic');
    %plot(F, S);
