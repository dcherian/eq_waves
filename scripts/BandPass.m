function [out] = BandPass(in, opt)

    % first higher freq cutoff
    opt.N = opt.cutoff(1);
    hi = FilterSeries(in, opt);
    opt.N = opt.cutoff(2);
    low = FilterSeries(in, opt);

    out = hi - low;

    if opt.debugflag
        figure;
        PlotSpectrum(in);
        hold on;
        PlotSpectrum(low*1);
        PlotSpectrum(hi*1);
        PlotSpectrum(out*1);
        linex(1./opt.cutoff);
    end
end
