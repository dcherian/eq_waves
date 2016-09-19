function [out] = BandPass(in, opt)

    % first higher freq cutoff
    opt.N = opt.cutoff(1);
    hi = FilterSeries(in, opt);
    opt.N = opt.cutoff(2);
    low = FilterSeries(in, opt);

    out = hi - low;

    if opt.debugflag
        figure;
        hin = PlotSpectrum(in);
        hin.LineWidth = 1;
        hold on;
        PlotSpectrum(low*1);
        PlotSpectrum(hi*1);
        hout = PlotSpectrum(out*1);
        legend('input', 'low-freq', 'high-freq', 'out');
        linex(1./opt.cutoff);
        title(num2str(opt.cutoff));
    end
end
