function [out] = BandPass(in, opt)

    if ~strcmpi(opt.window, 'butterworth')
        opt.N = opt.cutoff(1); % first higher freq cutoff
        hi = FilterSeries(in, opt);
        opt.N = opt.cutoff(2);
        low = FilterSeries(in, opt);

        out = hi - low;
    else
        opt.N = max(opt.cutoff); % only used to chuck points near
                                 % end of segment
        out = FilterSeries(in, opt);
    end

    if opt.debugflag
        figure;
        hin = PlotSpectrum(in);
        hin.LineWidth = 1;
        hold on;
        hout = PlotSpectrum(out*1);

        if ~strcmpi(opt.window, 'butterworth')
            PlotSpectrum(low*1);
            PlotSpectrum(hi*1);
        end
        legend('input', 'out', 'low-freq', 'high-freq');
        uistack(hout, 'top');
        linex(1./opt.cutoff);
        title(num2str(opt.cutoff));
    end
end
