% Filter a gappy time series using smooth_1d.m
function [out] = FilterSeries(in, opt)

    debugflag = 0;

    if ~isfield(opt, 'window')
        opt.window = 'rect';
        disp('Using default rectangular window');
    end

    if ~isfield(opt, 'halfdef')
        opt.halfdef = 'power';
    end

    [gapstart, gapend] = FindGaps(in);

    if strcmpi(opt.window, 'butterworth')
        % butterworth requires (desired freq)/(sampling freq/2)
        % the factor of 4 is because
        % a) 1./opt.cutoff = (desired freq)/2
        %      as designed for other windows.
        % b) (sampling freq/2) = 1/2 cpd.
        [b,a] = butter(4, sort(2./opt.cutoff/(1/2)), 'bandpass');
    end

    start = 1;
    out = nan(size(in));
    for ii=1:length(gapstart)
        range = start:gapstart(ii)-1;

        if isempty(range) | (length(range) < opt.N), continue; end

        % remove mean for each section and the filter.
        % This makes the filtering work better
        % and analysis less sensitive to choice of window.
        if strcmpi(opt.window, 'butterworth')
            out(range) = filter(b, a, in(range)-nanmean(in(range)));
        else
            out(range) = smooth_1d(in(range) - nanmean(in(range)), opt.N, ...
                                   opt.halfdef, opt.window);
        end

        % NaN out contaminated edges
        out(range(1) : min(range(1)+ceil(opt.N)+1, length(in))) = NaN;
        out( max(1,range(end)-ceil(opt.N)-1) : range(end)) = NaN;

        % set start for next iteration to be end of current gap.
        start = gapend(ii)+1;
    end

    if debugflag
        figure;
        PlotSpectrum(in);
        hold on;
        PlotSpectrum(out);
        linex(1/opt.N);
    end
end
