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

    if opt.debugflag
        figure;
        hax = gca;
        hold on;
    end

    [gapstart, gapend] = FindGaps(in);

    if strcmpi(opt.window, 'butterworth')
        % butterworth requires (desired freq)/(sampling freq/2)
        % the factor of 4 is because
        % a) 1./opt.cutoff = (desired freq)/2
        %      as designed for other windows.
        % b) (sampling freq/2) = 1/2 cpd.
        [b,a] = butter(1, sort(2./opt.cutoff/(1/2)), 'bandpass');
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
            segment = in(range) - nanmean(in(range));
            out(range) = filter(b, a, segment);

            if opt.debugflag
                % plot filtered time series and step response
                axes(hax);
                plot(range, out(range), '-', 'Color', hax.ColorOrder(1,:))
                [h,t] = stepz(b, a, length(range), 1);
                h = h./max(abs(h));
                plot(range, (h + flip(h)).*nanmax(segment)/2, '-', ...
                     'Color', hax.ColorOrder(3,:))
                keyboard;
            end

            % NaN out contaminated edges
            out(range(1) : min(range(1)+60, length(in))) = NaN;
            out( max(1,range(end)-60-1) : range(end)) = NaN;

            if opt.debugflag
                plot(range, out(range), 'k-');
            end
        else
            out(range) = smooth_1d(in(range) - nanmean(in(range)), opt.N, ...
                                   opt.halfdef, opt.window);

            % NaN out contaminated edges
            out(range(1) : min(range(1)+ceil(opt.N)+1, length(in))) = NaN;
            out( max(1,range(end)-ceil(opt.N)-1) : range(end)) = NaN;
        end

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
