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

    in = in - nanmean(in);

    nans = isnan(in);
    edges = diff(nans);
    gapstart = find(edges == 1) + 1;
    gapend = find(edges == -1);

    if isnan(in(1)), gapstart = [1 gapstart]; end
    if isnan(in(end)), gapend(end+1) = length(in); end

    if isempty(gapstart) & isempty(gapend) ...
            & isequal(nans, zeros(size(nans)))
        % input series has no gaps
        gapstart = length(in) + 1;
        gapend = gapstart;
    end

    assert(length(gapstart) == length(gapend), ...
           ['FilterSeries: gapstart and gapend are not same ' ...
            'size.']);

    start = 1;
    out = nan(size(in));
    for ii=1:length(gapstart)
        range = start:gapstart(ii)-1;

        if isempty(range) | (length(range) < opt.N), continue; end

        out(range) = smooth_1d(in(range), opt.N, ...
                               opt.halfdef, opt.window);

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
        PlotSpectrum(smoothed);
        linex(1/opt.N);
    end
end
