function [] = PlotBathy(etopo, lonrange, latrange, skip)

    x1 = min(find_approx(etopo.x, lonrange(1)));
    x2 = max(find_approx(etopo.x, lonrange(2)));

    xrange = [min([x1 x2]):skip:max([x1 x2])];
    yrange = [find_approx(etopo.y, latrange(1)):skip: ...
              find_approx(etopo.y, latrange(2))];

    etopo.z(etopo.z > 0) = NaN;

    if isempty(xrange)
        error('X range is empty!');
    end

    if isempty(yrange)
        error('Y range is empty');
    end

    contour(etopo.x(xrange), etopo.y(yrange), ...
            etopo.z(xrange, yrange)', 'k');
end