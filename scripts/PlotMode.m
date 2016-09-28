function [] = PlotMode(modename, mm, nn, hax)

    load([modename '.mat']);

    linewidth = 2;

    if modes.lon(mm) > 0
        lonstr = 'E';
    else
        lonstr = 'W';
    end

    if modes.lat(nn) < 0
        latstr = 'S';
    else
        latstr = 'N';
    end

    if ~exist('hax', 'var')
        figure;
        set(gcf, 'Position', [520 118 650 680]);
    else
        axes(hax);
    end

    hold on;
    % inferred mode from TAO data
    errorbar(modes.InferredMode{mm,nn}, ...
             modes.depth{mm,nn} * -1, ...
             modes.InferredModeError{mm,nn}, ...
             'horizontal', ...
             'Marker', '.', 'MarkerSize', 12, ...
             'LineStyle', 'none', 'LineWidth', linewidth, ...
             'DisplayName', 'T_{TAO/TRITON}');

    % 0 mean flow mode
    for ii=1:2
        plot(squeeze(modes.IdealTempMode(mm,nn,:,ii)) ...
             ./ max(abs(modes.IdealTempMode(mm,nn,:,ii))), ...
             modes.zTmode * -1, 'LineWidth', linewidth, ...
             'DisplayName', ['T_{bc' num2str(ii) '}']);
    end

    title(sprintf('(%3d%s, %1d%s) | %s', abs(modes.lon(mm)), lonstr, ...
                  abs(modes.lat(nn)), latstr, opt.filt.window));
    ylim([-700 0]);
    xlim([-0.2 1.3]);
    hleg = legend('Location', 'SouthEast');
    hleg.Box = 'off';
    linex(0);
end
