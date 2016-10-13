function [handles] = PlotMode(modename, mm, nn, plotopt, hax)

    linewidth = 1;
    capwidth = 12;
    linestylemode = {'--'; '-'; '-.'}; % line style for theoretical modes.

    stack = dbstack;
    if length(stack) > 1
        if strcmpi(stack(2).name, 'PlotModeMap')
            capwidth = 50;
        end
    end

    if ischar(modename)
        load([modename '.mat']);
    else
        % PlotModeMap is calling
        modes = modename{1};
        try
            % TestInference won't have this.
            data = modename{2};
        catch ME
        end
    end

    if ~exist('plotopt', 'var')
        plotopt.nmode = [1 2];
        plotopt.plotstd = 0;
        plotopt.plotcorr = 1;
    else
        if ~isfield(plotopt, 'nmode')
            plotopt.nmode = [1 2];
        end
        if ~isfield(plotopt, 'plotstd')
            plotopt.plotstd = 0;
        end
        if ~isfield(plotopt, 'plotcorr')
            plotopt.plotcorr = 1;
        end
    end

    if ~exist('hax', 'var')
        figure;
        set(gcf, 'Position', [520 118 650 680]);
        hax = gca;
        providedHax = 0;

        try % will fail for TestInference where this is N/A
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
        catch ME
        end
    else
        providedHax = 1;
    end

    hold on;
    % inferred mode from TAO data
    handles.herr = supererr(hax, modes.InferredMode{mm,nn}, ...
                            modes.depth{mm,nn} * -1, ...
                            abs(modes.InferredModeError{mm,nn}), ...
                            [], 'rI', capwidth, ...
                            'LineWidth', linewidth, ...
                            'DisplayName', 'T_{TAO/TRITON}');
    handles.herr = cut_nan(handles.herr(:,1));

    % 0 mean flow mode
    for ii=plotopt.nmode
        handles.hmode(ii) = ...
            plot(hax, squeeze(modes.IdealTempMode(mm,nn,:,ii)) ...
                 ./ max(abs(modes.IdealTempMode(mm,nn,:,ii))), ...
                 modes.zTmode * -1, 'LineWidth', linewidth, ...
                 'Color', 'k', 'LineStyle', linestylemode{ii}, ...
                 'DisplayName', ['T_{bc' num2str(ii) '}']);
    end

    % temp std
    if plotopt.plotstd
        handles.hstd = ...
            plot(hax, data.Tstd{mm,nn}./max(data.Tstd{mm,nn}), ...
                 data.depth{mm,nn} * -1,'k', 'LineWidth', linewidth, ...
                 'DisplayName', 'T_{std}');
    end

    % correlation coefficient
    if plotopt.plotcorr
        handles.hcorr = ...
            plot(hax, modes.corr{mm,nn}, modes.depth{mm,nn} * -1, ...
                 'b.', 'MarkerSize', 12, 'DisplayName', ...
                 'Corr. Coeff.');
    end

    if ~providedHax
        if isfield(modes, 'lon')
            title(sprintf('(%3d%s, %1d%s)', abs(modes.lon(mm)), lonstr, ...
                          abs(modes.lat(nn)), latstr));
        end
        ylim([-700 0]);
        xlim([-0.2 1.3]);
        handles.hleg = legend('Location', 'SouthEast');
        handles.hleg.Box = 'off';
        set(gca, 'XAxisLocation','Top');
        set(gcf, 'Renderer', 'painters');
    end

    handles.h0 = plot(hax, [0 0], hax.YLim, '--', 'Color', [1 1 1]*0.6, ...
                      'LineWidth', 1, 'LegendDisplay', 'off');

    uistack(handles.h0, 'bottom');
    uistack(handles.herr, 'top');
end
