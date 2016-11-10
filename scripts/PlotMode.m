%  [handles] = PlotMode(modename, mm, nn, plotopt, hax)
function [handles] = PlotMode(modename, mm, nn, plotopt, hax)

    if ~exist('mm', 'var'), mm = 1; end
    if ~exist('nn', 'var'), nn = 1; end

    linewidth = 1;
    capwidth = 12;
    linestylemode = {'--'; '-'; '-.'}; % line style for theoretical modes.

    red = [227,74,51]/255;
    green = [35,132,67]/255;

    stack = dbstack;
    if length(stack) > 1
        if strcmpi(stack(2).name, 'PlotModeMap')
            capwidth = 50;
        end
    end

    load('flat-bot-modes.mat');
    if ischar(modename)
        load([modename '.mat']);
    else
        if iscell(modename) & all(size(modename) == [2 1])
            % TestInference calling
            modes = modename{1};
            flatbot = modename{2};
        else
            modes = modename;
        end
    end

    if ~exist('plotopt', 'var')
        plotopt.nmode = [1 2];
        plotopt.plotstd = 0;
        plotopt.plotcorr = 1;
        plotopt.ploterr = 1;
        plotopt.plotOLS = 1;
        plotopt.plotWTLS = 1;
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
        if ~isfield(plotopt, 'ploterr')
            plotopt.ploterr = 1;
        end
        if ~isfield(plotopt, 'plotOLS')
            plotopt.plotOLS = 1;
        end
        if ~isfield(plotopt, 'plotWTLS')
            plotopt.plotWTLS = 0;
        end
    end

    if ~exist('hax', 'var')
        figure;
        set(gcf, 'Position', [520 118 650 680]);
        hax = gca;
        providedHax = 0;
    else
        providedHax = 1;
    end

    hax.PlotBoxAspectRatio = [1 1 1]; %[0.5 1.2 0.2294];

    hold on;
    % inferred mode from TAO data

    if plotopt.plotOLS
        if plotopt.ploterr
            handles.herro = ...
                supererr(hax, modes.InferredModeOLS{mm,nn}, ...
                         modes.depth{mm,nn} * -1, ...
                         abs(modes.InferredModeErrorOLS{mm,nn}), ...
                         [], 'I', capwidth, ...
                         'LineWidth', linewidth, 'Color', red, ...
                         'DisplayName', 'T_{TAO/TRITON} OLS');
            handles.herro = cut_nan(handles.herro(:,1));
        else
            handles.herro = plot(hax, modes.InferredModeOLS{mm,nn}, ...
                                 modes.depth{mm,nn} * -1, ...
                                 '.', 'Color', red, 'MarkerSize', 16, ...
                                 'DisplayName', 'T_{TAO/TRITON} OLS');
        end
    end

    if plotopt.plotWTLS
        if plotopt.ploterr
            handles.herrw = ...
                supererr(hax, modes.InferredModeWTLS{mm,nn}, ...
                         modes.depth{mm,nn} * -1, ...
                         abs(modes.InferredModeErrorWTLS{mm,nn}), ...
                         [], 'I', capwidth, ...
                         'LineWidth', linewidth, 'Color', green, ...
                         'DisplayName', 'T_{TAO/TRITON} WTLS');
            handles.herrw = cut_nan(handles.herrw(:,1));
        else
            handles.herrw = plot(hax, modes.InferredModeWTLS{mm,nn}, ...
                                 modes.depth{mm,nn} * -1, ...
                                 '.', 'Color', green, 'MarkerSize', 16, ...
                                 'DisplayName', 'T_{TAO/TRITON} WTLS');
        end
    end

    % 0 mean flow mode
    for ii=plotopt.nmode
        handles.hmode(ii) = ...
            plot(hax, squeeze(flatbot.IdealTempMode(mm,nn,:,ii)) ...
                 ./ max(abs(flatbot.IdealTempMode(mm,nn,:,ii))), ...
                 flatbot.zTmode * -1, 'LineWidth', linewidth, ...
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
            title(getTitleString(modes.lon(mm), modes.lat(nn)))
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

    if plotopt.plotOLS
        uistack(handles.herro, 'top');
    end
    if plotopt.plotWTLS
        uistack(handles.herrw, 'top');
    end
end
