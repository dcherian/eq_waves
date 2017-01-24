%  [handles] = PlotMode(modename, mm, nn, plotopt, hax)
function [handles] = PlotMode(modename, mm, nn, plotopt, hax)

    if ~exist('mm', 'var'), mm = 1; end
    if ~exist('nn', 'var'), nn = 1; end

    linewidth = 1;
    capwidth = 12;
    linestylemode = {'--'; '-'; '-.'}; % line style for theoretical modes.

    red = [217,95,2]/255; %[227,74,51]/255;
    blue = [117,112,179]/255;
    green = [27,158,119]/255; %[35,132,67]/255;

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
        [~, plotopt] = DefaultOptions;
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
        if ~isfield(plotopt, 'MarkWaterDepth')
            plotopt.MarkWaterDepth = 1;
        end
        if ~isfield(plotopt, 'plotPhaseLag')
            plotopt.plotPhaseLag = 1;
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

    %hax.PlotBoxAspectRatio = [1 1 1]; %[0.5 1.2 0.2294];
    axis fill;

    hold on;

    % inferred mode from TAO data
    if plotopt.plotOLS
        handles.hols = plot(hax, modes.InferredModeOLS{mm,nn}, ...
                            modes.depth{mm,nn} * -1, ...
                            '.', 'Color', red, 'MarkerSize', 16, ...
                            'DisplayName', 'T_{TAO/TRITON} OLS');
        if plotopt.ploterr
            handles.herrols = PlotErrorPatch(hax, ...
                                             modes.depth{mm,nn}, ...
                                             modes.InferredModeOLS{mm,nn}, ...
                                             modes.InferredModeErrorOLS{mm, nn},...
                                             red);
            % handles.herro = ...
            %     supererr(hax, modes.InferredModeOLS{mm,nn}, ...
            %              modes.depth{mm,nn} * -1, ...
            %              abs(modes.InferredModeErrorOLS{mm,nn}), ...
            %              [], 'I', capwidth, ...
            %              'LineWidth', linewidth, 'Color', red, ...
            %              'DisplayName', 'T_{TAO/TRITON} OLS');
            % handles.herro = cut_nan(handles.herro(:,1));
        end
    end

    if plotopt.plotPhaseLag
        phaselag = modes.phaselag{mm,nn};
        phasedepth = modes.depth{mm,nn}(~isnan(phaselag));
        phaselag = phaselag(~isnan(phaselag));
        handles.hlag = plot(hax, phaselag/180, -phasedepth, ...
                            '.-', 'Color', green, ...
                            'MarkerSize', 12, ...
                            'LineWidth', linewidth, ...
                            'DisplayName', 'Phase lag/180^o');
    end

    if plotopt.plotWTLS
        handles.hwtls = plot(hax, modes.InferredModeWTLS{mm,nn}, ...
                             modes.depth{mm,nn} * -1, ...
                             '.', 'Color', green, 'MarkerSize', 16, ...
                             'DisplayName', 'T_{TAO/TRITON} WTLS');

        if plotopt.ploterr
            handles.herrwtls = PlotErrorPatch(hax, ...
                                              modes.depth{mm,nn}, ...
                                              modes.InferredModeWTLS{mm,nn}, ...
                                              modes.InferredModeErrorWTLS{mm, nn},...
                                              green);
            % handles.herrw = ...
            %     supererr(hax, modes.InferredModeWTLS{mm,nn}, ...
            %              modes.depth{mm,nn} * -1, ...
            %              abs(modes.InferredModeErrorWTLS{mm,nn}), ...
            %              [], 'I', capwidth, ...
            %              'LineWidth', linewidth, 'Color', green, ...
            %              'DisplayName', 'T_{TAO/TRITON} WTLS');
            % handles.herrw = cut_nan(handles.herrw(:,1));
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
            plot(hax, modes.Tstd{mm,nn}./nanmax(modes.Tstd{mm,nn}), ...
                 modes.depth{mm,nn} * -1,'k.', 'LineWidth', linewidth, ...
                 'DisplayName', 'T_{std}');
    end

    % correlation coefficient
    if plotopt.plotcorr
        handles.hcorr = ...
            plot(hax, modes.corr{mm,nn}, modes.depth{mm,nn} * -1, ...
                 '.', 'MarkerSize', 12, 'Color', blue, ...
                 'DisplayName', 'Corr. Coeff.');
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

    if plotopt.MarkWaterDepth
        handles.hdepth = text(hax, 0.5, min(ylim), ...
                              [num2str(round(flatbot.etopoDepth(mm,nn))) ...
                            'm'], 'Units', 'data', ...
                              'Color', [1 1 1]*0.5, ...
                              'HorizontalAlignment', 'left', ...
                              'VerticalAlignment', 'bottom');
    end
    handles.h0 = plot(hax, [0 0], hax.YLim, '--', 'Color', [1 1 1]*0.6, ...
                      'LineWidth', 1, 'LegendDisplay', 'off');

    uistack(handles.h0, 'bottom');

    if plotopt.plotOLS
        uistack(handles.hols, 'top');
    end
    if plotopt.plotWTLS
        uistack(handles.hwtls, 'top');
    end
end

function [handle] = PlotErrorPatch(hax, depth, mode, error, color)

    depth = -abs(depth(~isnan(mode)));
    mode = cut_nan(mode);
    error = cut_nan(error);

    if ~isequal(size(depth), size(mode))
        depth = depth';
    end

    handle = patch(hax, ...
                   [mode-error; flip(mode+error); mode(1)-error(1)], ...
                   [depth; flip(depth); depth(1)], color, ...
                   'EdgeColor', 'None', 'FaceAlpha', 0.3, ...
                   'HandleVisibility', 'off');
end