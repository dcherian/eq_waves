% makes huge plot with a subplot for each point on the TAO array and plots
% theoretical and inferred mode

function [] = PlotModeMap(plotopt)

  load([plotopt.name '.mat']);

  hfig = figure; clf;
  hfig.Position = [0 0 1600 900];
  hfig.Renderer = 'painters';

  hash = githash([mfilename('fullpath') '.m']);
  insertAnnotation([plotopt.name ': ' hash]);

  xlimits = [-0.2 1];
  ylimits = sort([0 600]*-1);
  fontSize = [20 24 28];
  linewidth = 1;
  labelcolor = [1 1 1]*0.3;

  nlon = length(modes.lon);
  nlat = length(modes.lat);

  % nlat rows x nlon columns
  hax = packfig(nlat,nlon);

  for mm=1:nlon
      for nn=1:nlat
          ind500 = find_approx(modes.zTmode,500,1);

          % to maximize tension, subplot counts along row first
          % while sub2ind does column first. So use sub2ind on
          % on transposed size matrix.
          subplot_index = sub2ind([nlon nlat],mm,nn);

          ax = hax(subplot_index);
          axes(ax);
          ax.Color = 'none';
          hold on

          if mod(subplot_index,nlon) == 1
              ylabel([num2str(modes.lat(nn))  'N']);
              ax.YLabel.Units = 'normalized';
              ax.YLabel.Position(1) = -0.8;
          else
              set(gca,'YTickLabel',[]);
          end

          if (subplot_index >= nlon*nlat - nlon+1)
              if mm ~= ceil(nlon/2)
                  ax.XAxis.Color = [1 1 1];
                  ax.XLabel.Color = labelcolor;
              end
              if modes.lon(mm) < 0
                  lonstr = 'E';
              else
                  lonstr = 'W';
              end
              xlabel([num2str(abs(modes.lon(mm))) lonstr]);
          else
              set(gca,'XTickLabel',[]);
              ax.XAxis.Color = [1 1 1];
          end

          if isempty(modes.InferredMode{mm,nn})
              ax.YAxis.Color = [1 1 1];
              ax.YLabel.Color = labelcolor;

              ax.XAxis.Color = [1 1 1];
              ax.XLabel.Color = labelcolor;
              continue;
          end

          % for debugging subplot placement
          % text(0.5,0.5, num2str([mm nn]), 'Units', 'normalized');

          % inferred mode from TAO data
          plot(modes.InferredMode{mm,nn}, ...
               modes.depth{mm,nn} * -1, '.', ...
               'MarkerSize', 18, ...
               'LineWidth', linewidth, 'DisplayName', 'T_{TAO/TRITON}');

          % 0 mean flow mode
          for ii=plotopt.nmode
              plot(squeeze(modes.IdealTempMode(mm,nn,:,ii)) ...
                   ./ max(abs(modes.IdealTempMode(mm,nn,1:ind500,ii))), ...
                   modes.zTmode * -1, 'LineWidth', linewidth, ...
                   'DisplayName', ['T_{bc' num2str(ii) '}']);
          end

          % temp std
          if plotopt.plotstd
              plot(data.Tstd{mm,nn}./max(data.Tstd{mm,nn}), ...
                   data.depth{mm,nn} * -1,'k', 'LineWidth', linewidth, ...
                   'DisplayName', 'T_{std}');
          end

          plot([0 0], ylimits, '--', 'Color', [1 1 1]*0.6, ...
               'LineWidth', 1, 'LegendDisplay', 'off');
          xlim(xlimits); ylim(ylimits);

          ax.YTick = ylimits(1):200:0;
          ax.YTickLabels{1} = '';
          ax.XTick = [0 1];

          if subplot_index == 1
              hleg = legend('Location', 'NorthWest');
              hleg.Box = 'off';
              hleg.Position(1) = hax(1).Position(1);
              hleg.Position(2) = 0.15;
          end

          beautify(fontSize);
      end
  end

  linkaxes(hax, 'xy');

  hax(nlon*(nlat-1)+ceil(nlon/2)).XAxis.Axle.VertexData(1) = 0;

  figure(hfig)
  [ax,~] = suplabel([opt.name ' | ' ...
                     opt.filt.window ' [' num2str(sort(opt.filt.cutoff), ...
                                                  '%.1f ') ']'], 't');
  ax.YLabel.String = 'Z (m)';
  ax.YLabel.Visible = 'on';
  ax.YLabel.Position(1) = -0.05;
  ax.Title.FontSize = fontSize(3);

  %export_fig('-nocrop','-r150','../images/' opt.name '.png');
end

function [] = beautify(fontSizes)

    font_name = 'Fira Sans';
    font_name_axis = font_name;

    % Get required handles for current figure
    hAxis = evalin('caller','gca');

    % Get some more handles
    hXLabel = hAxis.XLabel;
    hYLabel = hAxis.YLabel;
    hTitle  = hAxis.Title;

    set(hAxis, ...
        'Color', 'none', ...
        'FontName'    , font_name_axis, ...
        'Box'         , 'off'     , ...
        'TickDir'     , 'out'     , ...
        'TickLength'  , [1 1] * .06, ... % IMPROVE THIS
        'FontWeight'  , 'normal', ...
        ...%'XMinorTick'  , 'on'      , ...
        ...%'YMinorTick'  , 'on'      , ...
        ...%'ZMinorTick'  , 'on'      , ...
        'XGrid'       , 'off'      , ...
        'YGrid'       , 'off'      , ...
        'ZGrid'       , 'off'      , ...
        'FontSize'    , fontSizes(1), ...
        'LineWidth'   , 1        );

    set([hXLabel, hYLabel]  , ...
        'FontName'   , font_name, ...
        'FontWeight' , 'normal', ...
        'FontSize'   , fontSizes(2)         );

    set(hTitle, ...
        'FontSize'   , fontSizes(3) , ...
        'FontWeight' , 'normal', ...
        'FontName'   , font_name);
end
