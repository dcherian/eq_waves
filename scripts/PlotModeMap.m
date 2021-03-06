% makes huge plot with a subplot for each point on the TAO array and plots
% theoretical and inferred mode

function [hax, supax] = PlotModeMap(plotopt, lonrange, latrange, modes, opt)

  if ~exist('modes', 'var')
      load([plotopt.name '-' plotopt.window '.mat']);
  end

  if ~exist('lonrange', 'var')
      lonrange = 1:length(modes.lon);
      latrange = 1:length(modes.lat);
  end

  hfig = figure;
  hfig.Color = [1 1 1];
  hfig.Position = [0 0 1600 900];
  hfig.Renderer = 'painters';

  hash = githash([mfilename('fullpath') '.m']);
  insertAnnotation([plotopt.name ': ' hash]);

  xlimits = [-0.3 1.2];
  ylimits = [-750 0];
  fontSize = [20 24 28];
  linewidth = 1;
  labelcolor = 'k'; [1 1 1]*0.3;
  linestylemode = {'--'; '-'; '-.'}; % line style for theoretical modes.

  nlon = length(lonrange);
  nlat = length(latrange);

  % nlat rows x nlon columns
  hax = packfig(nlat,nlon);

  for mm=lonrange
      for nn=latrange
          % to maximize tension, subplot counts along row first
          % while sub2ind does column first. So use sub2ind on
          % on transposed size matrix.
          subplot_index = sub2ind([nlon nlat], ...
                                  mm-lonrange(1)+1, ...
                                  nn-latrange(1)+1);

          ax = hax(subplot_index);
          ax.Color = 'none';
          ax.FontSize = fontSize(1);
          ax.FontName = 'Times';
          ax.YLabel.FontName = 'Times';
          ax.XLabel.FontName = 'Times';
          ax.Box = 'off';
          ax.TickDir = 'out';
          ax.TickLength = [1 1]*0.06;
          ax.NextPlot = 'add'; % hold on

          if mod(subplot_index,nlon) == 1
              if modes.lat(nn) > 0
                  latstr = 'N';
              else
                  latstr = 'S';
              end
              ax.YLabel.String = [num2str(abs(modes.lat(nn)))  latstr];
              ax.YLabel.Units = 'normalized';
              ax.YLabel.Position(1) = -0.8;
          else
              ax.YTickLabel = [];
          end

          if (subplot_index >= nlon*nlat - nlon+1)
              if mm ~= ceil(nlon/2)
                  ax.XAxis.Color = [1 1 1];
                  ax.XLabel.Color = labelcolor;
              end
              if modes.lon(mm) > 0
                  lonstr = 'E';
              else
                  lonstr = 'W';
              end
              ax.XLabel.String = [num2str(abs(modes.lon(mm))) lonstr];
          else
              ax.XTickLabel = [];
              ax.XAxis.Color = [1 1 1];
          end

          if isempty(modes.InferredModeOLS{mm,nn})
              ax.YAxis.Color = [1 1 1];
              ax.YLabel.Color = labelcolor;

              ax.XAxis.Color = [1 1 1];
              ax.XLabel.Color = labelcolor;
              continue;
          end

          % for debugging subplot placement
          % text(0.5,0.5, num2str([mm nn]), 'Units', 'normalized');

          ax.XLim = xlimits;
          ax.YLim = ylimits;
          ax.YTick = ylimits(1):200:0;
          ax.YTickLabels{1} = '';
          ax.XTick = [0 1];

          handles = PlotMode(modes, mm, nn, plotopt, ax, 1);

          ax.XLabel.Color = labelcolor;
          ax.YLabel.Color = labelcolor;

          if plotopt.MarkWaterDepth
              handles.hdepth.FontName = 'Times';
              handles.hdepth.FontSize = fontSize(1)-4;
              handles.hdepth.Position(2) = ylimits(1);
          end

          if subplot_index == 1
              hleg = legend(ax, 'Location', 'NorthWest');
              hleg.Box = 'off';
              hleg.Position(1) = hax(1).Position(1);
              hleg.Position(2) = 0.15;
          end
      end
  end

  if strcmpi(opt.filt.window, 'butterworth')
      cutoff = opt.filt.cutoff/2;
  else
      cutoff = opt.filt.cutoff;
  end

  [supax,~] = suplabel([opt.name ' | ' ...
                      opt.filt.window ' [' num2str(sort(cutoff), ...
                                                   '%.1f ') ']'], 't');
  supax.Color = 'none';
  supax.YLabel.String = 'Z (m)';
  supax.YLabel.Visible = 'on';
  supax.YLabel.Position(1) = -0.05;
  supax.Title.FontSize = fontSize(3);

  supax.Title.FontName = 'Times';
  supax.YLabel.FontName = 'Times';
  supax.XLabel.FontName = 'Times';

  linkaxes(hax, 'xy');
  hax(nlon*(nlat-1)+ceil(nlon/2)).XAxis.Axle.VertexData(1,1) = 0;
  hax(nlon*(nlat-1)+ceil(nlon/2)).XAxis.Axle.VertexData(1,2) = 1;

  %export_fig('-nocrop','-r150','../images/' opt.name '.png');
end
