% makes huge plot with a subplot for each point on the TAO array and plots
% theoretical and inferred mode

function [] = PlotModeMap(name)

  load([name '.mat']);

  hfig = figure; clf;
  set(hfig,'Position',[0 0 1600 900]);
  hash = githash([mfilename('fullpath') '.m']);
  insertAnnotation([name ': ' hash]);

  xlimits = [-1 1.2];
  ylimits = [0 600];
  fontSize = [20 24 28];

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
          subplot_index = sub2ind([nlon nlat],mm,nn)

          axes(hax(subplot_index));
          hold on

          if mod(subplot_index,nlon) == 1
              ylabel([num2str(modes.lat(nn))  'N']);
          else
              set(gca,'YTickLabel',[]);
          end
          if (subplot_index >= nlon*nlat - nlon+1)
              xlabel([num2str(modes.lon(mm)) 'W']);
          else
              set(gca,'XTickLabel',[]);
          end

          if isempty(modes.InferredMode{mm,nn})
              hax(subplot_index).YAxis.Color = [1 1 1];
              hax(subplot_index).XAxis.Color = [1 1 1];

              hax(subplot_index).XLabel.Color = hax(1).XLabel.Color;
              hax(subplot_index).YLabel.Color = hax(1).YLabel.Color;

              continue;
          end

          % for debugging subplot placement
          % text(0.5,0.5, num2str([mm nn]), 'Units', 'normalized');

          % temp std
          hstd = plot(data.Tstd{mm,nn}./max(data.Tstd{mm,nn}), ...
                      data.depth{mm,nn},'k');

          % 0 mean flow mode
          hideal = plot(squeeze(abs(modes.IdealTempMode(mm,nn,:))) ...
                        ./ max(abs( ...
                            modes.IdealTempMode(mm,nn,1:ind500))), ...
                        modes.zTmode);

          % inferred mode from TAO data
          hinfer = plot(modes.InferredMode{mm,nn}, ...
                        modes.depth{mm,nn}, 'o-');

          linex(0);
          xlim(xlimits); ylim(ylimits);
          set(gca,'YTick',[0 200 400]);
          revz;

          beautify(fontSize);
          set(gca,'TickLength',[0.03 0.03]);
          box off

          % set line width and ensure that beautify doesn't override
          linkprop([hstd hideal hinfer], 'LineWidth');
          hstd.LineWidth = 1;
          linkprop([hstd hideal hinfer], 'Tag');
          hstd.Tag = 'dcline';
      end
  end

  figure(hfig)
  [ax,h] = suplabel([name ' (' num2str(opt.filt.cutoff, '%.1f') ') | ' ...
                     opt.filt.window ' | ' ...
                     'Red = Theory, Blue = Inferred, ' ...
                     'Black = std(temp)'],'t');
  ax.YLabel.String = 'Z (m)';
  set(h,'FontSize',fontSize(3));

  linkaxes(ax, 'xy');

  %export_fig('-nocrop','-r150','../images/' name '.png');
end
