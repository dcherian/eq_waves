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

  hax = packfig(nlon,nlat);

  for mm=1:nlon
      for nn=1:nlat
          ind500 = find_approx(modes.zTmode,500,1);
          if isempty(modes.InferredMode{mm,nn}), continue; end
          subplot_index = sub2ind([nlon nlat],mm,nn);

          axes(hax(subplot_index));
          hold on

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

  % remove locations without data
  hax(4).delete;

  hax(46).YAxis.Visible = 'off';
  hax(46).XLabel.String = '140W';
  axes(hax(46)); beautify(fontSize);
  hax(46).XAxis.Color = [1 1 1];
  hax(46).XLabel.Color = hax(45).XLabel.Color;

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
